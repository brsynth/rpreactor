"""
Core code for firing rules
"""

import concurrent.futures
import logging
import sqlite3
import itertools
import csv
import collections

import pebble
from rdkit import Chem
from rdkit.Chem import AllChem

from rpreactor.chemical.standardizer import Standardizer
from rpreactor.rule.exceptions import ChemConversionError, RuleFireError, RuleConversionError


logger = logging.getLogger(__name__)


def _chunkify(it, size):
    it = iter(it)
    return iter(lambda: tuple(itertools.islice(it, size)), ())


class RuleBurner(object):
    """Apply any number of rules on any number of compounds.

    :param  with_hs:        bool, Enable explicit Hs when sanitizing chemicals
    :param  with_stereo:    bool, Keep stereochemistry (if any) when sanitizing chemicals
    """
    # TODO: rename class to 'RuleManager' or something?

    _SQL_CREATETABLE_MOLECULES = """
        CREATE TABLE IF NOT EXISTS molecules (
            id TEXT NOT NULL PRIMARY KEY,
            rd_mol BLOB NOT NULL,
            smiles TEXT NOT NULL,
            inchi TEXT NOT NULL,
            inchikey TEXT NOT NULL,
            is_computed INTEGER NOT NULL DEFAULT 0
        );
    """
    _SQL_CREATETABLE_RULES = """
        CREATE TABLE IF NOT EXISTS rules (
            id TEXT NOT NULL PRIMARY KEY,
            rd_rule BLOB NOT NULL,
            diameter INTEGER DEFAULT NULL
        );
    """
    _SQL_CREATETABLE_RESULTS = """
        CREATE TABLE IF NOT EXISTS results (
            sid TEXT NOT NULL,
            rid TEXT NOT NULL,
            pid TEXT,
            pgroup INTEGER,
            FOREIGN KEY (sid) REFERENCES molecules (id),
            FOREIGN KEY (rid) REFERENCES rules (id),
            FOREIGN KEY (pid) REFERENCES molecules (id)
        );
    """

    def __init__(self, database=None, with_hs=False, with_stereo=False):
        """Setting up everything needed for behavior decisions and firing rules."""
        # Set-up the database
        if database is None:
            database = ":memory:"
        self.db = sqlite3.connect(database)  # warning: only current thread will be able to use it
        self.db_path = database
        self.db.execute(self._SQL_CREATETABLE_MOLECULES)
        self.db.execute(self._SQL_CREATETABLE_RULES)
        self.db.execute(self._SQL_CREATETABLE_RESULTS)
        self.db.row_factory = sqlite3.Row
        self.db.commit()
        self._chemicals = None
        self._rules = None
        n_rules = self.db.execute("SELECT count(*) FROM rules").fetchone()['count(*)']
        n_mols = self.db.execute("SELECT count(*) FROM molecules").fetchone()['count(*)']
        logger.info(f"Connected to a database with {n_rules} rules and {n_mols} molecules (at '{database}').")
        # Sanitization
        self._with_hs = with_hs
        self._with_stereo = with_stereo
        if self._with_stereo:
            raise NotImplementedError("Stereo is not implemented at the time being.")

    @property
    def chemicals(self):
        """Returns a list of chemicals (molecules) identifier stored in the database."""
        if self._chemicals is None:
            self._chemicals = [x['id'] for x in self.db.execute("SELECT id FROM molecules").fetchall()]
        return self._chemicals

    @property
    def rules(self):
        """Returns a list of rules identifier stored in the database."""
        if self._rules is None:
            self._rules = [x['id'] for x in self.db.execute("SELECT id FROM rules").fetchall()]
        return self._rules

    @staticmethod
    def _task_fire(rd_rule, rd_mol, with_hs, with_stereo):
        """Apply one reaction rule on one chemical.

        We do not trust the task to terminate by itself, that's why we put it in a separate process that we can kill if
        need be.
        """
        try:
            ans = rd_rule.RunReactants((rd_mol,))
            rdmols, failed = RuleBurner._standardize_results(ans, with_hs, with_stereo)
            return RuleBurner._handle_results(rdmols)
        except Exception as e:
            raise RuleFireError(e) from e

    @staticmethod
    def _standardize_chemical(rdmol, with_hs, with_stereo, heavy=False):
        """Simple standardization of RDKit molecule."""
        params = {
            'OP_REMOVE_ISOTOPE': False,
            'OP_NEUTRALISE_CHARGE': False,
            'OP_REMOVE_STEREO': not with_stereo,
            'OP_COMMUTE_INCHI': True,
            'OP_KEEP_BIGGEST': False,
            'OP_ADD_HYDROGEN': with_hs,
            'OP_KEKULIZE': False,
            'OP_NEUTRALISE_CHARGE_LATE': True
        }
        if heavy:
            params['OP_REMOVE_ISOTOPE'] = True
            params['OP_NEUTRALISE_CHARGE'] = True
            params['OP_KEEP_BIGGEST'] = True
        return Standardizer(sequence_fun='sequence_tunable', params=params).compute(rdmol)

    @staticmethod
    def _standardize_results(tuple_tuple_rdmol, with_hs, with_stereo):
        """Perform sanitization and remove duplicates from reaction rule results.

        :param      tuple_tuple_rdmol:      tuple of tuple of RDKit Mol
        :returns    list_list_std:          list of list of standardized RDKit Mol
        :returns    list_idx_tuple_failed:  list of index of tuples that failed the standardization
        """
        uniq_depics = set()
        list_list_std = list()
        list_idx_tuple_failed = list()

        for idx_tuple, tuple_rdmol in enumerate(tuple_tuple_rdmol):
            try:
                list_std = list()
                list_inchikeys = list()
                # Standardize
                for rdmol in tuple_rdmol:
                    for rd_frag in Chem.GetMolFrags(rdmol, asMols=True, sanitizeFrags=False):
                        list_std.append(RuleBurner._standardize_chemical(rd_frag, with_hs, with_stereo))
                # Get InChIKey
                for rdmol in list_std:
                    inchikey = Chem.MolToInchiKey(rdmol)
                    if inchikey:
                        list_inchikeys.append(inchikey)
                    else:
                        msg = 'Product conversion to InChIKey raised an empty string'
                        logger.warning(ChemConversionError(msg))
                        raise ChemConversionError(msg)
                # Get unique depiction
                depic = '.'.join(sorted(list_inchikeys))
                # Store only if unique depiction never met
                if depic not in uniq_depics:
                    uniq_depics.add(depic)
                    list_list_std.append(list_std)
            except ChemConversionError as e:
                logger.warning("{}".format(e))
                list_idx_tuple_failed.append(idx_tuple)
                raise e
            except Exception as e:
                logger.warning("Cannot handle a tuple of result, skipped")
                logger.warning("{}".format(e))
                list_idx_tuple_failed.append(idx_tuple)

        return list_list_std, list_idx_tuple_failed

    @staticmethod
    def _handle_results(list_list_rdmol):
        """Generate InChIKey, InChI and SMILES from results.

        :param      list_list_rdmol:        list of list of RDKit Mol
        :returns    list_list_inchikeys:    list of list of InchiKeys
        :returns    list_list_inchis:       list of list of Inchis
        :returns    list_list_smiles:       list of list of SMILES
        """
        list_list_inchikeys = list()
        list_list_inchis = list()
        list_list_smiles = list()

        for list_rdmol in list_list_rdmol:
            try:
                list_inchikeys = list()
                list_inchis = list()
                list_smiles = list()
                for rdmol in list_rdmol:
                    # Get & check depictions
                    inchikey = Chem.MolToInchiKey(rdmol)  # TODO: this part could be optimized
                    inchi = Chem.MolToInchi(rdmol)
                    smiles = Chem.MolToSmiles(rdmol)
                    if not all([inchikey, inchi, smiles]):
                        raise ChemConversionError("Chemical conversion error")
                    # Store if we reach there
                    list_inchikeys.append(inchikey)
                    list_inchis.append(inchi)
                    list_smiles.append(smiles)
                # Store if we reach the end
                list_list_inchikeys.append(list_inchikeys)
                list_list_inchis.append(list_inchis)
                list_list_smiles.append(list_smiles)
            except ChemConversionError as e:
                logger.warning("{}".format(e))
                raise e
            except Exception as e:
                logger.warning("Cannot handle a tuple of result, skipped")
                logger.warning("{}".format(e))
        return list_list_rdmol, list_list_inchikeys, list_list_inchis, list_list_smiles  # Quick but dirty

    @staticmethod
    def _get_highest_int(untrusted_list):
        """Return the highest integer from a list of untrusted data that may or may not be cast to integers.

        Importantly, any value that can be cast to an integer will be considered. Default to 0.
        """
        highest = 0
        for item in untrusted_list:
            try:
                candidate = int(item)
                if candidate > highest:
                    highest = candidate
            except ValueError:
                pass  # silently ignore errors: those are probably not numbers
        return highest

    def _init_rdkit_rule(self, rsmarts):
        """Return RDKit reaction object."""
        try:
            rd_rule = AllChem.ReactionFromSmarts(rsmarts)
            rd_rule.Initialize()
        except Exception as e:
            raise RuleConversionError(e) from e
        return rd_rule

    def _init_rdkit_mol_from_inchi(self, inchi):
        """Return standardized RDKit molecule object."""
        try:
            rd_mol = Chem.MolFromInchi(inchi, sanitize=False)  # important: Sanitize = False
            rd_mol = RuleBurner._standardize_chemical(rd_mol, self._with_hs, self._with_stereo, heavy=True)
        except Exception as e:
            raise ChemConversionError(e) from e
        return rd_mol

    def _init_rdkit_mol_from_smiles(self, smiles):
        """Return standardized RDKit molecule object."""
        try:
            rd_mol = Chem.MolFromSmiles(smiles, sanitize=False)  # important: Sanitize = False
            rd_mol = RuleBurner._standardize_chemical(rd_mol, self._with_hs, self._with_stereo, heavy=True)
        except Exception as e:
            raise ChemConversionError(e) from e
        return rd_mol

    def _gen_rules(self, ids):
        """Generator of fully initialized/standardized RDKit reaction objects (rules)."""
        if ids is None:
            query = self.db.execute("select * from rules;")
        else:
            query = self.db.execute(f"select * from rules where rules.id in ({','.join(['?'] * len(ids))})", ids)
        for row in query:
            yield row['id'], Chem.rdChemReactions.ChemicalReaction(row['rd_rule'])

    def _gen_chemicals(self, ids):
        """Generator of fully initialized/standardized RDKit molecule objects (chemicals)."""
        if ids is None:
            query = self.db.execute("select * from molecules;")
        else:
            query = self.db.execute(f"select * from molecules where molecules.id in ({','.join(['?'] * len(ids))})", ids)
        for row in query:
            yield row['id'], Chem.Mol(row['rd_mol'])

    def _gen_couples(self, rule_mol):
        """Generator of fully initialized/standardized RDKit couples of rules and molecules."""
        # Make lists of all data that we need
        all_rid = set()
        all_cid = set()
        for rid, cid in rule_mol:
            all_rid.add(rid)
            all_cid.add(cid)
        logger.debug(rule_mol)
        # Gather everything (warning: watch out for memory limitation)
        rule_data = {k: v for k, v in self._gen_rules(list(all_rid))}
        mol_data = {k: v for k, v in self._gen_chemicals(list(all_cid))}
        # Finally, yield the data
        for rid, cid in rule_mol:
            try:
                yield rid, rule_data[rid], cid, mol_data[cid]
            except KeyError as err:
                logger.warning(f"It seems that data for {err} is not reachable. Is it a valid identifier?")

    def _gen_records(self, data, rdkit_func, blob_colname, other_colnames):
        """Helper generator of molecule or rule records for insertion in the database."""
        for key, value in data.items():
            # values of 'data' are either Dict-like or directly the input for rdkit_func
            if other_colnames:
                item = value[blob_colname]
                args = [value[x] for x in other_colnames]
            else:
                item = value
                args = []
            try:
                rd_item = rdkit_func(item)
                if blob_colname == "rd_mol":  # Special case: molecules always have smiles, inchi, inchikey
                    smiles = Chem.MolToSmiles(rd_item)
                    inchi = Chem.MolToInchi(rd_item)
                    inchikey = Chem.MolToInchiKey(rd_item)
                    args = [smiles, inchi, inchikey] + args
                record = (key, rd_item.ToBinary()) if not args else (key, rd_item.ToBinary(), *args)
                yield record
            except ChemConversionError as error:
                logger.error(f"Something went wrong converting chemical '{key}': {error}")
            except RuleConversionError as error:
                logger.error(f"Something went wrong converting rule '{key}': {error}")

    def _gen_precomputed_results(self, rule_mol):
        """Yields results found in the database."""
        query_results_str = """
        SELECT sid, rid, pid, pgroup, smiles, inchi, inchikey, rd_mol
        FROM results
        INNER JOIN molecules ON results.pid=molecules.id
        WHERE sid=? AND rid=?
        ORDER BY pgroup;
        """
        for rule, mol in rule_mol:
            result = None
            pgroup = None  # pgroup is a 0-based counter for each distinct solution
            for row in self.db.execute(query_results_str, [rule, mol]):
                if pgroup is None:             # init a new result object
                    result = {
                        'rule_id': rule,
                        'substrate_id': mol,
                        'product_list': [],
                        'product_inchikeys': [],
                        'product_inchis': [],
                        'product_smiles': [],
                    }
                if row['pgroup'] != pgroup:    # new group => a solution of one rule on one chemical
                    result['product_list'].append([Chem.Mol(row['rd_mol'])])
                    result['product_inchikeys'].append([row['inchikey']])
                    result['product_inchis'].append([row['inchi']])
                    result['product_smiles'].append([row['smiles']])
                elif row['pgroup'] == pgroup:  # metabolites sharing the same "pgroup" are coproducts of the same solution
                    result['product_list'][pgroup].append(Chem.Mol(row['rd_mol']))
                    result['product_inchikeys'][pgroup].append(row['inchikey'])
                    result['product_inchis'][pgroup].append(row['inchi'])
                    result['product_smiles'][pgroup].append(row['smiles'])
                pgroup = row['pgroup']  # remember last solution index
            if result:
                yield result

    def _gen_compute_results(self, rule_mol, commit, max_workers, timeout, chunk_size):
        """Yield new results using RDKit to apply a rule on a chemical."""
        with pebble.ProcessPool(max_workers=max_workers) as pool:
            # Prepare chunks of tasks
            # NB: it seems that pool.map does not avoid tasks to hold resources (memory) until they are consumed
            # even if a generator is used as input; so we use pool.schedule and we do our own chunks to avoid saturating
            # the RAM.
            logger.debug(f"Computing tasks in chunks of at most {chunk_size} couples (rule,  molecule) "
                         f"with {max_workers} workers and a per-task timeout of {timeout} seconds.")
            for chunk_idx, chunk in enumerate(_chunkify(rule_mol, chunk_size)):
                if chunk_idx > 0:
                    logger.debug(f"Working on task chunk #{chunk_idx+1}...")
                # Submit all the tasks for this chunk, and retrieve previous results
                all_running_tasks = []  # list of Future objects
                for rid, rd_rule, cid, rd_mol in self._gen_couples(chunk):
                    task = (rid, cid, pool.schedule(RuleBurner._task_fire,
                                                    args=(rd_rule, rd_mol, self._with_hs, self._with_stereo),
                                                    timeout=timeout))
                    all_running_tasks.append(task)
                # Gather the results
                for i, (rid, cid, future) in enumerate(all_running_tasks):
                    try:
                        rd_mol_list_list, inchikeys, inchis, smiles = future.result()
                        if rd_mol_list_list:  # silently discard tasks without a match
                            result = {
                                'rule_id': rid,
                                'substrate_id': cid,
                                'product_list': rd_mol_list_list,  # TODO: replace by list of ids?
                                'product_inchikeys': inchikeys,
                                'product_inchis': inchis,
                                'product_smiles': smiles,
                            }
                            if commit:
                                self._insert_result(rid, cid, rd_mol_list_list, inchikeys, inchis, smiles)
                            yield result
                    except concurrent.futures.TimeoutError:
                        logger.warning(f"Task {rid} on {cid} (#{i}) timed-out.")
                        # task['future'].cancel()  # NB: no need to cancel it, it's already canceled
                    except RuleFireError as error:
                        logger.error(f"Task {rid} on {cid} (#{i}) failed: {error}.")
                    except pebble.ProcessExpired as error:
                        logger.critical(f"Task {rid} on {cid} (#{i}) crashed unexpectedly: {error}.")
                # Attempt to free the memory
                del all_running_tasks

    def _insert_chemicals(self, data, rdkit_func):
        """Helper function to insert chemicals into the database."""
        # First, convert data to a Dict-like structure with key as identifiers
        if not isinstance(data, collections.Mapping):
            offset = RuleBurner._get_highest_int(self.chemicals) + 1 if self.chemicals else 0
            data = {k+offset: v for k, v in enumerate(data)}
        # Sniff the structure of data values (Dict-like or plain text)
        example = data[next(iter(data))]
        if isinstance(example, collections.Mapping):
            other_colnames = [x for x in example.keys() if x != 'rd_mol']
        else:
            other_colnames = []
        # Standardize items and generate the RDKit objects
        logger.debug(f"Inserting {len(data)} RDKit chemicals into the database.")
        cols_str = ','.join(['id', 'rd_mol', 'smiles', 'inchi', 'inchikey'] + other_colnames)
        values_str = ','.join(['?'] * (5 + len(other_colnames)))
        self.db.executemany(f"insert into molecules ({cols_str}) values ({values_str})",
                            self._gen_records(data, rdkit_func, 'rd_mol', other_colnames))
        self.db.commit()
        self._chemicals = None  # Reset the list of identifiers

    def _insert_rules(self, data, rdkit_func):
        """Helper function to insert rules into the database."""
        # First, convert data to a Dict-like structure with key as identifiers
        if not isinstance(data, collections.Mapping):
            offset = RuleBurner._get_highest_int(self.rules) + 1 if self.rules else 0
            data = {k+offset: v for k, v in enumerate(data)}
        # Sniff the structure of data values (Dict-like or plain text)
        example = data[next(iter(data))]
        if isinstance(example, collections.Mapping):
            other_colnames = [x for x in example.keys() if x != 'rd_rule']
        else:
            other_colnames = []
        # Standardize items and generate the RDKit objects
        logger.debug(f"Inserting {len(data)} RDKit rules into the database.")
        cols_str = ','.join(['id', 'rd_rule'] + other_colnames)
        values_str = ','.join(['?'] * (2 + len(other_colnames)))
        self.db.executemany(f"insert into rules ({cols_str}) values ({values_str})",
                            self._gen_records(data, rdkit_func, 'rd_rule', other_colnames))
        self.db.commit()
        self._rules = None  # Reset the list of identifiers

    def _insert_result(self, rid, cid, rd_mol_list_list, inchikeys, inchis, smiles):
        """Helper to insert results and newly found molecules."""
        # First, commit all the chemicals if they are not known
        next_valid_id = RuleBurner._get_highest_int(self.chemicals) + 1
        for idx1, rd_mol_list in enumerate(rd_mol_list_list):
            for idx2, rd_mol in enumerate(rd_mol_list):
                # Is this chemical known?
                this_inchi = inchis[idx1][idx2]
                ans = self.db.execute("select id from molecules where inchi=? limit 1;",
                                      [this_inchi]).fetchone()
                if ans:
                    this_id = ans['id']  # arbitrarily keep the first occurence
                    chemical_need_commit = False
                else:
                    this_id = next_valid_id
                    chemical_need_commit = True
                    next_valid_id += 1
                # Insert what needs to be inserted
                if chemical_need_commit:
                    record = (this_id, rd_mol.ToBinary(), smiles[idx1][idx2], this_inchi,
                              inchikeys[idx1][idx2], 1)
                    self.db.execute("insert into molecules values (?,?,?,?,?,?);", record)
                self.db.execute("insert into results values (?,?,?,?);", (rid, cid, this_id, idx1))
        self.db.commit()

    def insert_rsmarts(self, data):
        """Insert reaction rules defined as reaction SMARTS into the database."""
        self._insert_rules(data, self._init_rdkit_rule)

    def insert_inchi(self, data):
        """Insert molecules defined as InChI into the database."""
        self._insert_chemicals(data, self._init_rdkit_mol_from_inchi)

    def insert_smiles(self, data):
        """Insert molecules defined as InChI into the database."""
        self._insert_chemicals(data, self._init_rdkit_mol_from_smiles)

    def create_indexes(self):
        """Create SQL indexes on the database."""
        self.db.execute("CREATE INDEX IF NOT EXISTS idx_molecules ON molecules(id);")
        self.db.execute("CREATE INDEX IF NOT EXISTS idx_rules ON rules(id);")
        self.db.commit()

    def dump_to_sql(self, path):
        """Dump the database as a SQL file."""
        with open(path, 'w') as f:
            for line in self.db.iterdump():
                f.write(f'{line}\n')

    def drop_results(self):
        """Drop the results table and computed metabolites."""
        self.db.execute("DELETE FROM results;")
        self.db.execute("DELETE FROM molecules WHERE is_computed=1;")
        self.db.commit()

    def compute(self, rule_mol=None, commit=False, max_workers=1, timeout=60, chunk_size=1000):
        """Apply all rules on all chemicals and returns a generator over the results."""
        # "None" means all rules against all mol
        if rule_mol is None:
            rule_mol = [(rule, mol) for rule in self.rules for mol in self.chemicals]
        # First of all, yield precomputed results
        for result in self._gen_precomputed_results(rule_mol):
            rule_mol.remove((result['rule_id'], result['substrate_id']))  # TODO: check with list of list
            yield result
        # Only then, continue with non pre-computed results
        for result in self._gen_compute_results(rule_mol, commit, max_workers, timeout, chunk_size):
            yield result


def _create_db_from_retrorules_v1_0_5(path_retrosmarts_tsv, path_sqlite, with_hs, with_stereo):
    def helper_metabolite(store, cid, smiles):
        if cid not in store:
            store[cid] = smiles
        elif store[cid] != smiles:
            logger.warning(f"Metabolite {cid} is suspiciously associated to distinct SMILES. "
                           f"Only the first one will be considered: {store[cid]} and {smiles}")
    rules = {}
    metabolites = {}
    # Load all valuable data in-memory
    # NB: is this file small enough that we do not need to chunk it?
    with open(path_retrosmarts_tsv) as hdl:
        for row in csv.DictReader(hdl, delimiter='\t'):
            # each row contains 1 rule...
            # NB: rule identifier may be duplicated over several rows but must match the same RSMARTS
            rid = row["# Rule_ID"]
            rsmarts = row["Rule_SMARTS"]
            diameter = int(row["Diameter"])
            # TODO add Reaction_direction
            if rid not in rules:
                rules[rid] = {'rd_rule': rsmarts, 'diameter': diameter}  # must match database schema
            else:
                assert rules[rid]['rd_rule'] == rsmarts, f"UNEXPECTED: rule {rid} from {path_retrosmarts_tsv} has " \
                                                         f"mismatching RSMARTS: {rules[rid]} and {rsmarts}"
            # ... and 1 substrate ...
            helper_metabolite(metabolites, row["Substrate_ID"], row["Substrate_SMILES"])
            # ... and N products
            smiles_list = row["Product_SMILES"].split('.')
            for idx, cid in enumerate(row["Product_IDs"].split('.')):
                helper_metabolite(metabolites, cid, smiles_list[idx])
    # Create the database
    o = RuleBurner(database=path_sqlite, with_hs=with_hs, with_stereo=with_stereo)
    o.insert_rsmarts(rules)
    o.insert_smiles(metabolites)
    o.create_indexes()


def create_db_from_retrorules(path_retrosmarts_tsv, path_sqlite, with_hs=False, with_stereo=False, version="v1.0"):
    """Convert a RetroRules dataset to a rpreactor-ready sqlite3 database.

    All rules and all molecules will be extracted from the TSV file and imported into the database.
    For more information on RetroRules, see https://retrorules.org/.
    """
    if version.startswith("v1.0"):
        _create_db_from_retrorules_v1_0_5(path_retrosmarts_tsv, path_sqlite, with_hs, with_stereo)
