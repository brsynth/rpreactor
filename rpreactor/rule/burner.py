"""
Main module to handle reaction rules and to apply them on metabolites.
"""

import concurrent.futures
import logging
import sqlite3
import itertools
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
    """This class manage the connection to a SQLite database that contains standardized molecules and reaction rules.

    :param  db_path:        Path to the SQLite database (will be created if it doesn't exist). If None is provided
                            an in-memory database will be created. Default: None.
    :param  with_hs:        Enable explicit hydrogens when sanitizing chemicals. Default: False.
    :param  with_stereo:    Keep stereochemistry (if any) when sanitizing chemicals. Default: False.
    :type db_path:          str
    :type with_hs:          bool
    :type with_stereo:      bool
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
            diameter INTEGER DEFAULT NULL,
            usage INTEGER DEFAULT 0 CHECK (usage in (-1, 0, 1))
        );
    """
    _SQL_CREATETABLE_RESULTS = """
        CREATE TABLE IF NOT EXISTS results (
            sid TEXT NOT NULL,
            rid TEXT NOT NULL,
            pid TEXT,
            pstoichio INTEGER,
            pgroup INTEGER,
            PRIMARY KEY (sid, rid, pid, pgroup),
            FOREIGN KEY (sid) REFERENCES molecules (id),
            FOREIGN KEY (rid) REFERENCES rules (id),
            FOREIGN KEY (pid) REFERENCES molecules (id)
        );
    """
    _SQL_CREATETABLE_ATTEMPTS = """
        CREATE TABLE IF NOT EXISTS batch_attempts (
            sid TEXT NOT NULL,
            diameter INTEGER,
            usage INTEGER NOT NULL,
            PRIMARY KEY (sid, diameter, usage),
            FOREIGN KEY (sid) REFERENCES molecules (id)
        );
    """
    _SQL_CREATETABLE_CONFIG = """
    CREATE TABLE IF NOT EXISTS _config (
        key TEXT NOT NULL,
        value TEXT
    );
    """

    _CODE_RULEUSAGE = {'FORWARD': 1, 'RETRO': -1, 'BOTH': 0}

    _DEFAULT_HS = True
    _DEFAULT_STEREO = False

    def __init__(self, db_path=None, with_hs=None, with_stereo=None):
        """Setting up everything needed for behavior decisions and firing rules.

        The database will be created if it does not already exists. Parameters aside 'db_path' are only required
        if the database needs to be created.
        """
        # Set-up the database
        if db_path is None:
            db_path = ":memory:"
        self.db = sqlite3.connect(db_path)  # warning: only current thread will be able to use it
        self.db_path = db_path
        self.db.execute(self._SQL_CREATETABLE_MOLECULES)
        self.db.execute(self._SQL_CREATETABLE_RULES)
        self.db.execute(self._SQL_CREATETABLE_RESULTS)
        self.db.execute(self._SQL_CREATETABLE_CONFIG)
        self.db.execute(self._SQL_CREATETABLE_ATTEMPTS)
        self.db.row_factory = sqlite3.Row
        self.db.commit()
        # Check the database configuration
        db_config = self._read_database_config()
        input_state = {'with_hs': with_hs, 'with_stereo': with_stereo}
        if db_config:  # database already configured: we must check that everything is coherent with user input
            for k, v in input_state.items():
                if v is not None and db_config[k] != v:
                    raise ValueError(f"Impossible to initialize the database with parameter '{k}={v}': "
                                     f"database was already configured with '{k}={db_config[k]}'!")
        else:         # this is a new database: the configuration needs to be written
            self._set_database_config(with_hs, with_stereo)
            db_config = self._read_database_config()
        # Synchronize object and database configuration
        self._with_hs = db_config['with_hs']
        self._with_stereo = db_config['with_stereo']
        if self._with_stereo:
            raise NotImplementedError("Stereo is not implemented at the time being.")
        # Set up private attributes
        self._chemicals = None
        self._rules = None
        self._precomputed_count = 0
        self._newlycomputed_count = 0
        self._timeout_list = []
        self._errors_list = []

    def __str__(self):
        """String representation of the summary."""
        s = self.summary()
        txt = f"A RuleBurner object connected to a database with {s['database_rules_count']} rules, " \
              f"{s['database_chemical_count']} compounds, and {s['database_results_count']} results " \
              f"(at '{self.db_path}'). " \
              f"Rules {'use' if self._with_stereo else 'ignore'} stereo information, " \
              f"and {'use' if self._with_hs else 'ignore'} hydrogens. " \
              f"Last compute call yield {s['lastcompute_precomputed_count']+s['lastcompute_newlycomputed_count']} " \
              f"results ({s['lastcompute_precomputed_count']} precomputed, {s['lastcompute_newlycomputed_count']} new); " \
              f"{len(s['lastcompute_errors_list'])} errors were caught, " \
              f"and {len(s['lastcompute_timeout_list'])} timeouts were hit. More details are available in the logs."
        return txt

    @property
    def chemicals(self):
        """Returns a list of chemicals (molecules) identifiers currently used in the database."""
        if self._chemicals is None:
            self._chemicals = [x['id'] for x in self.db.execute("SELECT id FROM molecules").fetchall()]
        return self._chemicals

    @property
    def rules(self):
        """Returns a list of rules identifiers currently used in the database."""
        if self._rules is None:
            self._rules = [x['id'] for x in self.db.execute("SELECT id FROM rules").fetchall()]
        return self._rules

    def _read_database_config(self):
        """Read and return database configuration table."""
        state = {row[0]: row[1] for row in self.db.execute("SELECT key, value FROM _config").fetchall()}
        try:
            state['with_hs'] = bool(int(state['with_hs']))
            state['with_stereo'] = bool(int(state['with_stereo']))
        except KeyError:
            pass  # not configured yet
        return state

    def _set_database_config(self, with_hs, with_stereo):
        if with_hs is None:
            with_hs = self._DEFAULT_HS
            logger.warning(f"Parameter 'with_hs' was not provided: defaulting to 'with_hs={with_hs}'. "
                           "This may change in the future.")
        if with_stereo is None:
            with_stereo = self._DEFAULT_STEREO
            logger.warning(f"Parameter 'with_stereo' was not provided: defaulting to 'with_stereo={with_stereo}'. "
                           "This may change in the future.")
        self.db.executemany("INSERT INTO _config VALUES (?,?)", [('with_hs', with_hs), ('with_stereo', with_stereo)])
        self.db.commit()

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
            except Exception as e:
                logger.warning(f"Some rule product(s) could not be standardized and will be ignored. "
                               f"Caught a {type(e).__name__} error: {e}")
                list_idx_tuple_failed.append(idx_tuple)
        return list_list_std, list_idx_tuple_failed

    @staticmethod
    def _handle_results(list_list_rdmol):
        """Generate InChIKey, InChI and SMILES from results.

        :param      list_list_rdmol:        list of list of RDKit Mol
        :returns    list_list_rdmol:        list of list of RDKit Mol successfully converted
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
                    inchikey = Chem.MolToInchiKey(rdmol)
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
            except Exception as e:
                logger.warning(f"Some rule product(s) could not be represented as InChI, InChIKey or SMILES "
                               f"and will be ignored. Caught a {type(e).__name__} error: {e}")
        return list_list_rdmol, list_list_inchikeys, list_list_inchis, list_list_smiles  # Quick but dirty

    @staticmethod
    def _get_highest_int(untrusted_list):
        """Return the highest integer from a list of untrusted data that may or may not be cast to integers.

        Importantly, any value that can be cast to an integer will be considered. Default to 0.
        """
        if isinstance(untrusted_list, str):  # Prevent miss use of subscriptable
            raise TypeError(f'Type str is not a valid subscriptable type')
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
            msg = f"'{inchi}' is not a valid InChI."
            raise ChemConversionError(msg) from e
        return rd_mol

    def _init_rdkit_mol_from_smiles(self, smiles):
        """Return standardized RDKit molecule object."""
        try:
            rd_mol = Chem.MolFromSmiles(smiles, sanitize=False)  # important: Sanitize = False
            rd_mol = RuleBurner._standardize_chemical(rd_mol, self._with_hs, self._with_stereo, heavy=True)
        except Exception as e:
            msg = f"'{smiles}' is not a valid SMILES."
            raise ChemConversionError(msg) from e
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
                raise ValueError(f"{err} is not a valid identifier in task (rid={rid}, cid={cid}).")

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
        SELECT sid, rid, pid, pstoichio, pgroup, smiles, inchi, inchikey, rd_mol
        FROM results
        INNER JOIN molecules ON results.pid=molecules.id
        WHERE rid=? AND sid=?
        ORDER BY pgroup;
        """
        for rule, mol in rule_mol:
            result = None
            pgroup = None  # pgroup is a 0-based counter for each distinct solution
            for row in self.db.execute(query_results_str, [rule, mol]):
                stoichio = row['pstoichio']
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
                    result['product_list'].append([Chem.Mol(row['rd_mol'])] * stoichio)
                    result['product_inchikeys'].append([row['inchikey']] * stoichio)
                    result['product_inchis'].append([row['inchi']] * stoichio)
                    result['product_smiles'].append([row['smiles']] * stoichio)
                elif row['pgroup'] == pgroup:  # metabolites sharing the same "pgroup" are coproducts of the same solution
                    result['product_list'][-1] + [Chem.Mol(row['rd_mol'])] * stoichio
                    result['product_inchikeys'][-1] + [row['inchikey']] * stoichio
                    result['product_inchis'][-1] + [row['inchi']] * stoichio
                    result['product_smiles'][-1] + [row['smiles']] * stoichio
                pgroup = row['pgroup']  # remember last solution index
            if result is not None:
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
                # Submit all the tasks for this chunk
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
                        self._timeout_list.append((rid, cid))
                    except RuleFireError as error:
                        logger.error(f"Task {rid} on {cid} (#{i}) failed: {error}.")
                        self._errors_list.append((rid, cid))
                    except pebble.ProcessExpired as error:
                        logger.critical(f"Task {rid} on {cid} (#{i}) crashed unexpectedly: {error}.")
                        self._errors_list.append((rid, cid))
                # Attempt to free the memory
                del all_running_tasks

    def _insert_chemicals(self, data, rdkit_func):
        """Helper function to insert chemicals into the database."""
        n_before = self.db.execute("SELECT count(*) FROM molecules").fetchone()['count(*)']
        # First, convert data to a Dict-like structure with key as identifiers
        if not isinstance(data, collections.abc.Mapping):
            offset = RuleBurner._get_highest_int(self.chemicals) + 1 if self.chemicals else 0
            data = {k+offset: v for k, v in enumerate(data)}
        # Sniff the structure of data values (Dict-like or plain text)
        example = data[next(iter(data))]
        if isinstance(example, collections.abc.Mapping):
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
        n_after = self.db.execute("SELECT count(*) FROM molecules").fetchone()['count(*)']
        return n_after - n_before

    def _insert_rules(self, data, rdkit_func):
        """Helper function to insert rules into the database."""
        n_before = self.db.execute("SELECT count(*) FROM rules").fetchone()['count(*)']
        # First, convert data to a Dict-like structure with key as identifiers
        if not isinstance(data, collections.abc.Mapping):
            offset = RuleBurner._get_highest_int(self.rules) + 1 if self.rules else 0
            data = {k+offset: v for k, v in enumerate(data)}
        # Sniff the structure of data values (Dict-like or plain text)
        example = data[next(iter(data))]
        if isinstance(example, collections.abc.Mapping):
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
        n_after = self.db.execute("SELECT count(*) FROM rules").fetchone()['count(*)']
        return n_after - n_before

    def _insert_result(self, rid, cid, rd_mol_list_list, inchikeys, inchis, smiles):
        """Helper to insert results and newly found molecules."""
        # First, commit all the chemicals if they are not known
        next_valid_id = RuleBurner._get_highest_int(self.chemicals) + 1
        for idx_solution, rd_mol_list in enumerate(rd_mol_list_list):
            stoichio = {}      # <product chemical id>: <stoichiometry>
            seen_inchi = {}    # <product inchi>: <product chemical id>
            insert_results = set()
            for idx_coproduct, rd_mol in enumerate(rd_mol_list):
                # Is this product a known chemical in the database?
                this_inchi = inchis[idx_solution][idx_coproduct]
                if this_inchi in seen_inchi:  # product is duplicated: its stoichio is >1
                    stoichio[seen_inchi[this_inchi]] += 1
                    continue
                ans = self.db.execute("select id from molecules where inchi=? limit 1;", [this_inchi]).fetchone()
                if ans:  # product is already known in the database
                    this_id = ans['id']  # arbitrarily keep the first occurence
                    chemical_need_commit = False
                else:    # product is new
                    this_id = next_valid_id
                    chemical_need_commit = True
                    next_valid_id += 1
                # Remember this product in case its stoichio is >1
                seen_inchi[this_inchi] = this_id
                stoichio[this_id] = 1
                # Insert what needs to be inserted
                if chemical_need_commit:
                    record = (this_id, rd_mol.ToBinary(), smiles[idx_solution][idx_coproduct], this_inchi,
                              inchikeys[idx_solution][idx_coproduct], 1)
                    self.db.execute("insert into molecules values (?,?,?,?,?,?);", record)
                    self._chemicals.append(str(this_id))  # manual update of the list of used identifiers (kinda dangerous)
                insert_results.add((cid, rid, this_id, idx_solution))
            # Add the stoichiometry before to commit the results
            tmp = [(sid, rid, pid, stoichio[pid], pgroup) for sid, rid, pid, pgroup in insert_results]
            self.db.executemany("insert into results values (?,?,?,?,?);", tmp)
        self.db.commit()

    def insert_rsmarts(self, data):
        """Insert reaction rules defined as reaction SMARTS into the database.

        :param data: A Dict-like structure where keys are new identifiers and values SMARTS strings.
        :return: The number of inserted records.
        """
        return self._insert_rules(data, self._init_rdkit_rule)

    def insert_inchi(self, data):
        """Insert molecules defined as InChI into the database.

        :param data: A Dict-like structure where keys are new identifiers and values InChI strings.
        :return: The number of inserted records.
        """
        return self._insert_chemicals(data, self._init_rdkit_mol_from_inchi)

    def insert_smiles(self, data):
        """Insert molecules defined as SMILES into the database.

        :param data: A Dict-like structure where keys are new identifiers and values SMILES strings.
        :return: The number of inserted records.
        """
        return self._insert_chemicals(data, self._init_rdkit_mol_from_smiles)

    def create_indexes(self):
        """Create SQL indexes on the database.

        * Table `molecules`: `id`, `inchi` and `inchikey`.
        * Table `rules`: `id` and `(diameter, usage)`.
        * Table `results`: `(rid, sid)`.

        Important: Keep in mind that to benefit from a multi-index, the WHERE clause in an SQL queries must use
        the same column order as in the index.

        Warning: this will slow down inserts but should fasten computations.
        """
        # NB: index are auto-created for primary keys
        self.db.execute("CREATE INDEX IF NOT EXISTS idx_molecules_inchi ON molecules(inchi);")
        self.db.execute("CREATE INDEX IF NOT EXISTS idx_molecules_inchikey ON molecules(inchikey);")
        self.db.execute("CREATE INDEX IF NOT EXISTS idx_rules_diausa ON rules(diameter, usage);")
        self.db.execute("CREATE INDEX IF NOT EXISTS idx_results ON results(rid, sid);")
        self.db.commit()

    def dump_to_sql(self, path):
        """Dump the database as a SQL file.

        :param path: Path to the SQL file.
        :type path: str
        """
        with open(path, 'w') as f:
            for line in self.db.iterdump():
                f.write(f'{line}\n')

    def drop_results(self):
        """Drop the results table and computed metabolites."""
        self.db.execute("DELETE FROM batch_attempts;")
        self.db.execute("DELETE FROM results;")
        self.db.execute("DELETE FROM molecules WHERE is_computed=1;")
        self.db.commit()

    def get_chemical(self, mid):
        """Get a copy of chemical <mid> as a RDKit molecule object."""
        try:
            ans = Chem.Mol(self.db.execute("select rd_mol from molecules where id=?;", [mid]).fetchone()[0])
        except TypeError as err:
            raise ValueError(f"{mid} is not in the database.") from err
        return ans

    def get_rule(self, rid):
        """Get a copy of rule <rid> as a RDKit ChemReaction object."""
        try:
            tmp = self.db.execute("select rd_rule from rules where id=?;", [rid]).fetchone()
            ans = Chem.rdChemReactions.ChemicalReaction(tmp[0])
        except TypeError as err:
            raise ValueError(f"{rid} is not in the database.") from err
        return ans

    def list_chemicals_with_inchi(self, inchi):
        """List chemicals ids having this InChI."""
        mids = [x[0] for x in self.db.execute("select id from molecules where inchi=?;", [inchi])]
        return mids

    def list_chemicals_with_inchikey(self, inchikey):
        """List chemicals ids having this InChIKey."""
        mids = [x[0] for x in self.db.execute("select id from molecules where inchikey=?;", [inchikey])]
        return mids

    def list_rules(self, diameter, usage):
        """List rule ids at this diameter and usage.

        :param diameter: Rule diameter(s) to select.
        :type diameter: int or list of int
        :param usage: A list of rule usage to select among: 'FORWARD', 'RETRO', or 'BOTH'.
        :type diameter: list of str
        """
        # Checks and casts
        try:
            if any(x not in self._CODE_RULEUSAGE for x in usage):
                raise ValueError("Elements of 'usage' should be: 'FORWARD', 'RETRO', or 'BOTH'.")
        except TypeError as err:
            raise ValueError("'usage' should be a list.") from err
        if isinstance(diameter, int):
            diameter = [diameter]
        # Query
        usage_codes = [self._CODE_RULEUSAGE[x] for x in usage]
        query = f"""
        select id
        from rules
        where diameter in ({','.join(['?']*len(diameter))}) and usage in ({','.join(['?']*len(usage_codes))})
        """
        rule_list = [x[0] for x in self.db.execute(query, [*diameter, *usage_codes]).fetchall()]
        return rule_list

    def compute_batch_rules(self, mols, diameter, usage, commit=False, **kwargs):
        """Computes all rules at `diameter` and `usage` against each chemical identifier referenced in `mols`.

        In addition of collecting the right set of rules given `diameter` and `usage`, this method may commit to the
        database current query attempt. This has the advantage to allow for faster results retrieval if the same query
        was to be attempted again. Indeed, `compute()` may store postive results but negative results are not stored
        for performance reasons. Knowing which queries were already attempted avoid recomputing (negative) results.

        :param mols: List of chemical identifiers.
        :param diameter: Rule diameter(s) to select.
        :param usage: A list of rule usage to select among: 'FORWARD', 'RETRO', or 'BOTH'.
        :param commit: If true, the results will be commited to the database for later (faster) retrieval. Default: False.
        :param kwargs: Other parameters passed to the compute() call.
        :return: a generator over products predicted by applying all rules on molecules
        """
        already_attempted = {}  # already_attempted[(diameter,usage)] = [mol1, mol2, ...]
        # Checks and quick casts
        try:
            usage_codes = [self._CODE_RULEUSAGE[x] for x in usage]
        except KeyError as err:
            raise ValueError("Elements of 'usage' should be: 'FORWARD', 'RETRO', or 'BOTH'.") from err
        if isinstance(diameter, int):
            diameter = [diameter]
        # Check for previous attempts (not all of them have results)
        ph_mols = ','.join(['?'] * len(mols))
        ph_diameter = ','.join(['?'] * len(diameter))
        ph_usage = ','.join(['?'] * len(usage))
        query = f"""
        select sid, diameter, usage from batch_attempts
        where sid in ({ph_mols}) and diameter in ({ph_diameter}) and usage in ({ph_usage})
        """
        for row in self.db.execute(query, mols + diameter + usage_codes):
            key = (row['diameter'], row['usage'])  # usage code
            if key not in already_attempted:
                already_attempted[key] = []
            already_attempted[key].append(row['sid'])
        # Collect all tasks having pregenerated results
        # TODO: directly yield the results here?
        query = f"""
        select distinct sid, rid from results
        inner join rules on results.rid=rules.id
        where sid in ({ph_mols}) and diameter in ({ph_diameter}) and usage in ({ph_usage})
        """
        rule_mol = set([(x['rid'], x['sid']) for x in self.db.execute(query, mols + diameter + usage_codes)])
        # Add to the list of tasks everything that has not been attempted yet...
        new_attempts = []
        for this_diameter in diameter:
            for this_usage in usage:
                this_usage_code = self._CODE_RULEUSAGE[this_usage]
                relevant_rules = self.list_rules(this_diameter, [this_usage])
                try:
                    skip_sid = already_attempted[(this_diameter, this_usage_code)]
                except KeyError:
                    skip_sid = []
                for sid in mols:
                    if sid not in skip_sid:
                        new_attempts.append((sid, this_diameter, this_usage_code))
                        rule_mol |= set([(rule, sid) for rule in relevant_rules])
        # ... and eventually remember current attempts
        if commit and new_attempts:
            self.db.executemany("insert into batch_attempts values (?,?,?);", new_attempts)
            self.db.commit()
        return self.compute(rule_mol=list(rule_mol), commit=commit, **kwargs)

    def compute(self, rule_mol, commit=False, max_workers=1, timeout=60, chunk_size=1000):
        """Returns a generator over products predicted by applying rules on molecules.

        Important: rules and molecules must be inserted into the database beforehand and are referred to by their
        identifier. We assume each identifier refer to distinct rules and molecules: several identifiers referring to
        identical rules/molecules will be computed independently.

        :param rule_mol: Tasks to compute as in: [(<rule_id>, <mol_id>), ...]. If the magic argument "*" is provided,
            will compute all rules against all metabolites (mainly useful for debug purpose).
        :type rule_mol: list of tuple, or '*'
        :param commit: If true, the results be commited to the `results` table for later retrieval. Default: False.
        :type commit: bool
        :param max_workers: Maximum number of cores to use simultaneously to compute the tasks. Default: 1.
        :type max_workers: int
        :param timeout: Maximum time allowed for a single rule to apply on a single molecule (in seconds). Default: 60.
        :type timeout: float
        :param chunk_size: The maximum chunk size of tasks to hold in-memory at once. Increasing this number will yield
            better performacne but may saturate the memory. Default: 1000.
        :type chunk_size: int
        :return result: Dictionary (yield) providing result from one rule application
        """
        # Reset the summary attributes
        self._precomputed_count = 0
        self._newlycomputed_count = 0
        self._timeout_list = []
        self._errors_list = []
        # Deal with magic arguments '*' that means "all"
        if rule_mol == '*':
            rule_mol = [(rule, mol) for rule in self.rules for mol in self.chemicals]
        # First of all, yield precomputed results
        already_computed_tasks = []
        for result in self._gen_precomputed_results(rule_mol):
            already_computed_tasks.append((result['rule_id'], result['substrate_id']))
            self._precomputed_count += 1
            yield result
        # Remove those already-yield results from the query
        # NB: pregenerated results from the rules dataset are considered precomputed although there was no promiscuity
        for task in already_computed_tasks:
            try:
                rule_mol.remove(task)
            except ValueError:  # (rule_id, substrate_id) was not found... could it be a list?
                rule_mol.remove(list(task))
        # Only then, continue with non pre-computed results
        for result in self._gen_compute_results(rule_mol, commit, max_workers, timeout, chunk_size):
            self._newlycomputed_count += 1
            yield result

    def summary(self):
        """Returns a summary of the database state and last compute execution.

        :return: Dict
        """
        ans = {
            "database_path": self.db_path,
            "database_with_hs": self._with_hs,
            "database_with_stereo": self._with_stereo,
            "database_rules_count": len(self.rules),
            "database_chemical_count": len(self.chemicals),
            "database_results_count": self.db.execute("select count(*) from results").fetchone()[0],
            "lastcompute_precomputed_count": self._precomputed_count,
            "lastcompute_newlycomputed_count": self._newlycomputed_count,
            "lastcompute_timeout_list": self._timeout_list,
            "lastcompute_errors_list": self._errors_list,
        }
        return ans
