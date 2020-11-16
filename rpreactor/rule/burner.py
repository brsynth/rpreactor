"""
Core code for firing rules
"""

import concurrent.futures
import logging
import sqlite3
import itertools

import pebble
from rdkit import Chem
from rdkit.Chem import AllChem

from rpreactor.chemical.standardizer import Standardizer
from rpreactor.rule.exceptions import ChemConversionError, RuleFireError, RuleConversionError

# TODO follow best practise for logging


def _chunkify(it, size):
    it = iter(it)
    return iter(lambda: tuple(itertools.islice(it, size)), ())


class RuleBurner(object):
    """Apply any number of rules on any number of compounds.

    :param  with_hs:        bool, Enable explicit Hs when sanitizing chemicals
    :param  with_stereo:    bool, Keep stereochemistry (if any) when sanitizing chemicals
    """
    # TODO: rename class to 'RuleManager' or something?

    _SQL_CREATETABLE_MOLECULES = "CREATE TABLE IF NOT EXISTS molecules (id text, rd_mol blob)"
    _SQL_CREATETABLE_RULES = "CREATE TABLE IF NOT EXISTS rules (id text, rd_rule blob)"

    def __init__(self, database=None, with_hs=False, with_stereo=False):
        """Setting up everything needed for behavior decisions and firing rules."""
        # Set-up the database
        if database is None:
            database = ":memory:"
        self._db = sqlite3.connect(database)  # warning: only current thread will be able to use it
        self._db.execute(self._SQL_CREATETABLE_MOLECULES)
        self._db.execute(self._SQL_CREATETABLE_RULES)
        # TODO: add table "results" (watch out for indexes if we search for a previous result at each query!)
        self._db.row_factory = sqlite3.Row
        self._db.commit()
        n_rules = self._db.execute("SELECT count(*) FROM rules").fetchone()['count(*)']
        n_mols = self._db.execute("SELECT count(*) FROM molecules").fetchone()['count(*)']
        logging.info(f"Connected to a database with {n_rules} rules and {n_mols} molecules (at '{database}').")
        # Sanitization
        self._with_hs = with_hs
        self._with_stereo = with_stereo
        if self._with_stereo:
            raise NotImplementedError("Stereo is not implemented at the time being.")

    @staticmethod
    def _task_fire(rd_rule, rd_mol):
        """Apply one reaction rule on one chemical.

        We do not trust the task to terminate by itself, that's why we put it in a separate process that we can kill if
        need be.
        """
        try:
            ans = rd_rule.RunReactants((rd_mol,))
            rdmols, failed = RuleBurner._standardize_results(ans)
            return RuleBurner._handle_results(rdmols)
        except Exception as e:
            raise RuleFireError(e) from e

    @staticmethod
    def _standardize_chemical(rdmol, heavy=False):
        """Simple standardization of RDKit molecule."""
        params = {
            'OP_REMOVE_ISOTOPE': False,
            'OP_NEUTRALISE_CHARGE': False,
            'OP_REMOVE_STEREO': True, #not self._with_stereo, TODO
            'OP_COMMUTE_INCHI': True,
            'OP_KEEP_BIGGEST': False,
            'OP_ADD_HYDROGEN': True, #self._with_hs, TODO
            'OP_KEKULIZE': False,
            'OP_NEUTRALISE_CHARGE_LATE': True
        }
        if heavy:
            params['OP_REMOVE_ISOTOPE'] = True
            params['OP_NEUTRALISE_CHARGE'] = True
            params['OP_KEEP_BIGGEST'] = True
        return Standardizer(sequence_fun='sequence_tunable', params=params).compute(rdmol)

    @staticmethod
    def _standardize_results(tuple_tuple_rdmol):
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
                        list_std.append(RuleBurner._standardize_chemical(rd_frag))
                # Get InChIKey
                for rdmol in list_std:
                    inchikey = Chem.MolToInchiKey(rdmol)
                    if inchikey:
                        list_inchikeys.append(inchikey)
                    else:
                        msg = 'Product conversion to InChIKey raised an empty string'
                        logging.warning(ChemConversionError(msg))
                        raise ChemConversionError(msg)
                # Get unique depiction
                depic = '.'.join(sorted(list_inchikeys))
                # Store only if unique depiction never met
                if depic not in uniq_depics:
                    uniq_depics.add(depic)
                    list_list_std.append(list_std)
            except ChemConversionError as e:
                logging.warning("{}".format(e))
                list_idx_tuple_failed.append(idx_tuple)
                raise e
            except Exception as e:
                logging.warning("Cannot handle a tuple of result, skipped")
                logging.warning("{}".format(e))
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
                logging.warning("{}".format(e))
                raise e
            except Exception as e:
                logging.warning("Cannot handle a tuple of result, skipped")
                logging.warning("{}".format(e))
        return list_list_rdmol, list_list_inchikeys, list_list_inchis, list_list_smiles  # Quick but dirty

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
            rd_mol = RuleBurner._standardize_chemical(rd_mol, heavy=True)
        except Exception as e:
            raise ChemConversionError(e) from e
        return rd_mol

    def _gen_rules(self, ids):
        """Generator of fully initialized/standardized RDKit reaction objects (rules)."""
        for row in self._db.execute("select * from rules;").fetchall():
            yield row['id'], Chem.rdChemReactions.ChemicalReaction(row['rd_rule'])

    def _gen_chemicals(self, ids):
        """Generator of fully initialized/standardized RDKit molecule objects (chemicals)."""
        for row in self._db.execute("select * from molecules;").fetchall():
            yield row['id'], Chem.Mol(row['rd_mol'])

    def _gen_couples(self, rule_list, mol_list):
        """Generator of fully initialized/standardized RDKit couples of rules and molecules."""
        self._gen_rules(rule_list), self._gen_chemicals(mol_list)
        for rid, rd_rule in self._gen_rules(rule_list):
            for cid, rd_mol in self._gen_chemicals(mol_list):
                yield rid, rd_rule, cid, rd_mol

    def _insert_something(self, item_list, id_list, table_name, rdkit_func, chunk_size):
        """Helper function to insert something (metabolite or reaction) into the database."""
        # Consistency
        if id_list and len(item_list) != len(id_list):
            raise ValueError("Provided ID and item lists have different size: aborting.")
        # Standardize items and add them to the database
        n_blocks = (len(item_list)) // chunk_size   # 0-based
        logging.debug(f"Inserting {len(item_list)} RDKit objects in the database as {n_blocks+1} transactions of "
                      f"at most {chunk_size} elements.")
        for chunk_idx, chunk in enumerate(_chunkify(item_list, chunk_size)):
            listof_records = []
            if chunk_idx > 0:
                logging.debug(f"Working on transaction #{chunk_idx+1}...")
            for cindex, item in enumerate(chunk):
                idx = chunk_idx * chunk_size + cindex
                this_id = id_list[idx] if id_list else f"{idx}"
                try:
                    rd_item = rdkit_func(item)
                    listof_records.append((this_id, rd_item.ToBinary()))  # metabolite or reaction
                except ChemConversionError as error:
                    logging.error(f"Something went wrong converting chemical '{this_id}': {error}")
                except RuleConversionError as error:
                    logging.error(f"Something went wrong converting rule '{this_id}': {error}")
            self._db.executemany(f"insert into {table_name} values (?,?)", listof_records)
        self._db.commit()

    def insert_rsmarts(self, rsmarts_list, rid_list=None, chunk_size=10000):
        """Insert reaction rules defined as reaction SMARTS into the database."""
        self._insert_something(rsmarts_list, rid_list, "rules", self._init_rdkit_rule, chunk_size)

    def insert_inchi(self, inchi_list, cid_list=None, chunk_size=10000):
        """Insert molecules defined as InChI into the database."""
        self._insert_something(inchi_list, cid_list, "molecules", self._init_rdkit_mol_from_inchi, chunk_size)

    def dump_to_sql(self, path):
        """Dump the database as a SQL file."""
        with open(path, 'w') as f:
            for line in self._db.iterdump():
                f.write(f'{line}\n')

    def drop_results(self):
        """Drop the results table."""
        raise NotImplementedError  # TODO

    def compute(self, rule_list=None, mol_list=None, commit=False, max_workers=1, timeout=60, chunk_size=1000):
        """Apply all rules on all chemicals and returns a generator over the results."""
        # NB: parallelization will be useless, the time-consuming part is rule and chemical initialization
        with pebble.ProcessPool(max_workers=max_workers) as pool:
            # Prepare chunks of tasks
            # NB: it seems that pool.map does not avoid tasks to hold resources (memory) until they are consumed
            # even if a generator is used as input; so we use pool.schedule and we do our own chunks to avoid saturating
            # the RAM.
            logging.debug(f"Computing tasks in chunks of at most {chunk_size} couples (rule,  molecule) "
                          f"with {max_workers} workers and a per-task timeout of {timeout} seconds.")
            for chunk_idx, chunk in enumerate(_chunkify(self._gen_couples(rule_list, mol_list), chunk_size)):
                logging.debug(f"Working on task chunk #{chunk_idx+1}...")
                # Submit all the tasks for this chunk
                all_running_tasks = []  # list of Future objects
                for rid, rd_rule, cid, rd_mol in chunk:
                    task = (rid, cid, pool.schedule(RuleBurner._task_fire, args=(rd_rule, rd_mol), timeout=timeout))
                    all_running_tasks.append(task)
                # Gather the results
                for i, (rid, cid, future) in enumerate(all_running_tasks):
                    try:
                        rd_mol_list, inchikeys, inchis, smiles = future.result()
                        result = {
                            'rule_id': rid,
                            'substrate_id': cid,
                            'product_list': rd_mol_list,
                            'product_inchikeys': inchikeys,
                            'product_inchis': inchis,
                            'product_smiles': smiles,
                        }
                        if rd_mol_list:  # silently discard tasks without a match
                            yield result
                    except concurrent.futures.TimeoutError:
                        logging.warning(f"Task {rid} on {cid} (#{i}) timed-out.")
                        # task['future'].cancel()  # NB: no need to cancel it, it's already canceled
                    except RuleFireError as error:
                        logging.error(f"Task {rid} on {cid} (#{i}) failed: {error}.")
                    except pebble.ProcessExpired as error:
                        logging.critical(f"Task {rid} on {cid} (#{i}) crashed unexpectedly: {error}.")


def _create_db_from_retrorules_v1_0_5(path_retrosmarts_tsv, path_sqlite):
    raise NotImplementedError  # TODO


def create_db_from_retrorules(path_retrosmarts_tsv, path_sqlite, version="v1.0"):
    """Convert a RetroRules dataset to a rpreactor-ready sqlite3 database.

    For more information on RetroRules, see https://retrorules.org/.
    """
    if version.startswith("v1.0"):
        _create_db_from_retrorules_v1_0_5(path_retrosmarts_tsv, path_sqlite)
