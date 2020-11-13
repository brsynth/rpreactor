"""
Core code for firing rules
"""

import concurrent.futures
import logging
import sqlite3
import itertools as it

import pebble
from rdkit import Chem
from rdkit.Chem import AllChem

from rpreactor.chemical.standardizer import Standardizer
from rpreactor.rule.exceptions import ChemConversionError, RuleFireError, RuleConversionError


def _task_fire(rd_rule, rd_mol):
    """Apply a reaction a rule on a chemical."""
    try:
        return rd_rule.RunReactants((rd_mol,))
    except Exception as e:
        raise RuleFireError(e) from e


class RuleBurner(object):
    """Apply any number of rules on any number of compounds.

    :param  with_hs:        bool, Enable explicit Hs when sanitizing chemicals
    :param  with_stereo:    bool, Keep stereochemistry (if any) when sanitizing chemicals
    """

    def __init__( self, with_hs=False, with_stereo=False):
        """Setting up everything needed for behavior decisions and firing rules."""
        # Set-up the database
        # TODO select from file and fallback to :memory: if non selected
        self._db = sqlite3.connect(":memory:")
        c = self._db.cursor()
        c.execute("CREATE TABLE molecules (id text, mol text, rd_mol blob)")
        c.execute("CREATE TABLE rules (id text, rule text, rd_rule blob)")
        self._db.row_factory = sqlite3.Row
        self._db.commit()

        # Sanitization
        self._with_hs = with_hs
        self._with_stereo = with_stereo
        if self._with_stereo:
            raise NotImplementedError("Stereo is not implemented at the time being.")

    def _standardize_chemical(self, rdmol, heavy=False):
        """Simple standardization of RDKit molecule."""
        params = {
            'OP_REMOVE_ISOTOPE': False,
            'OP_NEUTRALISE_CHARGE': False,
            'OP_REMOVE_STEREO': not self._with_stereo,
            'OP_COMMUTE_INCHI': True,
            'OP_KEEP_BIGGEST': False,
            'OP_ADD_HYDROGEN': self._with_hs,
            'OP_KEKULIZE': False,
            'OP_NEUTRALISE_CHARGE_LATE': True
        }
        if heavy:
            params['OP_REMOVE_ISOTOPE'] = True
            params['OP_NEUTRALISE_CHARGE'] = True
            params['OP_KEEP_BIGGEST'] = True
        return Standardizer(sequence_fun='sequence_tunable', params=params).compute(rdmol)

    def _standardize_results(self, tuple_tuple_rdmol):
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
                        list_std.append(self._standardize_chemical(rd_frag))
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
        return list_list_inchikeys, list_list_inchis, list_list_smiles  # Quick but dirty

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
            rd_mol = self._standardize_chemical(rd_mol, heavy=True)
        except Exception as e:
            raise ChemConversionError(e) from e
        return rd_mol

    def _gen_rules(self):
        """Generator of fully initialized/standardized RDKit reaction objects (rules)."""
        for row in self._db.execute("select * from rules;").fetchall():
            yield row['id'], Chem.rdChemReactions.ChemicalReaction(row['rd_rule'])

    def _gen_chemicals(self):
        """Generator of fully initialized/standardized RDKit molecule objects (chemicals)."""
        for row in self._db.execute("select * from molecules;").fetchall():
            yield row['id'], Chem.Mol(row['rd_mol'])

    def _insert_something(self, item_list, id_list, where, rdkit_func, chunk_size):
        """Helper function to insert something (metabolite or reaction) into the database."""
        # Consistency
        if id_list and len(item_list) != len(id_list):
            logging.warning("provided ID and item lists have different size: IDs will be ignored")
            id_list = None
        # Standardize molecules and add them to the database
        logging.debug(f"Inserting {len(item_list)} RDKit objects in the database. Will display progress every 1000 conversions.")
        n_blocks = len(item_list) // chunk_size + 1
        for bindex, item_block in enumerate(it.tee(item_list, n_blocks)):
            listof_records = []
            for cindex, item in enumerate(item_block):
                cindex += bindex
                this_id = id_list[cindex] if id_list else f"{cindex}"
                if cindex % 1000 == 0:
                    logging.debug(f"Working on #{cindex} ({this_id})...")
                try:
                    rd_item = rdkit_func(item)
                    listof_records.append((this_id, item, rd_item.ToBinary()))
                except ChemConversionError as error:
                    logging.error(f"Something went wrong converting chemical '{this_id}': {error}")
                except RuleConversionError as error:
                    logging.error(f"Something went wrong converting rule '{this_id}': {error}")
            self._db.executemany(f"insert into {where} values (?,?,?)", listof_records)
        self._db.commit()

    def insert_inchi(self, inchi_list, cid_list=None, chunk_size=10000):
        """Register molecules defined as InChI."""
        self._insert_something(inchi_list, cid_list, "molecules", self._init_rdkit_mol_from_inchi, chunk_size)

    def insert_rsmarts(self, rsmarts_list, rid_list=None, chunk_size=10000):
        """Register reactions defined as resction SMARTS."""
        self._insert_something(rsmarts_list, rid_list, "rules", self._init_rdkit_rule, chunk_size)

    def compute(self, max_workers=1, timeout=60):
        """Apply all rules on all chemicals and returns a generator over the results."""
        # NB: parallelization will be useless, the time-consuming part is rule and chemical initialization
        with pebble.ProcessPool(max_workers=max_workers) as pool:
            # Submit all the tasks
            # NB: it seems that pool.map does not avoid tasks to hold resources (memory) until they are consumed
            # when a generator is used as input; so we use pool.schedule and we explicitly store parameters
            # (i.e. RDKit objects) for readability.
            all_running_tasks = []  # list of Future objects
            for rid, rd_rule in self._gen_rules():
                for cid, rd_mol in self._gen_chemicals():
                    task = (rid, cid, pool.schedule(_task_fire, args=(rd_rule, rd_mol), timeout=timeout))
                    all_running_tasks.append(task)
            # Gather the results
            logging.debug(f"Found {len(all_running_tasks)} tasks. Will display progress every 1000 tasks.")
            for i, (rid, cid, future) in enumerate(all_running_tasks):
                if i % 1000 == 0:
                    logging.debug(f"Working on task #{i} ({rid} on {cid})...")
                result = {
                    'rule_id': rid,
                    'substrate_id': cid,
                    'product_inchikeys': None,
                    'product_inchis': None,
                    'product_smiles': None,
                }
                try:
                    ans = future.result()
                    rdmols, failed = self._standardize_results(ans)
                    result['product_inchikeys'], result['product_inchis'], result['product_smiles'] = self._handle_results(rdmols)
                    yield result
                except concurrent.futures.TimeoutError:
                    logging.warning(f"Task {rid} on {cid} (#{i}) timed-out.")
                    # task['future'].cancel()  # NB: no need to cancel it, it's already canceled
                except RuleFireError as error:
                    logging.error(f"Task {rid} on {cid} (#{i}) failed: {error}.")
                except pebble.ProcessExpired as error:
                    logging.critical(f"Task {rid} on {cid} (#{i}) crashed unexpectedly: {error}.")
