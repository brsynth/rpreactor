"""
Core code for firing rules
"""

import os
import sys
import gzip
import json
import time
import signal
import logging

from rdkit import Chem
from rdkit.Chem import AllChem
import multiprocessing as mp

from rpreactor.chemical.standardizer import Standardizer
from .exceptions import ChemConversionError, RuleMatchError, RuleFireError, RuleConversionError


class RuleBurnerCore(object):
    """Apply one rule on one chemical.""" 

    def __init__(self, rd_rule, rd_mol):
        """Apply one rule on one chemical.
        
        Notice: no standardization is made on inputed chemicals and rules.
        
        :param  rd_rule:    RDKit reaction object, reactio rule to apply
        :param  rd_mol:     RDKit mol object, chemical
        """
        # Internal settings
        USE_CHIRALITY_IN_MATCH = False  # Default value anyway substrucre matching
        # Input
        self._rd_rule = rd_rule
        self._rd_mol = rd_mol
    
    def match(self):
        """Check if left reaction side match the chemical.
        
        returns:    bool, True if there is a match, else False
        """
        try:
            for reactant in self._rd_rule.GetReactants():
                if self._rd_mol.HasSubstructMatch(reactant, ):
                    return True
            return False
        except Exception as e:
            raise RuleMatchError(e) from e
        
    def fire(self):
        """Fire the rule on the chemical.
        
        returns:    tuple of tuple, list of results for each possible application.
        """
        try:
            return self._rd_rule.RunReactants((self._rd_mol,))
        except Exception as e:
            raise RuleFireError(e) from e


class RuleBurner(object):
    """Apply any number of rules on any number of compounds."""

    def __init__(
            self, rsmarts_list, inchi_list, rid_list=None,  cid_list=None,
            match_timeout=1, fire_timeout=1, ofile=None, compress=False,
            with_hs=False, with_stereo=False):
        """Setting up everything needed for behavor decisions and firing rules.

        :param  rsmarts_list:   list of reaction rule SMARTS
        :param  inchi_list:     list of inchis
        :param  rid_list:       list of reaction rule IDs
        :param  cid_list:       list of chemical IDs
        :param  match_timeout:  int, timeout execution for compound pre-matching
        :param  fire_timeout:   int, timeout execution for rule firing
        :param  ofile:          str, Output file to store results
        :param  compress:       bool, enable gzip compression on output
        :param  with_hs:        bool, Enable explicit Hs when sanitizing chemicals
        :param  with_stereo:    bool, Keep stereochemistry (if any) when sanitizing chemicals
        """

        # Internal settings
        self._INDENT_JSON = True
        self._TRY_MATCH = False

        # Input
        self._rsmarts_list = rsmarts_list
        self._rid_list = rid_list
        self._inchi_list = inchi_list
        self._cid_list = cid_list

        # Settings
        self._try_match = self._TRY_MATCH  # TODO: add option for that
        self._match_timeout = match_timeout
        self._fire_timeout = fire_timeout

        # Sanitization
        self._with_hs = with_hs
        self._with_stereo = with_stereo
        if self._with_stereo:
            raise NotImplementedError("Stereo is not implemented at the time being.")

        # Check for consistency between depictions and IDs
        if self._rid_list and len(self._rsmarts_list) != len(self._rid_list):
            logging.warning("ID and depiction rule lists have different size, compound IDs will be ignored")

        if self._cid_list and len(self._inchi_list) != len(self._cid_list):
            logging.warning("ID and depiction compounds lists have different size, rule IDs will be ignored")

        # Output
        self._json = list()
        self._compress = compress
        if not ofile:
            self._ofile = None
        else:
            pdir = os.path.abspath(os.path.dirname(ofile))
            os.makedirs(pdir, exist_ok=True)
            self._ofile = os.path.abspath(ofile)

        # A place to swim
        self._pool = None

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

    def _run_with_timeout(self, worker, kwargs, timeout=5):
        """Generic wrapper making use of multiprocessing to garantee effective timeout.

        :param  worker:     function to be called
        :param  kwargs:     dictionnary of args to be passed to the called function
        :param  timeout:    int, timeout
        :returns            depends of the worker function
        """
        if self._pool is None:
            self._pool = mp.Pool(processes=1)
        try:
            start_time = time.time()
            res = self._pool.apply_async(worker, kwds={'kwargs': kwargs})
            ans = res.get(timeout=timeout)
            end_time = time.time()
            exec_time = round(end_time - start_time, 4)
        except mp.TimeoutError as e:
            _kill(self._pool)
            self._pool = mp.Pool(processes=1)
            raise e
        except Exception as e:
            raise e
        return ans, exec_time

    def _jsonify(self, rid=None, cid=None,
                 has_match=None, match_timed_out=None,
                 match_exec_time=None, match_error=None,
                 fire_timed_out=None,
                 fire_exec_time=None, fire_error=None,
                 inchikeys_list=None, inchis_list=None, smiles_list=None):
        """Return the results as a JSON string.

        :param      rid:                str, reaction rule ID
        :param      cid:                str, substrate ID
        :param      has_match:          bool or None, True if there is match
        :param      match_timed_out:    bool, True if timeout reached
        :param      match_exec_time:    int, execution time for matching
        :param      match_error:        str, error message if any, else None
        :param      fire_timed_out:     bool, True if timeout reached
        :param      fire_exec_time:     bool, execution time for firing
        :param      fire_error:         str, error message if any, else None
        :param      inchikeys_list:     list of list, Inchikeys of products
        :param      inchis_list:        list of list, Inchis of products
        :param      smiles_list:        list of list, SMILES of products
        :returns    json_string:        JSON string
        """
        # General info
        data = {
            'rule_id': rid,
            'substrate_id': cid,
        }
        # Match info
        if self._try_match:
            data['match'] = has_match
            data['match_timed_out'] = match_timed_out
            data['match_exec_time'] = match_exec_time
            if match_error is not None:
                data['match_error'] = match_error
        # Fire info
        data['fire_timed_out'] = fire_timed_out
        data['fire_exec_time'] = fire_exec_time
        if fire_error is not None:
            data['fire_error'] = fire_error
        if (inchikeys_list is not None) and (len(inchikeys_list) > 0):
            data['product_inchikeys'] = inchikeys_list
            data['product_inchis'] = inchis_list
            data['product_smiles'] = smiles_list

        return json.dumps(obj=data, indent=self._INDENT_JSON)

    def write_json(self):
        """Write the JSON string."""
        # Handling file handler
        if self._ofile:
            if self._compress:
                ofh = gzip.open(self._ofile, 'wb', compresslevel=9)
            else:
                ofh = open(self._ofile, 'w')
        else:
            ofh = sys.stdout
        # Big string
        content = '[\n' + ','.join(self._json) + '\n]' + '\n'
        # Handling compression
        if self._ofile and self._compress:
            ofh.write(content.encode())
        else:
            ofh.write(content)
        ofh.close()

    def compute(self):
        """Rules under fire."""
        # TODO: simplify
        for rindex, rsmarts in enumerate(self._rsmarts_list):
            # Extract corresponding reaction rule ID if any
            if self._rid_list:
                rid = self._rid_list[rindex]
            else:
                rid = None
            # Get RDKit reaction object
            try:
                rd_rule = AllChem.ReactionFromSmarts(rsmarts)
                rd_rule.Initialize()
            except Exception as e:
                raise RuleConversionError(e) from e

            for cindex, inchi in enumerate(self._inchi_list):
                # Extract corresponding substrate ID if any
                if self._cid_list:
                    cid = self._cid_list[cindex]
                else:
                    cid = None
                # Get standardized RDKit mol
                try:
                    rd_mol = Chem.MolFromInchi(inchi, sanitize=False)  # Important: Sanitize = False
                    rd_mol = self._standardize_chemical(rd_mol, heavy=True)
                except Exception as e:
                    logging.warning(e)
                    raise ChemConversionError(e) from e
                # General args to used for both matching and firing
                kwargs = {
                        'rd_rule': rd_rule,
                        'rd_mol': rd_mol
                        }
                # Matching
                has_match = None
                match_exec_time = None
                match_timed_out = None
                match_error = None
                if self._try_match:
                    try:
                        has_match, match_exec_time = self._run_with_timeout(
                                worker=_worker_match, kwargs=kwargs,
                                timeout=self._match_timeout
                                )
                        match_timed_out = False
                    except mp.TimeoutError as e:
                        match_timed_out = True
                        match_error = str(e)
                    except Exception as e:
                        match_timed_out = False
                        match_error = str(e)
                # Firing
                fire_exec_time = None
                fire_timed_out = None
                fire_error = None
                inchikeys = None
                inchis = None
                smiles = None
                try:
                    ans, fire_exec_time = self._run_with_timeout(
                            worker=_worker_fire, kwargs=kwargs,
                            timeout=self._fire_timeout
                            )
                    rdmols, failed = self._standardize_results(ans)
                    inchikeys, inchis, smiles = self._handle_results(rdmols)
                    fire_timed_out = False
                except ChemConversionError as e:
                    fire_timed_out = False
                    fire_error = str(e)
                    logging.warning(e)
                except mp.TimeoutError as e:
                    fire_timed_out = True
                    fire_error = str(e)
                    logging.error('TIMEOUT: cid={}, rid={}'.format(cid, rid))
                    logging.error('TIMEOUT: original error={}'.format(e))
                except Exception as e:
                    fire_timed_out = False
                    fire_error = str(e)
                    logging.warning(e)
                # JSONify and store
                json_str = self._jsonify(
                        rid=rid,
                        cid=cid,
                        has_match=has_match,
                        match_timed_out=match_timed_out,
                        match_exec_time=match_exec_time,
                        match_error=match_error,
                        fire_timed_out=fire_timed_out,
                        fire_exec_time=fire_exec_time,
                        fire_error=fire_error,
                        inchikeys_list=inchikeys,
                        inchis_list=inchis,
                        smiles_list=smiles
                        )
                self._json.append(json_str)
        # Clean the pool
        if self._pool:
            self._pool.terminate()
            self._pool.join()


def _worker_match(kwargs):
    """Check if a chemical can be fired by a rule according to left side."""
    w = RuleBurnerCore(**kwargs)
    return w.match()


def _worker_fire(kwargs):
    """Apply a reaction a rule on a chemical."""
    r = RuleBurnerCore(**kwargs)
    return r.fire()


def _kill(pool):
    """Send SIGTERMs to kill all processes belonging to pool.

    Will not work on Windows OS.
    """
    # stop repopulating new child
    pool._state = mp.pool.TERMINATE
    pool._worker_handler._state = mp.pool.TERMINATE
    for p in pool._pool:
        os.kill(p.pid, signal.SIGKILL)
    # .is_alive() will reap dead process
    while any(p.is_alive() for p in pool._pool):
        pass
    # Get-lucky workaround: force releasing lock
    try:
        pool._inqueue._rlock.release()
    except ValueError as e:
        logging.error(e)
    pool.terminate()
