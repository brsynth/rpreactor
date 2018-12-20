r"""
Apply one rule on one substrate.
"""

import os
import sys
import json
import time
import rdkit
import signal
import logging
import argparse

from rdkit import Chem
from rdkit.Chem import AllChem

import multiprocessing as mp


class RuleBurnerError(Exception):
    """Home made exception."""
    pass
    
class RuleConversionError(RuleBurnerError):
    """Raised when something went wrong during SMARTS conversion to RDKit rxn object."""
    
    def __init__(self, msg):
        self._msg = msg

    def __str__(self):
        return "RULE-CONVERSION-ERROR: {}.".format(self._msg)

class RuleFireError(RuleBurnerError):
    """Raised when something went wrong when firing a rule."""

    def __init__(self, msg):
        self._msg = msg

    def __str__(self):
        return "RULE-FIRE-ERROR: {}.".format(self._msg)

class RuleMatchError(RuleBurnerError):
    """Raised when something went wrong when matching a rule."""

    def __init__(self, msg):
        self._msg = msg

    def __str__(self):
        return "RULE-MATCH-ERROR: {}.".format(self._msg)

class ChemConversionError(RuleBurnerError):
    """Raised when something went wrong during chemical conversion to RDKit mol object and sanitization."""
    
    def __init__(self, msg):
        self._msg = msg

    def __str__(self):
        return "CHEM-CONVERSION-ERROR: {}.".format(self._msg)

def worker_match(kwargs):
    """Check if a chemical can be fired by a rule according to left side."""
    w = RuleBurnerCore(**kwargs)
    return w.match()


def worker_fire(kwargs):
    """Apply a reaction a rule on a chemical."""
    r = RuleBurnerCore(**kwargs)
    return r.fire()


def kill(pool):
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
    print(pool)
    pool.terminate()


def worker_test(maxtime=5):
    time.sleep(maxtime)


class RuleBurnerCore(object):
    """Apply one rule on one chemical.""" 

    def __init__(self, rd_rule, rd_mol):
        """Apply one rule on one chemical.
        
        Notice: no standardization is made on inputed chemicals and rules.
        
        :param  rd_rule:    RDKit reaction object, reactio rule to apply
        :param  rd_mol:     RDKit mol object, chemical
        :param  timeout:    str, Reaction rule SMARTS
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
            self, rsmarts_list, csmiles_list, rid_list=None,  cid_list=None,
            match_timeout=1, fire_timeout=1, ofile=None):
        """Setting up everything needed for behavor decisions and firing rules.
        
        :param  rsmarts_list:   list of reaction rule SMARTS
        :param  csmiles_list:   list of compound SMILES
        :param  rid_list:       list of reaction rule IDs
        :param  cid_list:       list of compound IDs
        :param  match_timeout:  int, timeout execution for compound pre-matching
        :param  fire_timeout:   int, timeout execution for rule firing 
        :param  ofile:          str, Output file to store results
        """
        
        # Internal settigns
        self._STRIP_PARENTHESIS = True  # Surrounding parenthesis should be removed for left side matching
        self._INDENT_JSON = True
        
        # Input
        self._rsmarts_list = rsmarts_list
        self._rid_list = rid_list
        self._csmiles_list = csmiles_list
        self._cid_list = cid_list
    
        # Settings
        self._match_timeout = match_timeout
        self._fire_timeout = fire_timeout
        
        # Check for consistency between depictions and IDs
        try:
            if self._rid_list:
                assert len(self._rsmarts_list) == len(self._rid_list)
        except AssertionError as e:
            logging.warning("ID and depiction rule lists have different size, compound IDs will be ignored")

        try:
            if self._cid_list:
                assert len(self._csmiles_list) == len(self._cid_list)
        except AssertionError as e:
            logging.warning("ID and depiction compounds lists have different size, rule IDs will be ignored")
        
        # Output
        self._json = list()
        if not ofile:
            self._ofile = None
        else:
            pdir = os.path.abspath(os.path.dirname(ofile))
            os.makedirs(pdir, exist_ok=True)
            self._ofile = os.path.abspath(ofile)
        
        # A place to swim
        self._pool = None
                
    def _standardize_chemical(self, rd_mol):
        """Standardize a chemical using RDKit sanitize method.

        :param      rd_mol:  RDKit mol object
        :returns    nothing, in place operation
        """
        try:
            Chem.SanitizeMol(rd_mol)
        except ValueError as e:  # In case kekulize lead to some issues
            saniFlag = Chem.SanitizeFlags.SANITIZE_ALL
            saniFlag ^= Chem.SanitizeFlags.SANITIZE_KEKULIZE
            Chem.SanitizeMol(rd_mol, sanitizeOps=saniFlag)
            logging.warning('Partial sanization only')
        except Exception as e:
            raise e

    def _run_with_timeout(self, worker, kwargs, timeout=5):
        """Generic wrapper for using multiprocessing to garantee effective timeout.
        
        :param  worker:     function to be called
        :param  kwargs:     dictionnary of args to be passed to the called function
        :param  timeout:    int, timeout
        :param  pool:       multiprocessing Pool object, for optimization purpose
        :returns            depends of the worker function 
        """
        if self._pool is None:
            self._pool = mp.Pool(processes=1)
        try:
            start_time = time.time()
            res = self._pool.apply_async(worker, kwds={'kwargs': kwargs})
            ans = res.get(timeout=timeout)
            end_time = time.time()
            exec_time = round(end_time - start_time, 3)
        except mp.TimeoutError as e:
            kill(pool)
            self._pool = Pool(processes=1)
            raise e
        except Exception as e:
            raise e
        return ans, exec_time

    def _handle_results(self, tuple_tuple_raw):
        """Perform sanitization, InChI and SMILES generation and remove duplicates.
        
        :param      tuple_tuple_raw:        tuple of tuple of RDKit Mol
        :returns    tuple_tuple_inchis:     tuple of tuple of InChIs
        :returns    tuple_tuple_smiles:     tuple of tuple of SMILES
        """
        # Collect and standardize
        list_list_inchis = list()
        list_list_smiles = list()
        uniq_depics = set()
        for tuple_raw in tuple_tuple_raw:
            try: 
                list_inchis = list()
                list_smiles = list()
                list_std = list()
                # Standadize
                for rd_mol in tuple_raw:
                    for rd_frag in Chem.GetMolFrags(rd_mol, asMols=True, sanitizeFrags=False):
                        self._standardize_chemical(rd_frag)
                        list_std.append(rd_frag)
                # Get InChIs
                for rd_mol in list_std:
                    list_inchis.append(Chem.MolToInchi(rd_mol))
                # Get unique depiction
                depic = '.'.join(sorted(list_inchis))
                # Continue only if depiction never met
                if depic in uniq_depics:
                    continue
                # Get SMILES
                for rd_mol in list_std:
                    list_smiles.append(Chem.MolToSmiles(rd_mol))
                # Finally store those that reach the end
                list_list_inchis.append(list_inchis)
                list_list_smiles.append(list_smiles)
            except Exception as e:
                logging.warning("Cannot handle a tuple of result, skipped")
                logging.warning("{}".format(e))
        return list_list_smiles, list_list_inchis  # Quick and dirty


    def _jsonify(self, rsmarts, csmiles, rid=None, cid=None,
                 has_match=None, match_timed_out=None,
                 match_exec_time=None, match_error=None,
                 fire_timed_out=None, fire_exec_time=None, fire_error=None,
                 smiles_list=None, inchis_list=None):
        """Return the results as a json string."""
        
        data = {
            # General info
            'rule_id': rid,
            'rule_smarts': rsmarts,
            'substrate_id': cid,
            'substrate_smiles': csmiles,
            # Match info
            'match': has_match,
            'match_exec_time': match_exec_time,
            'match_timed_out': match_timed_out,
            'match_error': match_error,
            # Firing info
            'fire_timed_out': fire_timed_out,
            'fire_exec_time': fire_exec_time,
            'fire_error': fire_error,
            # Add product InChIs, SMILESs, ..
            'product_smiles': smiles_list,
            'product_inchis': inchis_list
        }
        
        # JSONify
        return json.dumps(obj=[data], indent=self._INDENT_JSON)

    def write_json(self):
        """Write the JSON string."""
        ofh = open(self._ofile, 'w') if self._ofile else sys.stdout
        ofh.write(','.join(self._json))
        ofh.write('\n')
        ofh.close()

    def compute(self):
        """Rules under fire."""
        
        for rindex, rsmarts in enumerate(self._rsmarts_list):
            if self._rid_list:
                rid = self._rid_list[rindex]
            else:
                rid = None
            # Get initialized RDKit rxn
            try:
                rd_rule = AllChem.ReactionFromSmarts(rsmarts)
                rd_rule.Initialize()
            except Exception as e:
                raise RuleConversionError(e) from e

            for cindex, csmiles in enumerate(self._csmiles_list):
                if self._cid_list:
                    cid = self._cid_list[cindex]
                else:
                    cid = None
                # Get standardized RDKit mol
                try:
                    rd_mol = Chem.MolFromSmiles(csmiles)
                    self._standardize_chemical(rd_mol)
                except Exception as e:
                    raise ChemConversionError(e) from e
                # General args
                kwargs = {
                        'rd_rule': rd_rule,
                        'rd_mol': rd_mol
                        }
                # Matching
                try:
                    has_match, match_exec_time = self._run_with_timeout(
                            worker=worker_match, kwargs=kwargs,
                            timeout=self._match_timeout
                            )
                    match_timed_out = False
                    match_error = None
                except mp.TimeoutError as e:
                    match_timed_out = True
                    match_error = str(e)
                except Exception as e:
                    match_timed_out = False
                    match_error = str(e)
                # Firing
                try:
                    ans, fire_exec_time = self._run_with_timeout(
                            worker=worker_fire, kwargs=kwargs,
                            timeout=self._fire_timeout
                            )
                    smiles, inchis = self._handle_results(ans)
                    fire_timed_out = False
                    fire_error = None
                except mp.TimeoutError as e:
                    fire_timed_out = True
                    fire_error = str(e)
                except Exception as e:
                    fire_timed_out = False
                    fire_error = str(e)
                # JSONify and store
                json_str = self._jsonify(
                        rsmarts=rsmarts, csmiles=csmiles, rid=rid, cid=cid,
                        has_match=has_match, match_timed_out=match_timed_out,
                        match_exec_time=match_exec_time, match_error=match_error,
                        fire_timed_out=fire_timed_out, fire_exec_time=fire_exec_time,
                        fire_error=fire_error,
                        smiles_list=smiles, inchis_list=inchis
                        )
                self._json.append(json_str)
                
        
        


def __cli():
    """Command line interface."""

    help = "Apply one rule on one substrate."

    parser = argparse.ArgumentParser(description=help)
    parser.add_argument('--rsmarts', help='Reaction rule SMARTS', required=True)
    parser.add_argument('--csmiles', help='Chemical SMILES depiction', required=True)
    parser.add_argument('--rid', help='Reaction rule ID, optional')
    parser.add_argument('--cid', help='Chemical ID, optional')
    parser.add_argument('--match_timeout', help='Rule matching timeout', default=5, type=int)
    parser.add_argument('--fire_timeout', help='Rule furing timeout', default=5, type=int)
    parser.add_argument('--ofile', help='Output file to store results. Default to STDOUT if none provided')
    args = parser.parse_args()

    # Logging
    logging.basicConfig(
            stream=sys.stderr, level=logging.INFO,
            datefmt='%d/%m/%Y %H:%M:%S',
            format='%(asctime)s -- %(levelname)s -- %(message)s'
            )

    r = RuleBurner(
            rsmarts_list=[args.rsmarts], csmiles_list=[args.csmiles],
            rid_list=[args.rid], cid_list=[args.cid],
            match_timeout=args.match_timeout, fire_timeout=args.fire_timeout,
            ofile=args.ofile
            )
    
    r.compute()
    r.write_json()

if __name__ == "__main__":
    __cli()
