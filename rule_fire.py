r"""
Apply one/many rules on one/many substrates.

Thomas Duigou, INRA, 2018
"""

import os
import sys
import gzip
import json
import time
import rdkit
import signal
import logging
import argparse
import multiprocessing as mp

from rdkit import Chem
from rdkit.Chem import AllChem



class RuleBurnerError(Exception):
    """Home made exception."""
    pass
    
class RuleConversionError(RuleBurnerError):
    """Raised when something went wrong during SMARTS conversion to RDKit rxn object."""
    
    def __init__(self, msg):
        self._msg = msg

    def __str__(self):
        return "RULE-CONVERSION-ERROR: {}".format(self._msg)

class RuleFireError(RuleBurnerError):
    """Raised when something went wrong when firing a rule."""

    def __init__(self, msg):
        self._msg = msg

    def __str__(self):
        return "RULE-FIRE-ERROR: {}".format(self._msg)

class RuleMatchError(RuleBurnerError):
    """Raised when something went wrong when matching a rule."""

    def __init__(self, msg):
        self._msg = msg

    def __str__(self):
        return "RULE-MATCH-ERROR: {}".format(self._msg)

class ChemConversionError(RuleBurnerError):
    """Raised when something went wrong during chemical conversion to RDKit mol object and sanitization."""
    
    def __init__(self, msg):
        self._msg = msg

    def __str__(self):
        return "CHEM-CONVERSION-ERROR: {}".format(self._msg)

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
    pool.terminate()


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
            match_timeout=1, fire_timeout=1, ofile=None, compress=False):
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
        self._compress = compress
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
            kill(self._pool)
            self._pool = mp.Pool(processes=1)
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
        list_list_inchis = list()
        list_list_smiles = list()
        uniq_depics = set()
        
        for tuple_raw in tuple_tuple_raw:
            try: 
                list_inchis = list()
                list_smiles = list()
                list_std = list()
                # Standardize
                for rd_mol in tuple_raw:
                    for rd_frag in Chem.GetMolFrags(rd_mol, asMols=True, sanitizeFrags=False):
                        self._standardize_chemical(rd_frag)
                        list_std.append(rd_frag)
                # Get InChIs
                for rd_mol in list_std:
                    inchi = Chem.MolToInchi(rd_mol)
                    if inchi: 
                        list_inchis.append(inchi)
                    else:
                        raise ChemConversionError("Product conversion to InChI raised an empty string")
                # Get unique depiction
                depic = '.'.join(sorted(list_inchis))
                # Continue only if depiction never met
                if depic in uniq_depics:
                    continue
                uniq_depics.add(depic)
                # Get SMILES
                for rd_mol in list_std:
                    list_smiles.append(Chem.MolToSmiles(rd_mol))
                # Finally store those that reach the end
                list_list_inchis.append(list_inchis)
                list_list_smiles.append(list_smiles)
            except ChemConversionError as e:
                logging.warning("{}".format(e))
                raise e
            except Exception as e:
                logging.warning("Cannot handle a tuple of result, skipped")
                logging.warning("{}".format(e))
        return list_list_smiles, list_list_inchis  # Quick but dirty


    def _jsonify(self, rsmarts, csmiles, rid=None, cid=None,
                 has_match=None, match_timed_out=None,
                 match_exec_time=None, match_error=None,
                 fire_timed_out=None, fire_exec_time=None, fire_error=None,
                 smiles_list=None, inchis_list=None):
        """Return the results as a JSON string.
        
        :param      rsmarts:            str, reaction rule string depiction
        :param      csmiles:            str, substrate string depiction
        :param      rid:                str, reaction rule ID
        :param      cid:                str, substrate ID
        :param      has_match:          bool or None, True if there is match
        :param      match_timed_out:    bool, True if timeout reached
        :param      match_exec_time:    int, execution time for matching
        :param      match_error:        str, error message if any, else None
        :param      fire_timed_out:     bool, True if timeout reached
        :param      fire_exec_time:     bool, execution time for firing
        :param      fire_error:         str, error message if any, else None
        :param      smiles_list:        list of list, SMILES of products
        :param      inchis_list:        list of list, SMILES of InChIs
        :returns    json_string:        JSON string
        """
        # General info
        data = {
            'rule_id': rid,
            # 'rule_smarts': rsmarts,
            'substrate_id': cid,
            # 'substrate_smiles': csmiles,
        }
        # Match info
        data['match'] = has_match
        data['match_timed_out'] = match_timed_out
        data['match_exec_time'] = match_exec_time
        if match_error is not None:
            data['match_error'] = match_error
        # Fire info
        if (smiles_list is None) or (len(smiles_list) > 0):
            data['fire'] = True
        else:
            data['fire'] = False
        data['fire_timed_out'] = fire_timed_out
        data['fire_exec_time'] = fire_exec_time
        if fire_error is not None:
            data['fire_error'] = fire_error
        if (smiles_list is None) or (len(smiles_list) > 0):
            data['product_smiles'] = smiles_list
            data['product_inchis'] = inchis_list

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

            for cindex, csmiles in enumerate(self._csmiles_list):
                # Extract corresponding substrate ID if any
                if self._cid_list:
                    cid = self._cid_list[cindex]
                else:
                    cid = None
                # Get standardized RDKit mol
                try:
                    rd_mol = Chem.MolFromSmiles(csmiles, sanitize=False)  # Important: Sanitize = False
                    self._standardize_chemical(rd_mol)
                except Exception as e:
                    raise ChemConversionError(e) from e
                # General args to used for both matching and firing
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
                    has_match = None
                    match_exec_time = None
                    match_timed_out = True
                    match_error = str(e)
                except Exception as e:
                    has_match = None
                    match_exec_time = None
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
                except ChemConversionError as e:
                    smiles = None
                    inchis = None
                    fire_timed_out = False
                    fire_error = str(e)
                except mp.TimeoutError as e:
                    fire_exec_time = None
                    smiles = None
                    inchis = None
                    fire_timed_out = True
                    fire_error = str(e)
                except Exception as e:
                    fire_exec_time = None
                    smiles = None
                    inchis = None
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

    help = "Apply rules on chemicals."

    def inline_mode(args):
        """Execution mode to be used when a single rule and a single chemical
        are porvided through CLI.
        """
        r = RuleBurner(
                rsmarts_list=[args.rsmarts], csmiles_list=[args.csmiles],
                rid_list=[args.rid], cid_list=[args.cid],
                match_timeout=args.match_timeout, fire_timeout=args.fire_timeout,
                ofile=args.ofile, compress=args.compress
                )
        r.compute()
        r.write_json()

    def infile_mode(args):
        """Execution mode to be used when rules and chemicals are provided 
        in CSV files.
        """
        
        rsmarts_list = list()
        rids_list = list()

        rsmiles_list = list()
        cids_list = list()
        
        import csv
        
        with open(args.rfile, 'r') as ifh:
            reader = csv.DictReader(ifh, delimiter='\t')
            for row in reader:
                rsmarts_list.append(row['Rule_SMARTS'].strip())
                rids_list.append(row['Rule_ID'].strip())

        with open(args.cfile, 'r') as ifh:
            reader = csv.DictReader(ifh, delimiter='\t')
            for row in reader:
                rsmiles_list.append(row['Substrate_SMILES'].strip())
                cids_list.append(row['Substrate_ID'].strip())
                
        r = RuleBurner(
                rsmarts_list=rsmarts_list, csmiles_list=rsmiles_list,
                rid_list=rids_list, cid_list=cids_list,
                match_timeout=args.match_timeout, fire_timeout=args.fire_timeout,
                ofile=args.ofile, compress=args.compress
                )
        r.compute()
        r.write_json()

    parser = argparse.ArgumentParser(description=help)
    parser.add_argument('--match_timeout', help='Rule matching timeout. Default: 1.', default=1, type=int)
    parser.add_argument('--fire_timeout', help='Rule furing timeout. Default: 1.', default=1, type=int)
    parser.add_argument('--ofile', help='Output file to store results. Default to STDOUT if none provided')
    parser.add_argument('--compress', action='store_true', help='Enable gzip compression (only when output to file).')
    subparsers = parser.add_subparsers(help='Input mode')

    parser_inline = subparsers.add_parser('inline', help='Get inputs from command line')
    parser_inline.set_defaults(func=inline_mode)
    parser_inline.add_argument('--rsmarts', help='Reaction rule SMARTS', required=True)
    parser_inline.add_argument('--csmiles', help='Chemical SMILES depiction', required=True)
    parser_inline.add_argument('--rid', help='Reaction rule ID, optional')
    parser_inline.add_argument('--cid', help='Chemical ID, optional')
    
    parser_file = subparsers.add_parser('infile', help='Get inputs from files')
    parser_file.set_defaults(func=infile_mode)
    parser_file.add_argument(
            '--rfile', required=True, help=' '.join([
                    'Reaction rule file.',
                    'Tab separated columns.',
                    'One reaction rule per line.',
                    'Mandatory column: Rule_SMARTS.',
                    'Optional column: Rule_ID.',
                    'Other columns will be ignored.'
                    ])
            )
    parser_file.add_argument(
            '--cfile', required=True, help=' '.join([
                    'Chemical file.',
                    'Tab separated columns.',
                    'One chemical per line.',
                    'Mandatory column: Substrate_SMILES.',
                    'Optional column: Substrate_ID.',
                    'Other columns will be ignored.'
                    ])
            )

    # Logging
    logging.basicConfig(
            stream=sys.stderr, level=logging.INFO,
            datefmt='%d/%m/%Y %H:%M:%S',
            format='%(asctime)s -- %(levelname)s -- %(message)s'
            )
    
    # Execute right mode
    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    __cli()
