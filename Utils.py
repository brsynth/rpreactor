"""
Set of methods to handle reaction I/Os
"""


import copy
import rdkit
import logging

from rdkit import Chem
from rdkit import RDLogger


RD_LOGGER = RDLogger.logger()
RD_LOGGER.setLevel(RDLogger.CRITICAL)  # Silent most of RDKit complains


class ChemConversionError(Exception):
    """Raised when something went wrong during chemical conversion to RDKit mol object."""
    
    def __init__(self, msg):
        self._msg = msg

    def __str__(self):
        return "CHEM-CONVERSION-ERROR: {}".format(self._msg)


def wild_stereo_removal(rdmol):
    """Wild stereo removal using back and forth Inchi depiction.
    
    :param      rdmol:      RDKit mol
    :returns    rdmol_new:  newly generated RDKit mol 
    """
    tmp_rdmol = copy.deepcopy(rdmol)
    Chem.RemoveStereochemistry(tmp_rdmol)
    return Chem.MolFromInchi(Chem.MolToInchi(tmp_rdmol))


def standardize_chemical(rdmol, add_hs=True, rm_stereo=True):
    """Standardize a chemical using RDKit sanitize method.

    :param      rdmol:      RDKit mol object
    :param      add_hs:     append Hs, bool (default: True)
    :param      rm_stereo:  remove stereo, bool (default: True)
    :returns    rdmol:      RDKit mol object
    """
    try:
        Chem.SanitizeMol(rdmol)
        if rm_stereo:  # Important: do this before adding Hs (else re-add Hs)
            rdmol = wild_stereo_removal(rdmol)
        if add_hs:
            rdmol = Chem.AddHs(rdmol)   
        return rdmol
    except Exception as e:
        logging.warning(e)
        raise e


def standardize_results(tuple_tuple_rdmol, add_hs=True, rm_stereo=True):
    """Perform sanitization and remove duplicates from reaction rule results.
    
    :param      tuple_tuple_rdmol:   tuple of tuple of RDKit Mol
    :param      add_hs:              append Hs, bool (default: True)
    :param      rm_stereo:           remove stereo, bool (default: True)
    :returns    list_list_std:       list of list of standardized RDKit Mol
    """
    uniq_depics = set()
    list_list_std = list()

    for tuple_rdmol in tuple_tuple_rdmol:
        try: 
            list_std = list()
            list_inchikeys = list()
            # Standardize
            for rdmol in tuple_rdmol:
                for rd_frag in Chem.GetMolFrags(rdmol, asMols=True, sanitizeFrags=False):
                    list_std.append(standardize_chemical(rd_frag, add_hs=add_hs, rm_stereo=rm_stereo))
            # Get Inchikeys
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
            # Stoer only if unique depiction never met
            if depic not in uniq_depics:
                uniq_depics.add(depic)
                list_list_std.append(list_std)
        except ChemConversionError as e:
            logging.warning("{}".format(e))
            raise e
        except Exception as e:
            logging.warning("Cannot handle a tuple of result, skipped")
            logging.warning("{}".format(e))
    return list_list_std


def handle_results(list_list_rdmol):
    """Generate InchiKey, Inchi and SMILES from results.
    
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
            list_std = list()
            for rdmol in list_rdmol:
                # Get & check depictions
                inchikey = Chem.MolToInchiKey(rdmol)  # DEBUG: this part could be optimized
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
