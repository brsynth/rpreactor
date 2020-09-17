#!/usr/bin/env python
"""
Everything to standardize chemicals (metabolites).

The idea is to use a Standardizer object to store a set of "standardization rules", and to (re)use this object
to standardize each chemical. We call "filters" those "standardization rules" to avoid confusing them with
reaction rules. Each filter is applied sequentially. For convenience, some pre-defined sequences of filters
are defined in the Standardizer class.

@author: Baudoin Delepine, 2016-2017
@author: Thomas Duigou, 2018-2019
"""

from rdkit.Chem import Cleanup, SanitizeMol, SanitizeFlags
from rdkit.Chem.AllChem import AssignStereochemistry
from .filters import Filters


class Standardizer(object):
    """Handle standardization of compound(s) through user-defined "filters".

    Some pre-defined sequences of filters are defined in this class.
    """

    def __call__(self, mol):
        """Calling a Standardizer object like a function is the same as calling its "compute" method.

        From:
            https://github.com/mcs07/MolVS/blob/master/molvs/standardize.py
        """
        return self.compute(mol)

    def __init__(self, sequence_fun=None, params=None):
        """Set up parameters for the standardization."""
        # Function to be used for standardizing compounds
        # Add you own function as method class
        if sequence_fun is None:
            self._sequence_fun = Standardizer.sequence_minimal
        elif callable(sequence_fun):     # guess: fun_filters is the function itself
            self._sequence_fun = sequence_fun
        elif type(sequence_fun) == str:  # guess: sequence_fun is the name of the method
            self._sequence_fun = getattr(Standardizer, sequence_fun)
        # Arguments to be passed to any custom standardization function
        self._params = params if params else None

    @staticmethod
    def sequence_minimal(mol):
        """Minimal standardization."""
        SanitizeMol(mol, sanitizeOps=SanitizeFlags.SANITIZE_ALL, catchErrors=False)
        AssignStereochemistry(mol, cleanIt=True, force=True, flagPossibleStereoCenters=True)  # Fix bug TD201904.01
        return mol

    @staticmethod
    def sequence_rr_legacy(mol):
        """Sequence of filters applied for the first version of RetroRules."""
        F = Filters()
        Cleanup(mol)
        SanitizeMol(mol, sanitizeOps=SanitizeFlags.SANITIZE_ALL, catchErrors=False)
        AssignStereochemistry(mol, cleanIt=True, force=True, flagPossibleStereoCenters=True)  # Fix bug TD201904.01
        mol = F.remove_isotope(mol)
        mol = F.neutralise_charge(mol)
        SanitizeMol(mol, sanitizeOps=SanitizeFlags.SANITIZE_ALL, catchErrors=False)
        mol = F.keep_biggest(mol)
        mol = F.add_hydrogen(mol, addCoords=True)
        mol = F.kekulize(mol)
        return mol

    @staticmethod
    def sequence_tunable(
            mol,
            OP_REMOVE_ISOTOPE=True, OP_NEUTRALISE_CHARGE=True,
            OP_REMOVE_STEREO=False, OP_COMMUTE_INCHI=False,
            OP_KEEP_BIGGEST=True, OP_ADD_HYDROGEN=True,
            OP_KEKULIZE=True, OP_NEUTRALISE_CHARGE_LATE=True
    ):
        """Tunable sequence of filters for standardization.

        Operations will made in the following order:
         1 RDKit Cleanup      -- always
         2 RDKIT SanitizeMol  -- always
         3 Remove isotope     -- optional (default: True)
         4 Neutralise charges -- optional (default: True)
         5 RDKit SanitizeMol  -- if 4 or 5
         6 Remove stereo      -- optional (default: False)
         7 Commute Inchi      -- if 6 or optional (default: False)
         8 Keep biggest       -- optional (default: True)
         9 RDKit SanitizeMol  -- if any (6, 7, 8)
        10 Add hydrogens      -- optional (default: True)
        11 Kekulize           -- optional (default: True)
        """
        F = Filters()
        # Always perform the basics..
        Cleanup(mol)
        SanitizeMol(mol, sanitizeOps=SanitizeFlags.SANITIZE_ALL, catchErrors=False)
        AssignStereochemistry(mol, cleanIt=True, force=True, flagPossibleStereoCenters=True)  # Fix bug TD201904.01
        #
        if OP_REMOVE_ISOTOPE:
            mol = F.remove_isotope(mol)
        if OP_NEUTRALISE_CHARGE:
            mol = F.neutralise_charge(mol)
        if any([OP_REMOVE_ISOTOPE, OP_REMOVE_ISOTOPE]):
            SanitizeMol(mol, sanitizeOps=SanitizeFlags.SANITIZE_ALL, catchErrors=False)
        #
        if OP_REMOVE_STEREO:
            mol = F.remove_stereo(mol)
            OP_COMMUTE_INCHI = True
        if OP_COMMUTE_INCHI:
            mol = F.commute_inchi(mol)
        if OP_KEEP_BIGGEST:
            mol = F.keep_biggest(mol)
        if any([OP_REMOVE_STEREO, OP_COMMUTE_INCHI, OP_KEEP_BIGGEST]):
            SanitizeMol(mol, sanitizeOps=SanitizeFlags.SANITIZE_ALL, catchErrors=False)
        #
        if OP_NEUTRALISE_CHARGE_LATE:
            mol = F.neutralise_charge(mol)
            SanitizeMol(mol, sanitizeOps=SanitizeFlags.SANITIZE_ALL, catchErrors=False)
        #
        if OP_ADD_HYDROGEN:
            mol = F.add_hydrogen(mol, addCoords=True)
        if OP_KEKULIZE:
            mol = F.kekulize(mol)
        #
        return mol

    def compute(self, mol):
        """Standardize the provided RDKit molecule."""
        if self._params is None:
            return self._sequence_fun(mol)
        else:
            return self._sequence_fun(mol, **self._params)