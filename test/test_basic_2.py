import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import pytest

from reactor.cli import RuleBurner
from reactor.Utils import standardize_chemical


class TestBasic2(object):


    def test_standardize_chemical_1(self):
        # Data
        violacein_smiles = 'OC1=NC(=C\\C1=C1/C(O)=NC2=CC=CC=C12)C1=CNC2=C1C=C(O)C=C2'
        violacein_mol = Chem.MolFromSmiles(violacein_smiles, sanitize=False)
        # Test simplest case
        std_mol_1 = standardize_chemical(violacein_mol, add_hs=False, rm_stereo=False)
        assert Chem.MolToSmiles(std_mol_1) == 'OC1=NC(c2c[nH]c3ccc(O)cc23)=C/C1=C1\\C(O)=Nc2ccccc21'
        # Test adding Hs
        std_mol_2 = standardize_chemical(violacein_mol, add_hs=True, rm_stereo=False)
        assert Chem.MolToSmiles(std_mol_2) == '[H]OC1=NC(c2c([H])n([H])c3c([H])c([H])c(O[H])c([H])c23)=C([H])/C1=C1\\C(O[H])=Nc2c([H])c([H])c([H])c([H])c21'
        # Test removing stereo
        std_mol_3 = standardize_chemical(violacein_mol, add_hs=False, rm_stereo=True)
        assert Chem.MolToSmiles(std_mol_3) == 'O=C1NC(c2c[nH]c3ccc(O)cc23)=CC1=C1C(=O)Nc2ccccc21'
        # Test adding Hs + removing stereo
        std_mol_4 = standardize_chemical(violacein_mol, add_hs=True, rm_stereo=True)
        assert Chem.MolToSmiles(std_mol_4) == '[H]Oc1c([H])c([H])c2c(c1[H])c(C1=C([H])C(=C3C(=O)N([H])c4c([H])c([H])c([H])c([H])c43)C(=O)N1[H])c([H])n2[H]'


    def test_standardize_chemical_2(self):
        # Data
        wrong_smiles = '[H]OC(=O)C([H])([H])C([H])([H])C([H])(N=C(O[H])C([H])([H])C([H])([H])C([H])(N=C(O[H])C([H])(OP(=O)(O[H])OC([H])([H])C([H])(O[H])C([H])(O[H])C([H])(O[H])C([H])([H])n1c2nc(=O)nc(O[H])c-2c([H])c2c([H])c([H])c(OP(=O)(OC([H])([H])C(C([H])([H])[H])(C([H])([H])[H])C([H])(O[H])C(=NC([H])([H])C([H])([H])C(=NC([H])([H])C([H])([H])SC(=O)C([H])([H])C([H])([H])C([H])([H])C([H])(C(=C([H])[H])C([H])([H])[H])C([H])([H])C(=O)O[H])O[H])O[H])OP(=O)(O[H])OC([H])([H])C3([H])OC([H])(n4[c]([H])n([H])[c]5[c](N([H])[H])[n][c]([H])[n][c]54)C([H])(O[H])C3([H])OP(=O)(O[H])O[H])c([H])c21)C([H])([H])[H])C(=O)O[H])C(=O)O[H]'
        # Test
        wrong_mol = Chem.MolFromSmiles(wrong_smiles, sanitize=False)
        with pytest.raises(Exception):
            standardize_chemical(wrong_mol)
