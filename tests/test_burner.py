"""
Test FireBurner class
"""

import pytest
from rdkit import Chem
from rpreactor.rule.burner import RuleBurner


# Data for tests
substate_inchi = 'InChI=1S/C3H6O3/c1-2(4)3(5)6/h2,4H,1H3,(H,5,6)'
reaction_smarts = '([#8&v2:1](-[#6&v4:2](-[#6&v4:3](-[#8&v2:4]-[#1&v1:5])=[#8&v2:6])(-[#6&v4:7](-[#1&v1:8])(-[#1&v1:9])-[#1&v1:10])-[#1&v1:11])-[#1&v1:12])>>([#15&v5](=[#8&v2])(-[#8&v2]-[#1&v1])(-[#8&v2]-[#1&v1])-[#8&v2:1]-[#6&v4:2](-[#6&v4:3](-[#8&v2:4]-[#1&v1:5])=[#8&v2:6])(-[#6&v4:7](-[#1&v1:8])(-[#1&v1:9])-[#1&v1:10])-[#1&v1:11].[#7&v3](=[#6&v4]1:[#7&v3]:[#6&v4](-[#8&v2]-[#1&v1]):[#6&v4]2:[#7&v3]:[#6&v4](-[#1&v1]):[#7&v3](-[#6&v4]3(-[#1&v1])-[#8&v2]-[#6&v4](-[#6&v4](-[#8&v2]-[#15&v5](=[#8&v2])(-[#8&v2]-[#1&v1])-[#8&v2]-[#15&v5](-[#8&v2]-[#1&v1:12])(=[#8&v2])-[#8&v2]-[#1&v1])(-[#1&v1])-[#1&v1])(-[#1&v1])-[#6&v4](-[#8&v2]-[#1&v1])(-[#1&v1])-[#6&v4]-3(-[#8&v2]-[#1&v1])-[#1&v1]):[#6&v4]:2:[#7&v3]:1-[#1&v1])-[#1&v1])'
tuple_product_inchikeys = ('CSZRNWHGZPKNKY-UHFFFAOYSA-N', 'QGWNDRXFNXRZMB-UHFFFAOYSA-N')
tuple_product_smiles = ('[H][O][C](=[O])[C]([H])([O][P](=[O])([O][H])[O][H])[C]([H])([H])[H]', '[H][N]=[c]1[n][c]([O][H])[c]2[n][c]([H])[n]([C]3([H])[O][C]([H])([C]([H])([H])[O][P](=[O])([O][H])[O][P](=[O])([O][H])[O][H])[C]([H])([O][H])[C]3([H])[O][H])[c]2[n]1[H]')
tuple_product_inchis = ('InChI=1S/C3H7O6P/c1-2(3(4)5)9-10(6,7)8/h2H,1H3,(H,4,5)(H2,6,7,8)', 'InChI=1S/C10H15N5O11P2/c11-10-13-7-4(8(18)14-10)12-2-15(7)9-6(17)5(16)3(25-9)1-24-28(22,23)26-27(19,20)21/h2-3,5-6,9,16-17H,1H2,(H,22,23)(H2,19,20,21)(H3,11,13,14,18)')


def dummy_worker(**kwargs):
    import time
    time.sleep(1)


def test_init():
    # Empty is OK
    rb = RuleBurner(rsmarts_list=[], inchi_list=[])  # Empty is OK
    rb.compute()


def test_standardize_chemical_1():
    rb = RuleBurner(rsmarts_list=[], inchi_list=[], with_hs=False)
    rdmol = Chem.MolFromSmiles('[H][O][C](=[O])[C]([H])([O][H])[C]([H])([H])[H]')
    rdmol_std_1 = rb._standardize_chemical(rdmol)
    assert Chem.MolToSmiles(rdmol_std_1) == 'CC(O)C(=O)O'
    # With H
    rb = RuleBurner(rsmarts_list=[], inchi_list=[], with_hs=True)
    rdmol_std_2 = rb._standardize_chemical(rdmol)
    assert Chem.MolToSmiles(rdmol_std_2, allHsExplicit=True) == '[H][O][C](=[O])[C]([H])([O][H])[C]([H])([H])[H]'


def test_standardize_chemical_2():
    # Data
    violacein_smiles = 'OC1=NC(=C\\C1=C1/C(O)=NC2=CC=CC=C12)C1=CNC2=C1C=C(O)C=C2'
    violacein_mol = Chem.MolFromSmiles(violacein_smiles, sanitize=False)
    # Test simplest case
    with pytest.raises(NotImplementedError):
        rb = RuleBurner(rsmarts_list=[], inchi_list=[], with_hs=False, with_stereo=True)
        std_mol_1 = rb._standardize_chemical(violacein_mol)
        assert Chem.MolToSmiles(std_mol_1) == 'OC1=NC(c2c[nH]c3ccc(O)cc23)=C/C1=C1\\C(O)=Nc2ccccc21'
    # Test adding Hs
    with pytest.raises(NotImplementedError):
        rb = RuleBurner(rsmarts_list=[], inchi_list=[], with_hs=True, with_stereo=True)
        std_mol_2 = rb._standardize_chemical(violacein_mol)
        assert Chem.MolToSmiles(std_mol_2) == '[H]OC1=NC(c2c([H])n([H])c3c([H])c([H])c(O[H])c([H])c23)=C([H])/C1=C1\\C(O[H])=Nc2c([H])c([H])c([H])c([H])c21'
    # Test removing stereo
    rb = RuleBurner(rsmarts_list=[], inchi_list=[], with_hs=False, with_stereo=False)
    std_mol_3 = rb._standardize_chemical(violacein_mol)
    assert Chem.MolToSmiles(std_mol_3) == 'O=C1NC(c2c[nH]c3ccc(O)cc23)=CC1=C1C(=O)Nc2ccccc21'
    # Test adding Hs + removing stereo
    rb = RuleBurner(rsmarts_list=[], inchi_list=[], with_hs=True, with_stereo=False)
    std_mol_4 = rb._standardize_chemical(violacein_mol)
    assert Chem.MolToSmiles(std_mol_4) == '[H]Oc1c([H])c([H])c2c(c1[H])c(C1=C([H])C(=C3C(=O)N([H])c4c([H])c([H])c([H])c([H])c43)C(=O)N1[H])c([H])n2[H]'


def test_standardize_chemical_3():
    # Data
    rb = RuleBurner(rsmarts_list=[], inchi_list=[])
    wrong_smiles = '[H]OC(=O)C([H])([H])C([H])([H])C([H])(N=C(O[H])C([H])([H])C([H])([H])C([H])(N=C(O[H])C([H])(OP(=O)(O[H])OC([H])([H])C([H])(O[H])C([H])(O[H])C([H])(O[H])C([H])([H])n1c2nc(=O)nc(O[H])c-2c([H])c2c([H])c([H])c(OP(=O)(OC([H])([H])C(C([H])([H])[H])(C([H])([H])[H])C([H])(O[H])C(=NC([H])([H])C([H])([H])C(=NC([H])([H])C([H])([H])SC(=O)C([H])([H])C([H])([H])C([H])([H])C([H])(C(=C([H])[H])C([H])([H])[H])C([H])([H])C(=O)O[H])O[H])O[H])OP(=O)(O[H])OC([H])([H])C3([H])OC([H])(n4[c]([H])n([H])[c]5[c](N([H])[H])[n][c]([H])[n][c]54)C([H])(O[H])C3([H])OP(=O)(O[H])O[H])c([H])c21)C([H])([H])[H])C(=O)O[H])C(=O)O[H]'
    # Test
    wrong_mol = Chem.MolFromSmiles(wrong_smiles, sanitize=False)
    with pytest.raises(Exception):
        rb._standardize_chemical(wrong_mol)


def test_standardize_results_1():
    tuple_tuple_raw = ((
                Chem.MolFromSmiles('[H][O][C](=[O])[C]([H])([O][P](=[O])([O][H])[O][H])[C]([H])([H])[H]'),
                Chem.MolFromSmiles('[H][N]=[c]1[n][c]([O][H])[c]2[n][c]([H])[n]([C]3([H])[O][C]([H])([C]([H])([H])[O][P](=[O])([O][H])[O][P](=[O])([O][H])[O][H])[C]([H])([O][H])[C]3([H])[O][H])[c]2[n]1[H]')
            ),(
                Chem.MolFromInchi('InChI=1S/C5H6N5O/c6-5-9-3-2(4(11)10-5)7-1-8-3/h1H,9H2,(H,7,8)(H2,6,10,11)')
            ))
    rb = RuleBurner(rsmarts_list=[], inchi_list=[], with_hs=True, with_stereo=False)
    tuple_tuple_rdmol, tuple_index_failed = rb._standardize_results(tuple_tuple_raw)
    assert len(tuple_tuple_rdmol) == 1
    assert tuple_index_failed == [1]


def test_jsonify():
    rb = RuleBurner(rsmarts_list=[], inchi_list=[])
    assert rb._jsonify(rid='RID', cid='CID').replace('\n', '') == """{ "rule_id": "RID", "substrate_id": "CID", "fire_timed_out": null, "fire_exec_time": null}"""


def test_handle_result():
    tuple_raw = (
            Chem.MolFromSmiles('[H][O][C](=[O])[C]([H])([O][P](=[O])([O][H])[O][H])[C]([H])([H])[H]'),
            Chem.MolFromSmiles('[H][N]=[c]1[n][c]([O][H])[c]2[n][c]([H])[n]([C]3([H])[O][C]([H])([C]([H])([H])[O][P](=[O])([O][H])[O][P](=[O])([O][H])[O][H])[C]([H])([O][H])[C]3([H])[O][H])[c]2[n]1[H]')
            )
    rb = RuleBurner(rsmarts_list=[], inchi_list=[], with_hs=True, with_stereo=False)
    tuple_tuple_rdmol, tuple_tuple_failed = rb._standardize_results(tuple_tuple_rdmol=(tuple_raw,), )
    inchikeys, inchis, smiles = rb._handle_results(list_list_rdmol=tuple_tuple_rdmol)
    # Check number products
    assert len(inchikeys) == len(inchis) == len(smiles) == 1   # Only one set of result
    assert len(inchikeys[0]) == len(inchis[0]) == len(smiles[0]) == 2  # 2 products
    # Check Inchikeys
    assert inchikeys[0][0] == 'CSZRNWHGZPKNKY-UHFFFAOYSA-N'
    assert inchikeys[0][1] == 'QGWNDRXFNXRZMB-UHFFFAOYSA-N'
    # Check Inchis
    assert inchis[0][0] == 'InChI=1S/C3H7O6P/c1-2(3(4)5)9-10(6,7)8/h2H,1H3,(H,4,5)(H2,6,7,8)'
    assert inchis[0][1] == 'InChI=1S/C10H15N5O11P2/c11-10-13-7-4(8(18)14-10)12-2-15(7)9-6(17)5(16)3(25-9)1-24-28(22,23)26-27(19,20)21/h2-3,5-6,9,16-17H,1H2,(H,22,23)(H2,19,20,21)(H3,11,13,14,18)'
    # Check SMILES #1
    assert smiles[0][0] == '[H]OC(=O)C([H])(OP(=O)(O[H])O[H])C([H])([H])[H]'
    rdmol = Chem.MolFromSmiles(smiles[0][0])
    rdmol = Chem.AddHs(rdmol)
    assert Chem.MolToSmiles(rdmol, allHsExplicit=True) == '[H][O][C](=[O])[C]([H])([O][P](=[O])([O][H])[O][H])[C]([H])([H])[H]'
    # Check SMILES #2
    assert smiles[0][1] == '[H]N=c1nc(O[H])c2nc([H])n(C3([H])OC([H])(C([H])([H])OP(=O)(O[H])OP(=O)(O[H])O[H])C([H])(O[H])C3([H])O[H])c2n1[H]'
    rdmol = Chem.MolFromSmiles(smiles[0][1])
    rdmol = Chem.AddHs(rdmol)
    assert Chem.MolToSmiles(rdmol, allHsExplicit=True) == '[H][N]=[c]1[n][c]([O][H])[c]2[n][c]([H])[n]([C]3([H])[O][C]([H])([C]([H])([H])[O][P](=[O])([O][H])[O][P](=[O])([O][H])[O][H])[C]([H])([O][H])[C]3([H])[O][H])[c]2[n]1[H]'


def test_compute():
    # Wrong reaction depiction should be caught and logged by RuleBurner
    rb = RuleBurner(rsmarts_list=['DUMMY'], inchi_list=[substate_inchi])
    rb.compute()
    # Wrong chemical depiction should be caught and logged by RuleBurner
    rb = RuleBurner(rsmarts_list=[reaction_smarts], inchi_list=['DUMMY'])
    rb.compute()
    # Timeout should be logged
    rb = RuleBurner(rsmarts_list=[reaction_smarts], inchi_list=[substate_inchi])
    rb.compute(timeout=0)
    assert ''.join(rb._json).find('"fire_timed_out": true')
    # OK
    rb = RuleBurner(rsmarts_list=[reaction_smarts], inchi_list=[substate_inchi])
    rb.compute()
    assert ''.join(rb._json).find('InChI=1S/C3H7O6P/c1-2(3(4)5)9-10(6,7)8/h2H,1H3,(H,4,5)(H2,6,7,8)')
    assert ''.join(rb._json).find('InChI=1S/C10H15N5O11P2/c11-10-13-7-4(8(18)14-10)12-2-15(7)9-6(17)5(16)3(25-9)1-24-28(22,23)26-27(19,20)21/h2-3,5-6,9,16-17H,1H2,(H,22,23)(H2,19,20,21)(H3,11,13,14,18)')

