"""
Test FireBurner class
"""

import pytest
import logging
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

@pytest.mark.parametrize("smarts,inchi,expected_inchi",
                         [(reaction_smarts, substate_inchi, ['InChI=1S/C3H7O6P/c1-2(3(4)5)9-10(6,7)8/h2H,1H3,(H,4,5)(H2,6,7,8)', 'InChI=1S/C10H15N5O11P2/c11-10-13-7-4(8(18)14-10)12-2-15(7)9-6(17)5(16)3(25-9)1-24-28(22,23)26-27(19,20)21/h2-3,5-6,9,16-17H,1H2,(H,22,23)(H2,19,20,21)(H3,11,13,14,18)'])
                          ])
def test_compute(smarts, inchi, expected_inchi):
    # Things that are expected to return results
    rb = RuleBurner(rsmarts_list=[reaction_smarts], inchi_list=[substate_inchi])
    result = [x for x in rb.compute()][0]
    assert all(x in result for x in ['rule_id', 'substrate_id', 'product_inchikeys', 'product_inchis', 'product_smiles']), result
    assert all(x in expected_inchi for x in result['product_inchis']), result


def test_compute_wrong_reaction(caplog):
    # Wrong reaction depiction should be caught and logged by RuleBurner
    with caplog.at_level(logging.ERROR):
        rb = RuleBurner(rsmarts_list=['DUMMY'], inchi_list=[substate_inchi])
        result = [x for x in rb.compute()]
    assert "RULE-CONVERSION-ERROR" in caplog.text, (result, caplog.text)


def test_compute_wrong_chemical(caplog):
    # Wrong chemical depiction should be caught and logged by RuleBurner
    with caplog.at_level(logging.ERROR):
        rb = RuleBurner(rsmarts_list=[reaction_smarts], inchi_list=['DUMMY'])
        result = [x for x in rb.compute()]
    assert "CHEM-CONVERSION-ERROR" in caplog.text, (result, caplog.text)


def test_compute_timeout(caplog):
    # Timeout should be logged
    rsmarts = "([#6&v4:1](=[#6&v4:2](-[#6&v4:3])-[#6&v4:4](-[#6&v4:5](-[#6&v4:6](=[#6&v4:7](-[#6&v4:8])-[#6&v4:9](-[#6&v4:10](-[#6&v4:11](=[#6&v4:12](-[#6&v4:13])-[#6&v4:14](-[#6&v4:15](-[#6&v4:16](=[#6&v4:17](-[#6&v4:18])-[#6&v4:19](-[#6&v4:20](-[#6&v4:21](=[#6&v4:22](-[#6&v4:23])-[#6&v4:24](-[#6&v4:25](-[#6&v4:26](=[#6&v4:27](-[#6&v4:28])-[#6&v4:29](-[#6&v4:30](-[#6&v4:31](=[#6&v4:32])-[#1&v1:33])(-[#1&v1:34])-[#1&v1:35])(-[#1&v1:36])-[#1&v1:37])-[#1&v1:38])(-[#1&v1:39])-[#1&v1:40])(-[#1&v1:41])-[#1&v1:42])-[#1&v1:43])(-[#1&v1:44])-[#1&v1:45])(-[#1&v1:46])-[#1&v1:47])-[#1&v1:48])(-[#1&v1:49])-[#1&v1:50])(-[#1&v1:51])-[#1&v1:52])-[#1&v1:53])(-[#1&v1:54])-[#1&v1:55])(-[#1&v1:56])-[#1&v1:57])-[#1&v1:58])(-[#1&v1:59])-[#1&v1:60])(-[#1&v1:61])-[#1&v1:62])(-[#6&v4:63](-[#6&v4:64](-[#6&v4:65](=[#6&v4:66](-[#6&v4:67](-[#6&v4:68](-[#6&v4:69](=[#6&v4:70](-[#6&v4:71](-[#6&v4:72](-[#6&v4:73](=[#6&v4:74](-[#6&v4:75](-[#6&v4:76](-[#6&v4:77](=[#6&v4:78](-[#6&v4:79](-[#6&v4:80](-[#6&v4:81](=[#6&v4:82](-[#6&v4:83](-[#8&v2:84]-[#15&v5:85])(-[#1&v1:86])-[#1&v1:87])-[#1&v1:88])-[#6&v4:89])(-[#1&v1:90])-[#1&v1:91])(-[#1&v1:92])-[#1&v1:93])-[#1&v1:94])-[#6&v4:95])(-[#1&v1:96])-[#1&v1:97])(-[#1&v1:98])-[#1&v1:99])-[#1&v1:100])-[#6&v4:101])(-[#1&v1:102])-[#1&v1:103])(-[#1&v1:104])-[#1&v1:105])-[#1&v1:106])-[#6&v4:107])(-[#1&v1:108])-[#1&v1:109])(-[#1&v1:110])-[#1&v1:111])-[#1&v1:112])-[#6&v4:113])(-[#1&v1:114])-[#1&v1:115])(-[#1&v1:116])-[#1&v1:117])-[#1&v1:118])>>([#15&v5:85]-[#8&v2:84]-[#6&v4:83](-[#6&v4:29](-[#6&v4:27](=[#6&v4:26](-[#1&v1:38])-[#1&v1])-[#6&v4:28])(-[#1&v1:36])-[#1&v1:37])(-[#1&v1:86])-[#1&v1:87].[#6&v4:30](-[#6&v4:31](=[#6&v4:32])-[#1&v1:33])(-[#8&v2]-[#15&v5](-[#8&v2]-[#1&v1])(-[#8&v2]-[#15&v5](-[#8&v2]-[#1&v1])(=[#8&v2])-[#8&v2]-[#1&v1])=[#8&v2])(-[#1&v1:34])-[#1&v1:35].[#6&v4:82](=[#6&v4:81](-[#6&v4:89])-[#6&v4:80](-[#6&v4:79](-[#8&v2]-[#15&v5](-[#8&v2]-[#1&v1])(-[#8&v2]-[#15&v5](-[#8&v2]-[#1&v1])(=[#8&v2])-[#8&v2]-[#1&v1])=[#8&v2])(-[#1&v1:92])-[#1&v1:93])(-[#1&v1:90])-[#1&v1:91])(-[#1&v1:88])-[#1&v1].[#6&v4:78](=[#6&v4:77](-[#6&v4:95])-[#6&v4:76](-[#6&v4:75](-[#8&v2]-[#15&v5](-[#8&v2]-[#1&v1])(-[#8&v2]-[#15&v5](-[#8&v2]-[#1&v1])(=[#8&v2])-[#8&v2]-[#1&v1])=[#8&v2])(-[#1&v1:98])-[#1&v1:99])(-[#1&v1:96])-[#1&v1:97])(-[#1&v1:94])-[#1&v1].[#6&v4:74](=[#6&v4:73](-[#6&v4:101])-[#6&v4:72](-[#6&v4:71](-[#8&v2]-[#15&v5](-[#8&v2]-[#1&v1])(-[#8&v2]-[#15&v5](-[#8&v2]-[#1&v1])(=[#8&v2])-[#8&v2]-[#1&v1])=[#8&v2])(-[#1&v1:104])-[#1&v1:105])(-[#1&v1:102])-[#1&v1:103])(-[#1&v1:100])-[#1&v1].[#6&v4:70](=[#6&v4:69](-[#6&v4:107])-[#6&v4:68](-[#6&v4:67](-[#8&v2]-[#15&v5](-[#8&v2]-[#1&v1])(-[#8&v2]-[#15&v5](-[#8&v2]-[#1&v1])(=[#8&v2])-[#8&v2]-[#1&v1])=[#8&v2])(-[#1&v1:110])-[#1&v1:111])(-[#1&v1:108])-[#1&v1:109])(-[#1&v1:106])-[#1&v1].[#6&v4:66](=[#6&v4:65](-[#6&v4:113])-[#6&v4:64](-[#6&v4:63](-[#8&v2]-[#15&v5](-[#8&v2]-[#1&v1])(-[#8&v2]-[#15&v5](-[#8&v2]-[#1&v1])(=[#8&v2])-[#8&v2]-[#1&v1])=[#8&v2])(-[#1&v1:116])-[#1&v1:117])(-[#1&v1:114])-[#1&v1:115])(-[#1&v1:112])-[#1&v1].[#6&v4:1](=[#6&v4:2](-[#6&v4:3])-[#6&v4:4](-[#6&v4:5](-[#8&v2]-[#15&v5](-[#8&v2]-[#1&v1])(-[#8&v2]-[#15&v5](-[#8&v2]-[#1&v1])(=[#8&v2])-[#8&v2]-[#1&v1])=[#8&v2])(-[#1&v1:59])-[#1&v1:60])(-[#1&v1:61])-[#1&v1:62])(-[#1&v1:118])-[#1&v1].[#6&v4:6](=[#6&v4:7](-[#6&v4:8])-[#6&v4:9](-[#6&v4:10](-[#8&v2]-[#15&v5](-[#8&v2]-[#1&v1])(-[#8&v2]-[#15&v5](-[#8&v2]-[#1&v1])(=[#8&v2])-[#8&v2]-[#1&v1])=[#8&v2])(-[#1&v1:54])-[#1&v1:55])(-[#1&v1:56])-[#1&v1:57])(-[#1&v1:58])-[#1&v1].[#6&v4:11](=[#6&v4:12](-[#6&v4:13])-[#6&v4:14](-[#6&v4:15](-[#8&v2]-[#15&v5](-[#8&v2]-[#1&v1])(-[#8&v2]-[#15&v5](-[#8&v2]-[#1&v1])(=[#8&v2])-[#8&v2]-[#1&v1])=[#8&v2])(-[#1&v1:49])-[#1&v1:50])(-[#1&v1:51])-[#1&v1:52])(-[#1&v1:53])-[#1&v1].[#6&v4:16](=[#6&v4:17](-[#6&v4:18])-[#6&v4:19](-[#6&v4:20](-[#8&v2]-[#15&v5](-[#8&v2]-[#1&v1])(-[#8&v2]-[#15&v5](-[#8&v2]-[#1&v1])(=[#8&v2])-[#8&v2]-[#1&v1])=[#8&v2])(-[#1&v1:44])-[#1&v1:45])(-[#1&v1:46])-[#1&v1:47])(-[#1&v1:48])-[#1&v1].[#6&v4:21](=[#6&v4:22](-[#6&v4:23])-[#6&v4:24](-[#6&v4:25](-[#8&v2]-[#15&v5](-[#8&v2]-[#1&v1])(-[#8&v2]-[#15&v5](-[#8&v2]-[#1&v1])(=[#8&v2])-[#8&v2]-[#1&v1])=[#8&v2])(-[#1&v1:39])-[#1&v1:40])(-[#1&v1:41])-[#1&v1:42])(-[#1&v1:43])-[#1&v1])"
    inchi = "InChI=1S/C70H116O7P2/c1-57(2)29-16-30-58(3)31-17-32-59(4)33-18-34-60(5)35-19-36-61(6)37-20-38-62(7)39-21-40-63(8)41-22-42-64(9)43-23-44-65(10)45-24-46-66(11)47-25-48-67(12)49-26-50-68(13)51-27-52-69(14)53-28-54-70(15)55-56-76-79(74,75)77-78(71,72)73/h29,31,33,35,37,39,41,43,45,47,49,51,53,55H,16-28,30,32,34,36,38,40,42,44,46,48,50,52,54,56H2,1-15H3,(H,74,75)(H2,71,72,73)"
    with caplog.at_level(logging.WARNING):
        rb = RuleBurner(rsmarts_list=[rsmarts], inchi_list=[inchi], with_hs=True)
        result = [x for x in rb.compute(timeout=0.1)]
    assert "timed-out" in caplog.text, (result, caplog.text)
