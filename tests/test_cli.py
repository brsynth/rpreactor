import pytest
from rpreactor import cli
import json


@pytest.fixture
def result_example():
    result_input = [{
        'rule_id': '0',
        'substrate_id': '0',
        'product_list':
            [['RDKIT_DUMMY', 'RDKIT_DUMMY']],
        'product_inchikeys':
            [['CSZRNWHGZPKNKY-UHFFFAOYSA-N', 'QGWNDRXFNXRZMB-UHFFFAOYSA-N']],
        'product_inchis':
            [['InChI=1S/C3H7O6P/c1-2(3(4)5)9-10(6,7)8/h2H,1H3,(H,4,5)(H2,6,7,8)', 'InChI=1S/C10H15N5O11P2/c11-10-13-7-4(8(18)14-10)12-2-15(7)9-6(17)5(16)3(25-9)1-24-28(22,23)26-27(19,20)21/h2-3,5-6,9,16-17H,1H2,(H,22,23)(H2,19,20,21)(H3,11,13,14,18)']],
        'product_smiles':
            [['[H]OC(=O)C([H])(OP(=O)(O[H])O[H])C([H])([H])[H]', '[H]N=c1nc(O[H])c2nc([H])n(C3([H])OC([H])(C([H])([H])OP(=O)(O[H])OP(=O)(O[H])O[H])C([H])(O[H])C3([H])O[H])c2n1[H]']]
    }]
    result_json = """
        [
         {
          "rule_id": "0",
          "substrate_id": "0",
          "product_inchikeys": [
           [
            "CSZRNWHGZPKNKY-UHFFFAOYSA-N",
            "QGWNDRXFNXRZMB-UHFFFAOYSA-N"
           ]
          ],
          "product_inchis": [
           [
            "InChI=1S/C3H7O6P/c1-2(3(4)5)9-10(6,7)8/h2H,1H3,(H,4,5)(H2,6,7,8)",
            "InChI=1S/C10H15N5O11P2/c11-10-13-7-4(8(18)14-10)12-2-15(7)9-6(17)5(16)3(25-9)1-24-28(22,23)26-27(19,20)21/h2-3,5-6,9,16-17H,1H2,(H,22,23)(H2,19,20,21)(H3,11,13,14,18)"
           ]
          ],
          "product_smiles": [
           [
            "[H]OC(=O)C([H])(OP(=O)(O[H])O[H])C([H])([H])[H]",
            "[H]N=c1nc(O[H])c2nc([H])n(C3([H])OC([H])(C([H])([H])OP(=O)(O[H])OP(=O)(O[H])O[H])C([H])(O[H])C3([H])O[H])c2n1[H]"
           ]
          ]
         }
        ]
    """
    return result_input, result_json


def test_jsonify(result_example):
    results_input, json_output = result_example
    assert json.loads(cli.jsonify(results_input)) == json.loads(json_output)


def test_write(tmp_path):
    text = 'Hello, world'
    ofile = tmp_path / 'ofile.json'
    cli.write(text, ofile=ofile)
    with open(ofile) as fh:
        assert fh.read() == text


def test_build_arg_parser(mocker):
    # No arguments
    args = ['prog']
    mocker.patch('sys.argv', args)
    parser = cli.__build_arg_parser()
    args = parser.parse_args()
    with pytest.raises(AttributeError):
        args.func()
    # inline mode
    args = ['prog', 'inline', '--rsmarts', 'DUMMY_RSMARTS', '--inchi', 'DUMMY_INCHI']
    mocker.patch('sys.argv', args)
    parser = cli.__build_arg_parser()
    parser.parse_args()
    # infile mode
    args = ['prog', 'infile', '--rfile', 'DUMMY_RFILE', '--cfile', 'DUMMY_CFILE']
    mocker.patch('sys.argv', args)
    parser = cli.__build_arg_parser()
    parser.parse_args()


def test_inline_mode(mocker):
    substrate_inchi = 'InChI=1S/C3H6O3/c1-2(4)3(5)6/h2,4H,1H3,(H,5,6)'
    reaction_smarts = '([#8&v2:1](-[#6&v4:2](-[#6&v4:3](-[#8&v2:4]-[#1&v1:5])=[#8&v2:6])(-[#6&v4:7](-[#1&v1:8])(-[#1&v1:9])-[#1&v1:10])-[#1&v1:11])-[#1&v1:12])>>([#15&v5](=[#8&v2])(-[#8&v2]-[#1&v1])(-[#8&v2]-[#1&v1])-[#8&v2:1]-[#6&v4:2](-[#6&v4:3](-[#8&v2:4]-[#1&v1:5])=[#8&v2:6])(-[#6&v4:7](-[#1&v1:8])(-[#1&v1:9])-[#1&v1:10])-[#1&v1:11].[#7&v3](=[#6&v4]1:[#7&v3]:[#6&v4](-[#8&v2]-[#1&v1]):[#6&v4]2:[#7&v3]:[#6&v4](-[#1&v1]):[#7&v3](-[#6&v4]3(-[#1&v1])-[#8&v2]-[#6&v4](-[#6&v4](-[#8&v2]-[#15&v5](=[#8&v2])(-[#8&v2]-[#1&v1])-[#8&v2]-[#15&v5](-[#8&v2]-[#1&v1:12])(=[#8&v2])-[#8&v2]-[#1&v1])(-[#1&v1])-[#1&v1])(-[#1&v1])-[#6&v4](-[#8&v2]-[#1&v1])(-[#1&v1])-[#6&v4]-3(-[#8&v2]-[#1&v1])-[#1&v1]):[#6&v4]:2:[#7&v3]:1-[#1&v1])-[#1&v1])'
    # Without ID
    args = ['prog', 'inline', '--rsmarts', reaction_smarts, '--inchi', substrate_inchi]
    mocker.patch('sys.argv', args)
    parser = cli.__build_arg_parser()
    args = parser.parse_args()
    cli.__inline_mode(args)
    # With IDs
    args = ['prog', 'inline', '--rsmarts', reaction_smarts, '--rid', 'rid', '--inchi', substrate_inchi, '--cid', 'cid']
    mocker.patch('sys.argv', args)
    parser = cli.__build_arg_parser()
    args = parser.parse_args()
    cli.__inline_mode(args)


def test_infile_mode(mocker, tmpdir):
    substrate_inchi = 'InChI=1S/C3H6O3/c1-2(4)3(5)6/h2,4H,1H3,(H,5,6)'
    reaction_smarts = '([#8&v2:1](-[#6&v4:2](-[#6&v4:3](-[#8&v2:4]-[#1&v1:5])=[#8&v2:6])(-[#6&v4:7](-[#1&v1:8])(-[#1&v1:9])-[#1&v1:10])-[#1&v1:11])-[#1&v1:12])>>([#15&v5](=[#8&v2])(-[#8&v2]-[#1&v1])(-[#8&v2]-[#1&v1])-[#8&v2:1]-[#6&v4:2](-[#6&v4:3](-[#8&v2:4]-[#1&v1:5])=[#8&v2:6])(-[#6&v4:7](-[#1&v1:8])(-[#1&v1:9])-[#1&v1:10])-[#1&v1:11].[#7&v3](=[#6&v4]1:[#7&v3]:[#6&v4](-[#8&v2]-[#1&v1]):[#6&v4]2:[#7&v3]:[#6&v4](-[#1&v1]):[#7&v3](-[#6&v4]3(-[#1&v1])-[#8&v2]-[#6&v4](-[#6&v4](-[#8&v2]-[#15&v5](=[#8&v2])(-[#8&v2]-[#1&v1])-[#8&v2]-[#15&v5](-[#8&v2]-[#1&v1:12])(=[#8&v2])-[#8&v2]-[#1&v1])(-[#1&v1])-[#1&v1])(-[#1&v1])-[#6&v4](-[#8&v2]-[#1&v1])(-[#1&v1])-[#6&v4]-3(-[#8&v2]-[#1&v1])-[#1&v1]):[#6&v4]:2:[#7&v3]:1-[#1&v1])-[#1&v1])'
    with open(f'{tmpdir}/cfile.tsv', 'w') as ofh:
        ofh.writelines([
            'inchi\tchem_id\n',
            f'{substrate_inchi}\ncid\n'
        ])
    with open(tmpdir/'rfile.tsv', 'w') as ofh:
        ofh.writelines([
            'rule_smarts\trule_id\n',
            f'{reaction_smarts}\nrid\n'
        ])
    args = ['prog', 'infile', '--rfile', f'{tmpdir}/rfile.tsv', '--cfile', f'{tmpdir}/cfile.tsv']
    mocker.patch('sys.argv', args)
    parser = cli.__build_arg_parser()
    args = parser.parse_args()
    cli.__infile_mode(args)
