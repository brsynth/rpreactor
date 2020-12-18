r"""
Apply one/many rules on one/many substrates.

Thomas Duigou, INRA, 2018
"""

import os
import sys
import json
import gzip
import logging
import argparse

from rpreactor.rule.burner import RuleBurner


def jsonify(results, indent=True):
    """Export results as JSON.

    :param results: results from  RuleBurner.compute()
    :type results: list of dict
    :param indent: indent the JSON for humans if True.
    :type indent: bool
    :returns json_str: JSON string representation
    :rtype: str
    """
    # Trim-out the RDKit objects under the "product_list" key
    results_no_mol = [{k: v for k, v in d.items() if k != 'product_list'} for d in results]
    return json.dumps(results_no_mol, indent=indent)


def write(text, ofile=None, compress=False):
    """Write text to STDOUT or file.

    :param text: character string to be written
    :type text: str
    :param ofile:
    :type ofile: str
    :param compress:
    :type compress: bool
    """
    if ofile:
        ofile = os.path.abspath(ofile)
        pdir = os.path.abspath(os.path.dirname(ofile))
        os.makedirs(pdir, exist_ok=True)
        if compress:
            ofh = gzip.open(ofile, 'wb', compresslevel=9)
            ofh.write(text.encode())
        else:
            ofh = open(ofile, 'w')
            ofh.write(text)
    else:
        ofh = sys.stdout
        ofh.write(text)
    ofh.close()


def __inline_mode(args):
    """Execution mode to be used when a single rule and a single chemical
    are provided through CLI.
    """
    r = RuleBurner(with_hs=args.with_hs, with_stereo=args.with_stereo)
    if args.rid:
        r.insert_rsmarts({args.rid: args.rsmarts})
    else:
        r.insert_rsmarts([args.rsmarts])
    if args.cid:
        r.insert_inchi({args.cid: args.inchi})
    else:
        r.insert_inchi([args.inchi])
    return [x for x in r.compute(rule_mol='*', timeout=args.fire_timeout)]


def __infile_mode(args):
    """
    Execution mode to be used when rules and chemicals are provided
    in CSV files.
    """
    import csv

    with open(args.rfile, 'r') as ifh:
        reader = csv.DictReader(ifh, delimiter='\t')
        rsmarts = {row['rule_id']: row['rule_smarts'] for row in reader}

    with open(args.cfile, 'r') as ifh:
        reader = csv.DictReader(ifh, delimiter='\t')
        inchis = {row['chem_id']: row['inchi'] for row in reader}

    r = RuleBurner(with_hs=args.with_hs, with_stereo=args.with_stereo)
    r.insert_rsmarts(rsmarts)
    r.insert_inchi(inchis)
    return [x for x in r.compute(rule_mol='*', timeout=args.fire_timeout)]


def __build_arg_parser(prog='python -m rpreactor.cli'):
    desc = "Apply rules on chemicals."

    parser = argparse.ArgumentParser(description=desc, prog=prog)
    parser.add_argument('--fire_timeout', help='Rule firing timeout (seconds). Default: 60.', default=60, type=int)
    parser.add_argument('--ofile', help='Output file to store results. Default to STDOUT if none provided')
    parser.add_argument('--compress', action='store_true', help='Enable gzip compression (only when output to file).')
    parser.add_argument('--with_hs', help='Enable explicit Hs when sanitizing chemicals. Default to False.',
                        default=False, type=lambda x: (str(x).lower() == 'true'))
    parser.add_argument('--with_stereo', help='Keep stereochemistry (if any) when sanitizing chemicals. Default to '
                                              'False.',
                        default=False, type=lambda x: (str(x).lower() == 'true'))

    subparsers = parser.add_subparsers(help='Input mode')

    parser_inline = subparsers.add_parser('inline', help='Get inputs from command line')
    parser_inline.set_defaults(func=__inline_mode)
    parser_inline.add_argument('--rsmarts', help='Reaction rule SMARTS', required=True)
    # parser_inline.add_argument('--csmiles', help='Chemical SMILES depiction', required=True)
    parser_inline.add_argument('--inchi', help='Chemical inchi depiction', required=True)
    parser_inline.add_argument('--rid', help='Reaction rule ID, optional')
    parser_inline.add_argument('--cid', help='Chemical ID, optional')

    parser_file = subparsers.add_parser('infile', help='Get inputs from files')
    parser_file.set_defaults(func=__infile_mode)
    parser_file.add_argument(
            '--rfile', required=True, help=' '.join([
                    'Reaction rule file.',
                    'Tab separated columns.',
                    'One reaction rule per line.',
                    'Mandatory column: rule_smarts, rule_id.',
                    'Other columns will be ignored.'
                    ])
            )
    parser_file.add_argument(
            '--cfile', required=True, help=' '.join([
                    'Chemical file.',
                    'Tab separated columns.',
                    'One chemical per line.',
                    'Mandatory column: inchi, chem_id.',
                    'Other columns will be ignored.'
                    ])
            )

    return parser


def __cli():
    logging.basicConfig(
            stream=sys.stderr, level=logging.INFO,
            datefmt='%d/%m/%Y %H:%M:%S',
            format='%(asctime)s -- %(levelname)s -- %(message)s'
            )
    parser = __build_arg_parser()
    args = parser.parse_args()
    try:
        results = args.func(args)
        results_json = jsonify(results)
        write(results_json, ofile=args.ofile, compress=args.compress)
    except AttributeError:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    __cli()
