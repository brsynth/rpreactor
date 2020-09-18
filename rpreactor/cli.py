r"""
Apply one/many rules on one/many substrates.

Thomas Duigou, INRA, 2018
"""

import sys
import logging
import argparse

from rpreactor.rule.burner import RuleBurner


def __cli():
    """Command line interface."""

    help = "Apply rules on chemicals."

    def inline_mode(args):
        """Execution mode to be used when a single rule and a single chemical
        are provided through CLI.
        """
        r = RuleBurner(
                rsmarts_list=[args.rsmarts], inchi_list=[args.inchi],
                rid_list=[args.rid], cid_list=[args.cid],
                match_timeout=args.match_timeout, fire_timeout=args.fire_timeout,
                ofile=args.ofile, compress=args.compress, 
                with_hs=args.with_hs, with_stereo=args.with_stereo
                )
        r.compute()
        r.write_json()

    def infile_mode(args):
        """
        Execution mode to be used when rules and chemicals are provided
        in CSV files.
        """

        rsmarts_list = list()
        rids_list = list()

        inchi_list = list()
        cids_list = list()

        import csv

        with open(args.rfile, 'r') as ifh:
            reader = csv.DictReader(ifh, delimiter='\t')
            for row in reader:
                rsmarts_list.append(row['rule_smarts'].strip())
                rids_list.append(row['rule_id'].strip())

        with open(args.cfile, 'r') as ifh:
            reader = csv.DictReader(ifh, delimiter='\t')
            for row in reader:
                inchi_list.append(row['inchi'].strip())
                cids_list.append(row['chem_id'].strip())

        r = RuleBurner(
                rsmarts_list=rsmarts_list, inchi_list=inchi_list,
                rid_list=rids_list, cid_list=cids_list,
                match_timeout=args.match_timeout, fire_timeout=args.fire_timeout,
                ofile=args.ofile, compress=args.compress,
                with_hs=args.with_hs, with_stereo=args.with_stereo
                )
        r.compute()
        r.write_json()

    parser = argparse.ArgumentParser(description=help)
    parser.add_argument('--match_timeout', help='Rule matching timeout. Default: 1.', default=1, type=int)
    parser.add_argument('--fire_timeout', help='Rule firing timeout. Default: 1.', default=1, type=int)
    parser.add_argument('--ofile', help='Output file to store results. Default to STDOUT if none provided')
    parser.add_argument('--compress', action='store_true', help='Enable gzip compression (only when output to file).')
    parser.add_argument('--with_hs', help='Enable explicit Hs when sanitizing chemicals. Default to False.',
                        default=False, type=lambda x: (str(x).lower() == 'true'))
    parser.add_argument('--with_stereo', help='Keep stereochemistry (if any) when sanitizing chemicals. Default to '
                                              'False.',
                        default=False, type=lambda x: (str(x).lower() == 'true'))

    subparsers = parser.add_subparsers(help='Input mode')

    parser_inline = subparsers.add_parser('inline', help='Get inputs from command line')
    parser_inline.set_defaults(func=inline_mode)
    parser_inline.add_argument('--rsmarts', help='Reaction rule SMARTS', required=True)
    # parser_inline.add_argument('--csmiles', help='Chemical SMILES depiction', required=True)
    parser_inline.add_argument('--inchi', help='Chemical inchi depiction', required=True)
    parser_inline.add_argument('--rid', help='Reaction rule ID, optional')
    parser_inline.add_argument('--cid', help='Chemical ID, optional')

    parser_file = subparsers.add_parser('infile', help='Get inputs from files')
    parser_file.set_defaults(func=infile_mode)
    parser_file.add_argument(
            '--rfile', required=True, help=' '.join([
                    'Reaction rule file.',
                    'Tab separated columns.',
                    'One reaction rule per line.',
                    'Mandatory column: rule_smarts.',
                    'Optional column: rule_id.',
                    'Other columns will be ignored.'
                    ])
            )
    parser_file.add_argument(
            '--cfile', required=True, help=' '.join([
                    'Chemical file.',
                    'Tab separated columns.',
                    'One chemical per line.',
                    'Mandatory column: inchi.',
                    'Optional column: chem_id.',
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
