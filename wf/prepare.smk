r"""
Prepare set of files.

Example (to be excuted from project root folder):
snakemake --directory . --snakefile bin/rule_fire/wf/prepare.snake -p --configfile bin/rule_fire/wf/config.yml > wf_prepare_mnx20190307.log 2>&1


Thomas Duigou, INRA, 2018-2019
"""

 
JOB_DIR = config['JOB_DIR']
IN_RULES = config['IN_RULES']
RADIUS_TO_KEEP = config['RADIUS_TO_KEEP']
MAX_CUTOFF = config['MAX_CUTOFF']
MIN_CUTOFF = config['MIN_CUTOFF']
FILTER_RULES = config['FILTER_RULES'] if 'FILTER_RULES' in config else []


rule all:
    input:
        '{job_dir}/rules'.format(job_dir=JOB_DIR),
        '{job_dir}/chemicals'.format(job_dir=JOB_DIR)
        

rule select_rules:
    input:
        rule_file = IN_RULES
    output:
        rule_dir = directory('{job_dir}/rules')
    params:
        radius_to_keep = RADIUS_TO_KEEP,
        rule_filter_files = FILTER_RULES
    run:
        import csv
        all_rules = dict()
        rules_to_filter = set()
        rules_cnt = dict()
        os.makedirs(output[0], exist_ok=True)
        # Get list of rules to filter in any case
        if len(params['rule_filter_files']) > 0:
            for ifilter in params['rule_filter_files']:
                print(ifilter)
                with open(ifilter, 'r') as ifh:
                    reader = csv.DictReader(ifh, delimiter='\t')
                    for row in reader:
                        rules_to_filter.add(row['Rule_ID'])
        # Get info
        with open(input['rule_file'], 'r') as ifh:
            reader = csv.DictReader(ifh, delimiter='\t')
            for row in reader:
                rule_hash = row['# Rule_ID'].split('-')[2]
                rule_smarts = row['Rule_SMARTS']
                # Radius checking
                radius = int(row['Diameter'])//2
                if radius not in params['radius_to_keep']:
                    continue
                if radius not in all_rules.keys():
                    all_rules[radius] = dict()
                    rules_cnt[radius] = dict()
                # Store rules
                if rule_hash in all_rules[radius]:
                    assert all_rules[radius][rule_hash] == rule_smarts
                    rules_cnt[radius][rule_hash] += 1
                else:
                    all_rules[radius][rule_hash] = rule_smarts
                    rules_cnt[radius][rule_hash] = 1
        # List of rules
        for radius in sorted(list(all_rules.keys())):
            rule_written = set()
            # Write wishlist
            opath = os.path.join(output[0], 'r{:02d}.tsv'.format(radius))
            with open(opath, 'w') as ofh:
                fields = ['Rule_ID', 'Rule_SMARTS']
                writer = csv.DictWriter(ofh, fieldnames=fields, quoting=csv.QUOTE_NONE, delimiter='\t')
                writer.writeheader()
                for rule_id, rule_smarts in all_rules[radius].items():
                    if (rule_id not in rules_to_filter) and (rule_id not in rule_written):
                        writer.writerow({'Rule_ID': rule_id, 'Rule_SMARTS': rule_smarts})
                        rules_to_filter.add(rule_id)
                        rule_written.add(rule_id)
            # Write stats
            opath = os.path.join(output[0], 'stats_r{:02d}.tsv'.format(radius))
            with open(opath, 'w') as ofh:
                fields = ['Rule_ID', 'Count', 'Written']
                writer = csv.DictWriter(ofh, fieldnames=fields, quoting=csv.QUOTE_NONE, delimiter='\t')
                writer.writeheader()
                for rule_id, count in sorted(rules_cnt[radius].items(), key=lambda x: x[1], reverse=True):
                    writer.writerow({'Rule_ID': rule_id, 'Count': count, 'Written': rule_id in rule_written})

rule select_chemicals:
    input:
        rule_file = IN_RULES
    output:
        chem_out = '{job_dir}/chems_wishlist.tsv',
        chem_stat = '{job_dir}/chems_stat.tsv'
    params:
        max_cutoff = MAX_CUTOFF,
        min_cutoff = MIN_CUTOFF
    log:
        '{job_dir}/select_chemicals.log'
    run:
        import csv
        import math
        import rdkit
        import logging
        from rdkit import Chem
        from rdkit import RDLogger  # To remove InChI warnings
        logging.basicConfig(filename=log[0])
        rd_logger = RDLogger.logger()
        all_smiles = dict()
        all_chems = dict()
        all_radius = set()
        
        # Get info
        id2inchikey = dict()  # Chemical MNX ID => InchiKey
        with open(input['rule_file'], 'r') as ifh:
            reader = csv.DictReader(ifh, delimiter='\t')
            for row in reader:
                chem_id = row['Substrate_ID']
                chem_smiles = row['Substrate_SMILES']
                # Check radius
                radius = str(int(row['Diameter'])//2)
                if radius not in all_radius:
                    all_radius.add(radius)
                # Get SMILES
                if chem_id in all_smiles.keys():
                    assert all_smiles[chem_id]['smiles'] == chem_smiles
                    all_smiles[chem_id]['count'] += 1
                else:
                    all_smiles[chem_id] = dict()
                    all_smiles[chem_id]['smiles'] = chem_smiles
                    all_smiles[chem_id]['count'] = 1
        # Gather similar chemicals based on InchiCkey
        for chem_id in all_smiles.keys():
            smiles_ori = all_smiles[chem_id]['smiles']
            count = all_smiles[chem_id]['count']
            # Standardize and get depictions
            rdmol = Chem.MolFromSmiles(smiles_ori)
            Chem.SanitizeMol(rdmol)
            rd_logger.setLevel(RDLogger.ERROR)  # Decrease Inchi verbosity
            inchi = Chem.MolToInchi(rdmol)
            inchikey = Chem.MolToInchiKey(rdmol)
            rd_logger.setLevel(RDLogger.WARNING)
            # Double check conversion
            try:
                rdmol_check = Chem.MolFromSmiles(Chem.MolToSmiles(rdmol))
                Chem.SanitizeMol(rdmol_check)
                rdmol_check = Chem.AddHs(rdmol_check)
                smiles_check = Chem.MolToSmiles(rdmol_check, allHsExplicit=True)
                assert smiles_ori == smiles_check
            except AssertionError:
                logging.warning('Back and forth produces different SMILES')
                logging.warning('First SMILES: {}'.format(smiles_ori))
                logging.warning('Secon SMILES: {}'.format(smiles_check))
            # Store
            if inchikey in all_chems.keys():
                try:
                    assert all_chems[inchikey]['inchi'] == inchi
                    assert all_chems[inchikey]['smiles'] == smiles_ori
                except AssertionError:
                    logging.warning('Same InchiKey but different SMILES...')
                    logging.warning('INCHIKEY: {}'.format(inchikey))
                    logging.warning('SMILES 1: {}'.format(all_chems[inchikey]['smiles']))
                    logging.warning('SMILES 2: {}'.format(smiles_ori))
            else:
                all_chems[inchikey] = dict()
                all_chems[inchikey]['inchi'] = inchi
                all_chems[inchikey]['smiles'] = smiles_ori
                all_chems[inchikey]['xref'] = set()
                all_chems[inchikey]['count'] = 0
                all_chems[inchikey]['average_use'] = None
            all_chems[inchikey]['xref'].add(chem_id)
            all_chems[inchikey]['count'] += count
        # Compute average usage
        for inchikey, info in all_chems.items():
            average_use = info['count'] / len(all_radius)
            all_chems[inchikey]['average_use'] = int(math.ceil(average_use))
        # Select chemicals to keep (based on average use)
        chems_to_keep = set()
        for inchikey, info in all_chems.items():
            if (info['average_use'] <= params['max_cutoff']) and (info['average_use'] >= params['min_cutoff']):
                chems_to_keep.add(inchikey)
        # Write wishlist of chemicals
        with open(output['chem_out'], 'w') as ofh:
            fields = ['Chemical_ID', 'Chemical_SMILES', 'Chemical_INCHI', 'Chemical_INCHIKEY', 'Chemical_XREF']
            writer = csv.DictWriter(ofh, fieldnames=fields, quoting=csv.QUOTE_NONE, delimiter='\t')
            writer.writeheader()
            for inchikey, info in all_chems.items():
                if inchikey in chems_to_keep:
                    writer.writerow({
                            'Chemical_ID': inchikey,
                            'Chemical_SMILES': info['smiles'],
                            'Chemical_INCHI': info['inchi'],
                            'Chemical_INCHIKEY': inchikey,
                            'Chemical_XREF': ','.join(info['xref'])
                            })
        # Write chemical count stats
        with open(output['chem_stat'], 'w') as ofh:
            fields = ['Chemical_ID', 'Count', 'Chemical_XREF', 'Written']
            writer = csv.DictWriter(ofh, fieldnames=fields, quoting=csv.QUOTE_NONE, delimiter='\t')
            writer.writeheader()
            for inchikey, info in sorted(all_chems.items(), key=lambda x: x[1]['average_use'], reverse=True):
                writer.writerow({
                        'Chemical_ID': inchikey,
                        'Count': info['average_use'],
                        'Chemical_XREF': ','.join(info['xref']),
                        'Written': chem_id in chems_to_keep
                        })

rule split:
    input:
        chem_list = '{job_dir}/chems_wishlist.tsv'
    output:
        outdir = directory('{job_dir}/chemicals')
    run:
        import os
        import csv
        os.makedirs(output['outdir'], exist_ok=True)

        with open(input['chem_list'], 'r') as ifh:
            reader = csv.reader(ifh, delimiter='\t')
            fields = next(reader)

        with open(input['chem_list'], 'r') as ifh:
            reader = csv.DictReader(ifh, delimiter='\t')
            for row in reader:
                chem_id = row['Chemical_ID']
                chem_path = os.path.join(output['outdir'], chem_id + '.tsv')
                with open(chem_path, 'w') as ofh:
                    writer = csv.DictWriter(ofh, fieldnames=fields, quoting=csv.QUOTE_NONE, delimiter='\t')
                    writer.writeheader()
                    writer.writerow(row)