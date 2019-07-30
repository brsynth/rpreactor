r"""
Prepare set of files.

Example (to be excuted from project root folder):
snakemake -s wf/prepare.smk --configfile wf/config.yml -p -F

Thomas Duigou, INRA, 2018-2019
"""

 
JOB_DIR = config['JOB_DIR']
IN_RULES = config['IN_RULES']
RADIUS_TO_KEEP = config['RADIUS_TO_KEEP']
MIN_OCC_CUTOFF = config['MIN_OCC_CUTOFF']
MAX_OCC_CUTOFF = config['MAX_OCC_CUTOFF']
MIN_MW_CUTOFF = config['MIN_MW_CUTOFF']
MAX_MW_CUTOFF = config['MAX_MW_CUTOFF']
FILTER_RULES = config['FILTER_RULES'] if 'FILTER_RULES' in config else []


rule all:
    input:
        '{job_dir}/rules'.format(job_dir=JOB_DIR),
        '{job_dir}/stats'.format(job_dir=JOB_DIR),
        '{job_dir}/chemicals'.format(job_dir=JOB_DIR)
        

rule select_rules:
    input:
        rule_file = IN_RULES
    output:
        rule_dir = directory('{job_dir}/rules'),
        stat_dir = directory('{job_dir}/stats')
    params:
        radius_to_keep = RADIUS_TO_KEEP,
        rule_filter_files = FILTER_RULES
    run:
        import os
        import csv
        all_rules = dict()
        rules_to_filter = set()
        rules_cnt = dict()
        os.makedirs(output[0], exist_ok=True)
        os.makedirs(output[1], exist_ok=True)
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
        # Write
        for radius in sorted(list(all_rules.keys())):
            rule_written = set()
            # Wish list
            opath = os.path.join(output['rule_dir'], 'r{:02d}.tsv'.format(radius))
            with open(opath, 'w') as ofh:
                fields = ['rule_id', 'rule_smarts']
                writer = csv.DictWriter(ofh, fieldnames=fields, quoting=csv.QUOTE_NONE, delimiter='\t')
                writer.writeheader()
                for rule_id, rule_smarts in all_rules[radius].items():
                    if (rule_id not in rules_to_filter) and (rule_id not in rule_written):
                        writer.writerow({'rule_id': rule_id, 'rule_smarts': rule_smarts})
                        rules_to_filter.add(rule_id)
                        rule_written.add(rule_id)
            # Stats
            opath = os.path.join(output['stat_dir'], 'r{:02d}.tsv'.format(radius))
            with open(opath, 'w') as ofh:
                fields = ['rule_id', 'count', 'written']
                writer = csv.DictWriter(ofh, fieldnames=fields, quoting=csv.QUOTE_NONE, delimiter='\t')
                writer.writeheader()
                for rule_id, count in sorted(rules_cnt[radius].items(), key=lambda x: x[1], reverse=True):
                    writer.writerow({'rule_id': rule_id, 'count': count, 'written': rule_id in rule_written})

rule select_chemicals:
    input:
        rule_file = IN_RULES,
        inchi_file = config['IN_CHEMICALS']
    output:
        chem_out = '{job_dir}/chem_wish_list.tsv',
        chem_stat = '{job_dir}/chem_stats.tsv'
    params:
        min_occ_cutoff = MIN_OCC_CUTOFF,
        max_occ_cutoff = MAX_OCC_CUTOFF,
        min_mw_cutoff = MIN_MW_CUTOFF,
        max_mw_cutoff = MAX_MW_CUTOFF
    run:
        import csv
        import math

        from rdkit import RDLogger
        from rdkit.Chem import MolFromInchi
        from rdkit.Chem.Descriptors import MolWt

        rd_logger = RDLogger.logger()
        rd_logger.setLevel(RDLogger.ERROR)

        def get_molecular_weight(inchi):
            try:
                return round(MolWt(MolFromInchi(inchi)), 2)
            except Exception:
                return 0

        # Get inchi info
        chem_id_2_inchikey = {}
        chem_from_inchikey = {}
        with open(input['inchi_file']) as ifh:
            reader = csv.DictReader(ifh, delimiter='\t')
            for row in reader:
                if row['inchikey'] not in chem_from_inchikey:
                    chem_from_inchikey[row['inchikey']] = {
                        'inchikey': row['inchikey'],
                        'inchi': row['inchi'],
                        'xref': {row['cid']},
                        'count': 0,
                        'mw': get_molecular_weight(row['inchi'])
                    }
                else:
                    assert chem_from_inchikey[row['inchikey']]['inchi'] == row['inchi']
                    chem_from_inchikey[row['inchikey']]['xref'].add(row['cid'])
                if row['cid'] not in chem_id_2_inchikey:
                    chem_id_2_inchikey[row['cid']] = row['inchikey']
                else:
                    assert chem_id_2_inchikey[row['cid']] == row['inchikey']

        # Get chemical used in rules
        all_radius = set()
        with open(input['rule_file']) as ifh:
            reader = csv.DictReader(ifh, delimiter='\t')
            for row in reader:
                # Chemical
                try:
                    assert row['Substrate_ID'] in chem_id_2_inchikey
                except KeyError as e:
                    print(e)
                    continue
                inchikey = chem_id_2_inchikey[row['Substrate_ID']]
                chem_from_inchikey[inchikey]['count'] += 1
                # Check radius
                radius = str(int(row['Diameter'])//2)
                if radius not in all_radius:
                    all_radius.add(radius)

        # Compute per diameter average usage of chemicals
        for item in chem_from_inchikey.values():
            item['average_count'] = int(math.ceil(item['count'] / len(all_radius)))

        # Select chemicals to keep
        to_keep = set()
        for inchikey, item in chem_from_inchikey.items():
            if (
                    params['max_occ_cutoff'] >= item['average_count'] >= params['min_occ_cutoff'] and
                    params['max_mw_cutoff'] >= item['mw'] >= params['min_mw_cutoff']
            ):
                to_keep.add(inchikey)

        # Write the wish list
        with open(output['chem_out'], 'w') as ofh:
            fields = ['chem_id', 'inchi', 'xref', 'molecular_weight', 'average_count']
            writer = csv.DictWriter(ofh, fieldnames=fields, quoting=csv.QUOTE_NONE, delimiter='\t')
            writer.writeheader()
            for item in chem_from_inchikey.values():
                if item['inchikey'] in to_keep:
                    writer.writerow({
                            'chem_id': item['inchikey'],
                            'inchi': item['inchi'],
                            'xref': ','.join(item['xref']),
                            'average_count': item['average_count'],
                            'molecular_weight': item['mw']
                            })

        # Write stats
        with open(output['chem_stat'], 'w') as ofh:
            fields = ['chem_id', 'written', 'average_count', 'molecular_weight', 'xref']
            writer = csv.DictWriter(ofh, fieldnames=fields, quoting=csv.QUOTE_NONE, delimiter='\t')
            writer.writeheader()
            for item in sorted(chem_from_inchikey.values(), key=lambda x: x['average_count'], reverse=True):
                writer.writerow({
                        'chem_id': item['inchikey'],
                        'written': item['inchikey'] in to_keep,
                        'average_count': item['average_count'],
                        'molecular_weight': item['mw'],
                        'xref': ','.join(item['xref'])
                        })


rule split:
    input:
        chem_list = '{job_dir}/chem_wish_list.tsv'
    output:
        outdir = directory('{job_dir}/chemicals')
    params:
        debug = True
    run:
        import os
        import csv
        os.makedirs(output['outdir'], exist_ok=True)

        with open(input['chem_list']) as ifh:
            reader = csv.reader(ifh, delimiter='\t')
            fields = next(reader)

        with open(input['chem_list']) as ifh:
            reader = csv.DictReader(ifh, delimiter='\t')
            cnt = 0
            for row in reader:
                cnt += 1
                chem_id = row['chem_id']
                chem_path = os.path.join(output['outdir'], chem_id + '.tsv')
                with open(chem_path, 'w') as ofh:
                    writer = csv.DictWriter(ofh, fieldnames=fields, quoting=csv.QUOTE_NONE, delimiter='\t')
                    writer.writeheader()
                    writer.writerow(row)
                    if params.debug and cnt > 20:
                        break
