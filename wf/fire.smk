r"""
Prepare set of files.

Simple example (from project root folder):
snakemake --directory . --snakefile bin/rule_fire/wf/fire.snake -p --jobs 6 --configfile bin/rule_fire/wf/config.yml --keep-going > wf_fire_mnx20190307.log 2>&1

Cluster example:
screen
unset PYTHONPATH
source activate pyrule
snakemake --directory . --snakefile bin/rule_fire/wf/fire.snake --drmaa " -j yes -q short.q -l h_rt=1:00:00" --drmaa-log-dir logs_cluster --jobs 120 --jobname snakejob.{rule}_{wildcards.inchikey}_r{wildcards.radius}_{jobid}.sh --jobscript wf/sge.sh --latency-wait 120 --keep-going --configfile bin/rule_fire/wf/config.yml > wf_fire_mnx20190112.log 2>&1

Thomas Duigou, INRA, 2018-2019
"""

JOB_DIR = config['JOB_DIR']

# Get InchiKeys and radius from folder content
INCHIKEYS, = glob_wildcards(os.path.join(JOB_DIR, "chemicals/{inchikey}.tsv"))
RADIUS, = glob_wildcards(os.path.join(JOB_DIR, "rules/r{radius}.tsv"))

# Cluster: prevent specific rules to be executed from cluster nodes
localrules: all, init

rule all:
    input:
        expand('{job_dir}/res/{inchikey}_r{radius}.json.gz', inchikey=INCHIKEYS, radius=RADIUS, job_dir=JOB_DIR)

rule init:
    output:
        logdir = directory('{job_dir}/logs')
    shell:
        """
        mkdir -p {output.logdir}
        """

rule fire:
    input:
        init_done = '{job_dir}/logs',
        chem_file = '{job_dir}/chemicals/{inchikey}.tsv',
        rule_file = '{job_dir}/rules/r{radius}.tsv'
    output:
        res_file = '{job_dir}/res/{inchikey}_r{radius}.json.gz'
    log:
        '{job_dir}/logs/fire_{inchikey}_r{radius}.log'
    shell:
        """
        python reactor/cli.py \
            --ofile {output.res_file} \
            --compress \
            infile --rfile {input.rule_file} --cfile {input.chem_file} \
            > {log} 2>&1
        """
