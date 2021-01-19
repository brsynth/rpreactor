"""
Toolbox to get the most out of rpreactor rules.
"""

import logging
import csv

from rpreactor.rule.burner import RuleBurner


logger = logging.getLogger(__name__)


def _create_db_from_retrorules_v1_0_5(path_retrosmarts_tsv, db_path, with_hs, with_stereo):
    def helper_metabolite(store, cid, smiles):
        if cid not in store:
            store[cid] = smiles
        elif store[cid] != smiles:
            logger.warning(f"Metabolite {cid} is suspiciously associated to distinct SMILES. "
                           f"Only the first one will be considered: {store[cid]} and {smiles}")
    rules = {}
    metabolites = {}  # both substrate and products
    results = set()   # all "obvious" results that directly come from the reaction database (no promiscuity)
    pgroup = {}       # index of a solution rid vs. cid
    # Load all valuable data in-memory
    # NB: is this file small enough that we do not need to chunk it?
    with open(path_retrosmarts_tsv) as hdl:
        for row in csv.DictReader(hdl, delimiter='\t'):
            # each row is a reaction rule automatically generated from a known metabolic reaction
            # each row contains 1 rule...
            # NB: rule identifier may be duplicated over several rows but must match the same RSMARTS
            rid = row["# Rule_ID"]
            rsmarts = row["Rule_SMARTS"]
            diameter = int(row["Diameter"])
            direction = row["Rule_usage"]  # -1, 0, 1 ==> retro, both, forward
            if direction == "both":
                direction = 0
            elif direction == "retro":
                direction = -1
            elif direction == "forward":
                direction = 1
            else:
                raise ValueError(f"Found an unexpected direction for rule {rid}: {direction}")
            if rid not in rules:
                # warning: keys must match database schema
                rules[rid] = {'rd_rule': rsmarts, 'diameter': diameter, 'direction': direction}
            else:
                assert rules[rid]['rd_rule'] == rsmarts, f"UNEXPECTED: rule {rid} from {path_retrosmarts_tsv} has " \
                                                         f"mismatching RSMARTS: {rules[rid]} and {rsmarts}"
            # ... and 1 substrate ...
            sid = row["Substrate_ID"]
            helper_metabolite(metabolites, sid, row["Substrate_SMILES"])
            # ... and N coproducts
            smiles_list = row["Product_SMILES"].split('.')
            tmp_results = []
            pid_stoichio = {}
            for idx, pid in enumerate(row["Product_IDs"].split('.')):
                if pid in pid_stoichio:
                    pid_stoichio[pid] += 1
                else:
                    pid_stoichio[pid] = 1
                    helper_metabolite(metabolites, pid, smiles_list[idx])
                    tmp_results.append((rid, sid, pid))  # TODO: bug, there can be several solutions for the same rule x mol couple
            # each row is also a distinct solution of "1 rule applied on 1 metabolite"... but there can be many
            # especially at low diameters!
            if (rid, sid) not in pgroup:
                pgroup[(rid, sid)] = -1
            else:
                pgroup[(rid, sid)] -= 1
            # record the results of this row (1 by distinct product)
            for rid, sid, pid in tmp_results:
                results.add((sid, rid, pid, pid_stoichio[pid], pgroup[(rid, sid)]))
    # Create the database
    o = RuleBurner(db_path=db_path, with_hs=with_hs, with_stereo=with_stereo)
    o.insert_rsmarts(rules)
    o.insert_smiles(metabolites)
    o.db.executemany("INSERT INTO results VALUES (?,?,?,?,?);", list(results))
    o.create_indexes()


def create_db_from_retrorules(path_retrosmarts_tsv, db_path, with_hs=False, with_stereo=False, version="v1.0"):
    """Convert a RetroRules dataset to a rpreactor-ready sqlite3 database.

    All rules and all molecules will be extracted from the TSV file and imported into the database.
    For more information on RetroRules, see https://retrorules.org/.
    """
    if version.startswith("v1.0"):
        _create_db_from_retrorules_v1_0_5(path_retrosmarts_tsv, db_path, with_hs, with_stereo)
