# Karl Gemayel
# Georgia Institute of Technology
#
# Created: 3/17/20

import logging
import argparse
import os
import pandas as pd
from typing import *

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import sbsp_log  # runs init in sbsp_log and configures logger

# Custom imports
from sbsp_container.genome_list import GenomeInfoList, GenomeInfo
from sbsp_general import Environment

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from sbsp_general.general import get_value
from sbsp_general.labels import Labels
from sbsp_general.labels_comparison_detailed import LabelsComparisonDetailed
from sbsp_io.labels import read_labels_from_file

parser = argparse.ArgumentParser("Description of driver.")

parser.add_argument('--pf-genome-list', required=True, help="List containing genome information")
parser.add_argument('--pf-gcfid-to-pd-sbsp', required=True, help="CSV file containing GCFID to SBSP run directory")

parser.add_argument('--pd-work', required=False, default=None, help="Path to working directory")
parser.add_argument('--pd-data', required=False, default=None, help="Path to data directory")
parser.add_argument('--pd-results', required=False, default=None, help="Path to results directory")
parser.add_argument("-l", "--log", dest="loglevel", choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                    help="Set the logging level", default='WARNING')

parsed_args = parser.parse_args()

# ------------------------------ #
#           Main Code            #
# ------------------------------ #

# Load environment variables
my_env = Environment(pd_data=parsed_args.pd_data,
                     pd_work=parsed_args.pd_work,
                     pd_results=parsed_args.pd_results)

# Setup logger
logging.basicConfig(level=parsed_args.loglevel)
logger = logging.getLogger("logger")  # type: logging.Logger


def os_join(*args):
    # type: (List[Any]) -> str
    return os.path.join(*args)


def get_stats_a_from_b_3p(labels_a, labels_b, output):
    # type: (Labels, Labels, Dict[str, Any]) -> None

    lcd = LabelsComparisonDetailed(labels_a, labels_b)

    a_name = labels_a.name if labels_a.name is not None else "a"
    b_name = labels_b.name if labels_b.name is not None else "b"

    output.update({
        a_name: len(labels_a),
        b_name: len(labels_b),
        "Number 3p match: {} from {}".format(a_name, b_name): len(lcd.match_3p("a")),
        "Percentage 3p match: {} from {}".format(a_name, b_name): len(lcd.match_3p("a"))
    })


def get_stats_a_from_b_3p_by_upstream(labels_a, labels_b, output):
    # type: (Labels, Labels, Dict[str, Any]) -> None

    a_name = labels_a.name if labels_a.name is not None else "a"
    b_name = labels_b.name if labels_b.name is not None else "b"


def get_stats_sn_sp(labels_a, labels_b, output):
    # type: (Labels, Labels, Dict[str, Any]) -> None

    lcd = LabelsComparisonDetailed(labels_a, labels_b)

    a_name = labels_a.name if labels_a.name is not None else "a"
    b_name = labels_b.name if labels_b.name is not None else "b"

    a_total = len(labels_a)
    b_total = len(labels_b)

    match_3p = len(lcd.match_3p("a"))
    match_3p_5p = len(lcd.match_3p_5p("a"))

    sp = 0 if a_total == 0 else match_3p / float(a_total)
    sn = 0 if match_3p == 0 else match_3p_5p / float(match_3p)

    output.update({
        a_name: a_total,
        b_name: b_total,
        "Number 3p match: {} from {}".format(a_name, b_name): match_3p,
        "Percentage 3p match: {} from {}".format(a_name, b_name): sp,
        "Number 5p-3p match: {} from {}".format(a_name, b_name): match_3p_5p,
        "Percentage 5p-3p match: {} from {}".format(a_name, b_name): sn,
    })


def analyze_predictions_on_verified_genes(env, gi, pd_sbsp, **kwargs):
    # type: (Environment, GenomeInfo, str, Dict[str, Any]) -> Dict[str, Any]
    pd_gcfid = os_join(env["pd-data"], gi.name)

    pf_sbsp = os_join(pd_sbsp, "accuracy", "{}.gff".format(gi.name))
    pf_gms2 = os_join(pd_gcfid, "runs", "gms2", "gms2.gff")
    pf_verified = os_join(pd_gcfid, "verified.gff")
    pf_ncbi = os_join(pd_gcfid, "ncbi.gff")

    kwargs_labels = {"ignore_frameshifted": True, "ignore_partial": True}

    labels_sbsp = read_labels_from_file(pf_sbsp, name="SBSP", **kwargs_labels)
    labels_verified = read_labels_from_file(pf_verified, name="Verified", **kwargs_labels)
    labels_gms2 = read_labels_from_file(pf_gms2, name="GMS2", **kwargs_labels)
    labels_ncbi = read_labels_from_file(pf_ncbi, name="NCBI", **kwargs_labels)

    labels_sbsp_eq_gms2 = LabelsComparisonDetailed(labels_sbsp, labels_gms2).match_3p_5p("a")
    labels_sbsp_eq_gms2.name = "GMS2=SBSP"

    stats = dict()

    # Stats: Verified genes from NCBI annotation
    get_stats_a_from_b_3p(labels_verified, labels_ncbi, stats)

    # Stats: Verified genes from GMS2 annotation
    get_stats_a_from_b_3p(labels_verified, labels_gms2, stats)

    # Stats: Verified genes from SBSP annotation
    get_stats_a_from_b_3p(labels_verified, labels_sbsp, stats)

    # Stats: Verified genes from SBSP annotation
    get_stats_a_from_b_3p_by_upstream(labels_verified, labels_ncbi, stats)

    # Stats: SBSP Accuracy on verified
    get_stats_sn_sp(labels_verified, labels_sbsp, stats)

    # Stats: GMS2=SBSP Accuracy on verified
    get_stats_sn_sp(labels_verified, labels_sbsp_eq_gms2, stats)

    return stats


def analyze_predictions_on_verified_genes_for_genome_list(env, gil, gcfid_to_pd_sbsp, **kwargs):
    # type: (Environment, GenomeInfoList, Dict[str, str], Dict[str, Any]) -> None

    fn_prefix = get_value(kwargs, "fn_prefix", "", default_if_none=True)

    info_per_gcfid = dict()

    for gi in gil:
        gcfid = gi.name
        try:
            pd_sbsp = gcfid_to_pd_sbsp[gcfid]
            info_per_gcfid[gcfid] = analyze_predictions_on_verified_genes(env, gi, pd_sbsp, **kwargs)
            info_per_gcfid[gcfid]["Genome"] = gi.attributes["name_txt"]
            info_per_gcfid["gi"] = gi
        except KeyError:
            logger.warning("Couldn't get SBSP directory for: {}".format(gcfid))

    list_stats = [info_per_gcfid[x] for x in info_per_gcfid.keys()]

    df = pd.DataFrame(list_stats)

    df.to_csv(os_join(env["pd-work"], "{}summary.csv".format(fn_prefix)))


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None

    gil = GenomeInfoList.init_from_file(args.pf_genome_list)
    df = pd.read_csv(args.pf_gcfid_to_pd_sbsp)

    gcfid_to_pd_sbsp = {
        x["gcfid"]: x["pd-sbsp"] for _, x in df.iterrows()
    }

    analyze_predictions_on_verified_genes_for_genome_list(env, gil, gcfid_to_pd_sbsp)


if __name__ == "__main__":
    main(my_env, parsed_args)
