# Karl Gemayel
# Georgia Institute of Technology
#
# Created: 1/15/20
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
from sbsp_general.labels_comparison_detailed import LabelsComparisonDetailed
from sbsp_io.general import read_rows_to_list
from sbsp_io.labels import read_labels_from_file
from sbsp_viz.general import FigureOptions
from sbsp_viz.labels_venn import venn_diagram_5prime

parser = argparse.ArgumentParser("Description of driver.")

parser.add_argument('--pf-genome-list', required=True, help="List containing genome information")
parser.add_argument('--pf-gcfid-to-pd-sbsp', required=True, help="CSV file containing GCFID to SBSP run directory")
parser.add_argument('--pf-output-summary', required=True, help="Output summary file")

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


def compare_gms2_sbsp_ncbi(env, pf_gms2, pf_sbsp, pf_ncbi, **kwargs):
    # type: (Environment, str, str, str, Dict[str, Any]) -> Dict[str, Any]

    venn_title = get_value(kwargs, "venn_title", None)
    pf_venn = get_value(kwargs, "pf_venn", os.path.join(env["pd-work"], "venn.pdf"))

    labels_gms2 = read_labels_from_file(pf_gms2, name="GMS2")
    labels_sbsp = read_labels_from_file(pf_sbsp, name="SBSP")
    labels_ncbi = read_labels_from_file(pf_ncbi, name="NCBI")

    lcd = LabelsComparisonDetailed(labels_gms2, labels_sbsp,
                                   name_a="gms2",
                                   name_b="sbsp")

    labels_gms2_sbsp_3p_5p = lcd.intersection("a")

    lcd_2 = LabelsComparisonDetailed(labels_gms2_sbsp_3p_5p, labels_ncbi,
                                     name_a="gms2_sbsp",
                                     name_b="ncbi")

    labels_gms2_sbsp_ncbi_3p_5p = lcd_2.intersection("a")

    venn_diagram_5prime(labels_gms2, labels_sbsp, labels_ncbi, FigureOptions(
        title=venn_title,
        save_fig=pf_venn
    ))

    return {
        "GMS2": len(labels_gms2),
        "SBSP": len(labels_sbsp),
        "NCBI": len(labels_ncbi),
        "GMS2=SBSP": len(labels_gms2_sbsp_3p_5p),
        "GMS2=SBSP=NCBI": len(labels_gms2_sbsp_ncbi_3p_5p)
    }


def compare_gms2_sbsp_ncbi_for_genome_list(env, gil, gcfid_to_pd_sbsp, pf_output_summary, **kwargs):
    # type: (Environment, GenomeInfoList, Dict[str, str], str, Dict[str, Any]) -> None
    list_summary = list()
    for gi in gil:
        pd_genome = os.path.join(env["pd-data"], gi.name)
        pf_gms2 = os.path.join(pd_genome, "runs", "gms2", "gms2.gff")
        pf_sbsp = os.path.join(gcfid_to_pd_sbsp[gi.name], "accuracy", "{}.gff".format(gi.name))
        pf_ncbi = os.path.join(pd_genome, "ncbi.gff")

        name = gi.attributes["name"] if "name" in gi.attributes else gi.name
        ancestor = gi.attributes["ancestor"] if "ancestor" in gi.attributes else ""

        out = compare_gms2_sbsp_ncbi(env, pf_gms2, pf_sbsp, pf_ncbi,
                               venn_title="{}, {}".format(name, ancestor),
                               pf_venn="venn_{}.pdf".format(gi.name))

        out["GCFID"] = gi.name
        out["Name"] = name
        out["Ancestor"] = ancestor

        list_summary.append(out)

    if len(list_summary) > 0:

        ordered_header = ["GCFID", "Name", "Ancestor", "GMS2", "SBSP", "NCBI", "GMS2=SBSP", "GMS2=SBSP=NCBI"]
        remaining_header = sorted(
            set(list_summary[0].keys()).difference(ordered_header)
        )
        header = ordered_header + remaining_header

        df = pd.DataFrame(list_summary, columns=header)
        df.to_csv(pf_output_summary, index=False)


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None

    gil = GenomeInfoList.init_from_file(args.pf_genome_list)
    df = pd.read_csv(args.pf_gcfid_to_pd_sbsp)

    gcfid_to_pd_sbsp = {
        x["gcfid"]: x["pd-sbsp"] for _, x in df.iterrows()
    }

    compare_gms2_sbsp_ncbi_for_genome_list(
        env,
        gil, gcfid_to_pd_sbsp, args.pf_output_summary
    )


if __name__ == "__main__":
    main(my_env, parsed_args)
