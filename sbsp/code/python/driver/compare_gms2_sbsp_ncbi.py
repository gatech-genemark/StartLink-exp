# Karl Gemayel
# Georgia Institute of Technology
#
# Created: 1/15/20
import logging
import argparse
import os
from typing import *

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import sbsp_log  # runs init in sbsp_log and configures logger

# Custom imports
from sbsp_general import Environment

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from sbsp_general.general import get_value
from sbsp_general.labels_comparison_detailed import LabelsComparisonDetailed
from sbsp_io.labels import read_labels_from_file
from sbsp_viz.general import FigureOptions
from sbsp_viz.labels_venn import venn_diagram_5prime

parser = argparse.ArgumentParser("Description of driver.")

parser.add_argument('--pf-gms2', required=True, help="GMS2 prediction")
parser.add_argument('--pf-sbsp', required=True, help="SBSP prediction")
parser.add_argument('--pf-ncbi', required=True, help="NCBI annotation")

parser.add_argument('--venn-title', required=False, default=None, help="Title for venn diagram")
parser.add_argument('--pf-venn', required=False, default=None, help="Name of venn diagram file")

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
    # type: (Environment, str, str, str, Dict[str, Any]) -> None

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

    out = "gms2,sbsp,ncbi,gms2_sbsp,gms2_sbsp_ncbi"
    out += "\n{},{},{},{},{}".format(
        len(labels_gms2), len(labels_sbsp), len(labels_ncbi), len(labels_gms2_sbsp_3p_5p),
        len(labels_gms2_sbsp_ncbi_3p_5p)
    )

    print(out)

    venn_diagram_5prime(labels_gms2, labels_sbsp, labels_ncbi, FigureOptions(
        title=venn_title,
        save_fig=pf_venn
    ))


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None
    compare_gms2_sbsp_ncbi(env, args.pf_gms2, args.pf_sbsp, args.pf_ncbi, venn_title=args.venn_title,
                           pf_venn=args.pf_venn)


if __name__ == "__main__":
    main(my_env, parsed_args)
