# Karl Gemayel
# Georgia Institute of Technology
#
# Created: 1/15/20

import logging
import argparse
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
from sbsp_general.labels_comparison_detailed import LabelsComparisonDetailed
from sbsp_io.labels import read_labels_from_file

parser = argparse.ArgumentParser("Description of driver.")

parser.add_argument('--pf-gms2', required=True, help="GMS2 prediction")
parser.add_argument('--pf-sbsp', required=True, help="SBSP prediction")
parser.add_argument('--pf-ncbi', required=True, help="NCBI annotation")

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

    labels_gms2 = read_labels_from_file(pf_gms2)
    labels_sbsp = read_labels_from_file(pf_sbsp)
    labels_ncbi = read_labels_from_file(pf_ncbi)

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


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None
    compare_gms2_sbsp_ncbi(env, args.pf_gms2, args.pf_sbsp, args.pf_ncbi)


if __name__ == "__main__":
    main(my_env, parsed_args)
