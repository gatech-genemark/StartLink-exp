# Karl Gemayel
# Georgia Institute of Technology
#
# Created: 04/29/2020
import logging
import argparse
from typing import *

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import sbsp_log  # runs init in sbsp_log and configures logger

# Custom imports
from sbsp_container.genome_list import GenomeInfoList
from sbsp_general import Environment

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #

parser = argparse.ArgumentParser("Collect result statistics per query.")

parser.add_argument('--pf-genome-list', required=True, help="File containing genome list")

parser.add_argument('--dn-sbsp', default="sbsp", help="Directory name for SBSP run (under genome/runs/)")
parser.add_argument('--dn-gms2', default="gms2", help="Directory name for GMS2 run (under genome/runs/)")
parser.add_argument('--dn-prodigal', default="prodigal", help="Directory name for Prodigal run (under genome/runs/)")

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


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None

    gil = GenomeInfoList.init_from_file(args.pf_genome_list)
    gcfid_to_pd_sbsp


if __name__ == "__main__":
    main(my_env, parsed_args)
