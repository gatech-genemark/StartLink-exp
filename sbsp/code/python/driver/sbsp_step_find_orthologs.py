# Karl Gemayel
# Georgia Institute of Technology
#
# Created: 12/16/19

import logging
import argparse
from typing import *

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import sbsp_log  # runs init in sbsp_log and configures logger

# Custom imports
from sbsp_alg.ortholog_finder import get_orthologs_from_files
from sbsp_general import Environment
from sbsp_alg.sbsp_steps import sbsp_step_get_orthologs

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #

parser = argparse.ArgumentParser("Description of driver.")

parser.add_argument('--pf-q-list', required=True, help="Path to file containing query genome GCFIDs")
parser.add_argument('--pf-t-list', required=True, help="Path to file containing target genome GCFIDs")

parser.add_argument('--fn-q-labels', required=True, help="Name of labels file for query genomes")
parser.add_argument('--fn-t-labels', required=True, help="Name of labels file for target genomes")

parser.add_argument('--clean', required=False, default=False, action="store_true", help="If set, clean up after BLAST")

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

    get_orthologs_from_files(env, args.pf_q_list, args.pf_t_list, args.pf_output,
                             fn_q_labels=args.fn_q_labels,
                             fn_t_labels=args.fn_t_labels,
                             clean=args.clean)


if __name__ == "__main__":
    main(my_env, parsed_args)
