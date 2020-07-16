# Karl Gemayel
# Georgia Institute of Technology
#
# Created: 3/17/20
import logging
import argparse
from typing import *

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import sbsp_log  # runs init in sbsp_log and configures logger

# Custom imports
from sbsp_general import Environment
from sbsp_general.blast import gen_cmd_create_blast_database

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from sbsp_general.general import run_shell_cmd

parser = argparse.ArgumentParser("Build Diamond Blast database from sequences file.")

parser.add_argument('--pf-sequences', required=True, help="Fasta file with specific header format. See "
                                                          "documentation for more information.")

parser.add_argument('--pf-db', required=True, help="Path to output file.")

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

    run_shell_cmd(gen_cmd_create_blast_database(args.pf_sequences, args.pf_db, "nucl", True))


if __name__ == "__main__":
    main(my_env, parsed_args)
