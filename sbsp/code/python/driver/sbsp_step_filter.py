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
from sbsp_alg.filtering import filter_orthologs
from sbsp_general import Environment

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from sbsp_options.msa import MSAOptions

parser = argparse.ArgumentParser("Description of driver.")

parser.add_argument('--pf-data', required=True, help="CSV Data file")
parser.add_argument('--pf-output', required=True, help="Output file containing features")

parser.add_argument('--pf-msa-options', required=True, help="Configuration file containing MSA options")

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
    msa_options = MSAOptions.init_from_dict(env, vars(args))

    filter_orthologs(env, args.pf_data, args.pf_output,
                     msa_options=msa_options)


if __name__ == "__main__":
    main(my_env, parsed_args)
