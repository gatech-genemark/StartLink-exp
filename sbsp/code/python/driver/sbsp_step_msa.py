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
from sbsp_alg.msa import run_sbsp_msa
from sbsp_general import Environment

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from sbsp_options.sbsp import SBSPOptions

parser = argparse.ArgumentParser("Runs the SBSP MSA step.")

parser.add_argument('--pf-data', required=True, help="CSV Data file")
parser.add_argument('--pf-output', required=True, help="Output file containing features")

parser.add_argument('--pf-msa-options', required=True, help="Configuration file containing MSA options")

parser.add_argument('--dn-msa-output', required=False, default=None, type=Union[str],
                    help="If set, multiple sequence alignments will be written to files in directory with "
                         "this name (in working directory)")

parser.add_argument('--msa-output-start', required=False, default=0, type=Union[int],
                    help="Counter start of output filenames")

parser.add_argument('--upstream-length-nt', required=False, default=None, type=Union[int],
                    help="The number of upstream nucleotides to use in MSA")
parser.add_argument('--downstream-length-nt', required=False, default=None, type=Union[int],
                    help="The number of downstream nucleotides to use in MSA")

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
    msa_options = SBSPOptions.init_from_dict(env, vars(args))

    run_sbsp_msa(env, args.pf_data, args.pf_output,
                 msa_options=msa_options,
                 dn_msa_output=args.dn_msa_output,
                 msa_output_start=args.msa_output_start,
                 upstream_length_nt=args.upstream_length_nt,
                 downstream_length_nt=args.downstream_length_nt)


if __name__ == "__main__":
    main(my_env, parsed_args)
