# Karl Gemayel
# Georgia Institute of Technology
#
# Created:

import logging
import argparse
import pandas as pd
from typing import *

# noinspection PyUnresolvedReferences
import pathmagic                        # add path to custom library

# Custom library imports
import sbsp_general
from sbsp_general import Environment
from sbsp_argparse.msa import add_msa_options
from sbsp_options.msa import MSAOptions
import sbsp_alg.feature_generation

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #

parser = argparse.ArgumentParser("Description of driver.")

parser.add_argument('--pf-data', required=True, help="File containing path to MSA file(s), under column pf-msa-output")
parser.add_argument('--pf-features', required=True, help="Output file with features")

add_msa_options(parser)

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
logger = logging.getLogger("logger")                    # type: logging.Logger


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None

    msa_options = MSAOptions.init_from_dict(env, vars(args))

    df_data = pd.read_csv(args.pf_data, header=0)
    df_features = sbsp_alg.feature_generation.generate_features_for_msa_from_df(df_data, msa_options,
                                                                                max_number_downstream=3)

    df_features.to_csv(args.pf_features, index=False)


if __name__ == "__main__":
    main(my_env, parsed_args)
