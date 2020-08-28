
# Karl Gemayel
# Georgia Institute of Technology
#
# Created: 3/17/20
import logging
import argparse

import pandas as pd

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import sbsp_log  # runs init in sbsp_log and configures logger

# Custom imports
from sbsp_general import Environment

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from sbsp_general.GMS2Noncoding import GMS2Noncoding
from sbsp_general.MotifModel import MotifModel
from sbsp_general.shelf import relative_entropy
from sbsp_io.objects import load_obj
import sbsp_viz.sns as sns

parser = argparse.ArgumentParser("Description of driver.")


parser.add_argument('--pf-data', required=True, help="File containing motifs")

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
    df = load_obj(args.pf_data)     # type: pd.DataFrame
    df.reset_index(inplace=True)
    df = df[df["GENOME_TYPE"] == "group-a"].copy()

    df["RE"] = df[["RBS_MAT", "NON_MAT"]].apply(
        lambda r: relative_entropy(
            MotifModel(r["RBS_MAT"], None),
            GMS2Noncoding(r["NON_MAT"])
        ), axis=1
    )

    sns.jointplot(df, "GC", "RE")
    sns.kdeplot(df, "GC", "RE")


if __name__ == "__main__":
    main(my_env, parsed_args)
