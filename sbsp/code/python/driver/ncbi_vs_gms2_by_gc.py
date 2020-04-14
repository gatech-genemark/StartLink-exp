# Karl Gemayel
# Georgia Institute of Technology
#
# Created: 3/17/20
import logging
import argparse
import pandas as pd
from typing import *

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import sbsp_log  # runs init in sbsp_log and configures logger

# Custom imports
from sbsp_general import Environment
import sbsp_viz.sns as sns

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from sbsp_general.shelf import next_name
from sbsp_viz.general import FigureOptions

parser = argparse.ArgumentParser("Plot the 5' end different between NCBI and GMS2 as GC increases.")

parser.add_argument('--pf-data', required=True)

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

    df = pd.read_csv(args.pf_data)
    df["Err(GMS2,NCBI)"] = 100 - 100 * df["Identical"] / df["Found"]
    df = df[df["Err(GMS2,NCBI)"] < 30].copy()
    df.loc[df["Group"] == "D2", "Group"] = "D"

    sns.lmplot(df, "GC", "Err(GMS2,NCBI)", hue="Type", figure_options=FigureOptions(
        xlabel="Genome GC",
        save_fig=next_name(env["pd-work"]),
        ylim=[0, None],
    ),
               legend_loc="best",
               sns_kwargs={"scatter_kws": {"s": 5, "alpha": 0.3}, "lowess": True, "scatter": True,
                           "aspect": 1.5}
               )

    sns.lmplot(df, "GC", "Err(GMS2,NCBI)", hue="Group", figure_options=FigureOptions(
        xlabel="Genome GC",
        save_fig=next_name(env["pd-work"]),
        ylim=[0, None],
    ),
               legend_loc="best",
               sns_kwargs={"scatter_kws": {"s": 5, "alpha": 0.3}, "lowess": True, "scatter": True,
                           "aspect": 1.5}
               )


if __name__ == "__main__":
    main(my_env, parsed_args)
