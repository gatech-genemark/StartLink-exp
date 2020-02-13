# Karl Gemayel
# Georgia Institute of Technology
#
# Created: 2/13/20

import logging
import argparse
import os

import pandas as pd
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
from sbsp_viz.general import FigureOptions
from sbsp_viz.scatter import scatter

parser = argparse.ArgumentParser("Description of driver.")

parser.add_argument('--pf-data', required=True, help="CSV data file")

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


def compare_distance_local_vs_global(env, df, **kwargs):
    # type: (Environment, pd.DataFrame, Dict[str, Any]) -> None
    pd_work = env["pd-work"]

    pf_distance = os.path.join(pd_work, "distance_local_vs_global.pdf")
    pf_alignment_length = os.path.join(pd_work, "alignment_length_local_vs_global.pdf")
    pf_ungapped_alignment_length = os.path.join(pd_work, "ungapped_alignment_length_local_vs_global.pdf")

    pf_diff_distance_vs_ratio_length = os.path.join(pd_work, "diff_distance_vs_ratio_length.pdf")
    pf_diff_distance_vs_ratio_ungapped_length = os.path.join(pd_work, "diff_distance_vs_ratio_ungapped_length.pdf")

    # compare kimura local vs global
    scatter(df, "global_distance", "local_distance", figure_options=FigureOptions(
        title="Distance by local vs global alignment",
        xlabel="Global",
        ylabel="Local",
        save_fig=pf_distance,
        balanced=True
    ))

    # compare alignment length of local vs global
    scatter(df, "global_length", "local_length", figure_options=FigureOptions(
        title="Alignment length of local vs global",
        xlabel="Global",
        ylabel="Local",
        save_fig=pf_alignment_length,
        balanced=True
    ))

    # compare ungapped alignment length of local vs global

    scatter(df, "global_length_without_gaps", "local_length_without_gaps", figure_options=FigureOptions(
        title="Ungapped alignment length of local vs global",
        xlabel="Global",
        ylabel="Local",
        save_fig=pf_ungapped_alignment_length,
        balanced=True
    ))

    # compare difference in alignment length versus difference in local/global
    df["diff_distance"] = df["global_distance"] - df["local_distance"]

    df["ratio_ungapped_length"] = df["local_length_without_gaps"] / df["global_length_without_gaps"]
    df["ratio_length"] = df["local_length"] / df["global_length"]

    scatter(df, "ratio_length", "diff_distance", figure_options=FigureOptions(
        title="Difference in distance vs ratio of alignment lengths",
        xlabel="Ratio of lengths",
        ylabel="Difference in distance",
        save_fig=pf_diff_distance_vs_ratio_length,
    ))

    scatter(df, "ratio_ungapped_length", "diff_distance", figure_options=FigureOptions(
        title="Difference in distance vs ratio of ungapped alignment lengths",
        xlabel="Ratio of ungapped lengths",
        ylabel="Difference in distance",
        save_fig=pf_diff_distance_vs_ratio_ungapped_length,
    ))


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None

    df = pd.read_csv(args.pf_data)
    compare_distance_local_vs_global(env, df)


if __name__ == "__main__":
    main(my_env, parsed_args)
