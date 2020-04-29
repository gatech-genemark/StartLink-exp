# Karl Gemayel
# Georgia Institute of Technology
#
# Created: 3/17/20
import logging
import argparse
import math

import pandas as pd
from typing import *

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import sbsp_log  # runs init in sbsp_log and configures logger

# Custom imports
from sbsp_general import Environment
import sbsp_viz.sns as sns
from sbsp_general.shelf import next_name
from sbsp_viz.colormap import ColorMap as CM

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from sbsp_viz.general import FigureOptions, save_figure

parser = argparse.ArgumentParser("Description of driver.")

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


def plot_for_block(df):
    # type: (pd.DataFrame) -> None
    import matplotlib.pyplot as plt
    import seaborn

    genomes = sorted(set(df["Genome"]))

    num_genome = len(genomes)

    fig, axes = plt.subplots(2, math.ceil(num_genome / 2), sharey="all", sharex="all")
    axes = axes.ravel()
    for i, g in enumerate(genomes):
        ax = axes[i]

        for x, df_group in df[df["Genome"] == g].groupby("region"):
            seaborn.kdeplot(df_group["score"], label=x, ax=ax)

        ax.set_title(g)

    plt.legend()
    plt.show()


def plot_for_5prime(df, bins=20):
    # type: (pd.DataFrame) -> None
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    import seaborn

    genomes = sorted(set(df["Genome"]))
    regions = sorted(set(df["region"]))
    num_genome = len(genomes)
    import numpy as np

    fig, axes = plt.subplots(2, math.ceil(num_genome/2), sharey="all", sharex="all")
    axes = axes.ravel()
    for i, g in enumerate(genomes):
        ax = axes[i]

        for x in regions:
            df_group = df[(df["Genome"] == g) & (df["region"] == x)]
        # for x, df_group in df[df["Genome"] == g].groupby("region"):
        #     seaborn.distplot(df_group["score"], label=x if i == 0 else None, ax=ax,
        #                      norm_hist=False, kde=False, bins=30)

            hist_values, bin_edges = np.histogram(df_group["score"], bins=bins)
            hist_values = 100* hist_values / sum(hist_values)

            # bin_edges = bin_edges[:len(bin_edges)-1]
            bin_edges = bin_edges[1:]
            seaborn.lineplot(bin_edges, hist_values, markers=False, ax=ax, label=x if i == 0 else None, legend=False)

            # ax.hist(df_group["score"], bins=30, normed=True, histtype="step")
            # seaborn.barplot()

        ax.set_title(r"\textit{{{}}}".format(str(g)))
        ax.set_xlabel(None)
        ax.set_ylim([0, None])
        y_max = ax.get_ylim()[1]
        ax.axvline(0.5, 0, y_max, color="grey", linestyle="dashed")

    plt.subplots_adjust(bottom=0.17)
    fig.legend(loc="lower center", ncol=len(regions))

    fig.add_subplot(111,frameon=False)
    plt.tick_params(top=False, bottom=False, left=False, right=False, which="both",
                    labelbottom=False, labeltop=False, labelleft=False, labelright=False)
    plt.xlabel("Score",labelpad=25)
    plt.ylabel("Percentage per group",labelpad=30)

    # figure_options = FigureOptions(
    #     xlabel="Score", ylabel="Frequency", save_fig=next_name(my_env["pd-work"])
    # )
    plt.savefig(next_name(my_env["pd-work"]))
    plt.show()


def fix_names(r):
    # type: (pd.Series) -> None
    return "{}. {}".format(
        r["Genome"][0], r["Genome"].split("_")[1]
    )


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None
    df = pd.read_csv(args.pf_data)
    df["Genome"] = df.apply(fix_names, axis=1)

    df[df["score"] < 0] = 0

    df.replace({"region": {"Dc": "Down close", "Df": "Down far", "Uc": "Up close", "Uf": "Up far"}}, inplace=True)

    df_block = df[df["type"] == "block"]
    df_5prime = df[df["type"] == "5prime"]

    # sns.kdeplot(df, "score", None, hue="region")

    plot_for_5prime(df_5prime, bins=40)
    plot_for_5prime(df_block, bins=10)



if __name__ == "__main__":
    main(my_env, parsed_args)
