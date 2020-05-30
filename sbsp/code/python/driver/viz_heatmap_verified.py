# Karl Gemayel
# Georgia Institute of Technology
#
# Created: 3/17/20
import logging
import argparse
import math

import pandas as pd
import matplotlib.pyplot as plt
from typing import *

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import sbsp_log  # runs init in sbsp_log and configures logger

# Custom imports
from sbsp_general import Environment
from sbsp_general.shelf import next_name
from sbsp_viz.colormap import ColorMap as CM
import sbsp_viz.sns as sns
# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from sbsp_general.general import os_join

parser = argparse.ArgumentParser("Description of driver.")

parser.add_argument('--pf-data', required=True, help="Path to accuracy file")

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

def plot_fix_min_move_max(df):
    genomes = sorted(set(df["Genome"]))

    df = df[df["Min"] == 0.1]

    df = df.sort_values("Max", axis=0)
    fig, axes = plt.subplots(2, math.ceil(len(genomes) / 2), sharex="all", sharey="all")
    axes = axes.ravel()
    lines = list()
    for i in range(len(genomes)):
        name = genomes[i]
        ax = axes[i]
        df_tmp = df[df["Genome"] == name]

        s = ax.plot(df_tmp["Max"], df_tmp["Sensitivity"], label="Sensitivity")
        s = ax.plot(df_tmp["Max"], df_tmp["Coverage"], label="Coverage")
        ax.set_title(r"\textit{{{}. {}}}".format(name[0], name.split()[1]), style="italic")
        ax.set_ylim([49,101])
        lines.append(s)
        if i % 2 == 0:
            ax.set_ylabel("Percentage")
        if i >= 2:
            ax.set_xlabel("Maximum Kimura")

    plt.subplots_adjust(bottom=0.17)
    fig.legend(loc="lower center", labels=["Accuracy", "Coverage"], ncol=2)
    fig.suptitle("SBSP Performance for Kimura thresholds [0.1, x]")
    plt.savefig(next_name(my_env["pd-work"]))
    plt.show()


def plot_fix_max_move_min(df):
    genomes = sorted(set(df["Genome"]))

    df = df[df["Max"] == 0.5]

    df = df.sort_values("Min", axis=0)
    fig, axes = plt.subplots(2, math.ceil(len(genomes) / 2), sharex="all", sharey="all")
    axes = axes.ravel()
    lines = list()
    for i in range(len(genomes)):
        name = genomes[i]
        ax = axes[i]
        df_tmp = df[df["Genome"] == name]

        s = ax.plot(df_tmp["Min"], df_tmp["Sensitivity"], label="Sensitivity")
        s = ax.plot(df_tmp["Min"], df_tmp["Coverage"], label="Coverage")
        ax.set_title(r"\textit{{{}. {}}}".format(name[0], name.split()[1]), style="italic")
        ax.set_ylim([49, 101])
        lines.append(s)
        if i % 2 == 0:
            ax.set_ylabel("Percentage")
        if i >=2:
            ax.set_xlabel("Minimum Kimura")

    plt.subplots_adjust(bottom=0.17)
    fig.suptitle("SBSP Performance for Kimura thresholds [x, 0.5]")
    fig.legend(loc="lower center", labels=["Accuracy", "Coverage"], ncol=2)
    plt.savefig(next_name(my_env["pd-work"]))
    plt.show()

def plot_move_consecutive_blocks(df):
    # type: (pd.DataFrame) -> None

    genomes = sorted(set(df["Genome"]))
    fig, axes = plt.subplots(2, math.ceil(len(genomes) / 2), sharex="all", sharey="all")
    axes = axes.ravel()

    all_kimura_values = sorted(set(df["Max"]).union(set(df["Min"])))

    for i in range(len(genomes)):
        name = genomes[i]
        ax = axes[i]

        # filter df only by those with consecutive block pair
        list_df = list()
        for j in range(1, len(all_kimura_values)):
            low = all_kimura_values[j-1]
            high = all_kimura_values[j]

            df_tmp = df[(df["Min"] == low) & (df["Max"] == high)]
            list_df.append(df_tmp)

        df_tmp = pd.concat(list_df, sort=False)
        df_tmp["Average"] = (df_tmp["Max"] + df_tmp["Min"]) / 2.0

        df_tmp = df_tmp[df_tmp["Genome"] == name]

        s = ax.plot(df_tmp["Average"], df_tmp["Sensitivity"], label="Sensitivity")
        s = ax.plot(df_tmp["Average"], df_tmp["Coverage"], label="Coverage")
        ax.set_title(r"\textit{{{}. {}}}".format(name[0], name.split()[1]), style="italic")
        ax.set_ylim([49, 101])
        if i % 2 == 0:
            ax.set_ylabel("Percentage")
        if i >= 2:
            ax.set_xlabel("Average Kimura")

    plt.subplots_adjust(bottom=0.17)
    fig.suptitle("SBSP Performance for small blocks of Kimura")
    fig.legend(loc="lower center", labels=["Accuracy", "Coverage"], ncol=2)
    plt.savefig(next_name(my_env["pd-work"]))
    plt.show()




def main(env, args):
    # type: (Environment, argparse.Namespace) -> None
    df = pd.read_csv(args.pf_data)
    df["Coverage"] = 100 * df["Found"] / df["Verified"]
    import seaborn


    genomes = sorted(set(df["Genome"]))
    plt.rcParams['figure.autolayout'] = False
    fig, axes = plt.subplots(2, math.ceil(len(genomes)/2), sharex="all", sharey="all")
    axes = axes.ravel()
    for i in range(len(genomes)):
        name = genomes[i]
        ax = axes[i]
        df_tmp = df[df["Genome"] == name]
        s = ax.scatter(df_tmp["Min"], df_tmp["Max"], c=df_tmp["Sensitivity"])
        ax.set_xlim([-0.1, 0.85])
        ax.set_ylim([-0.1, 0.85])

        # s = seaborn.scatterplot("Min", "Max", "Sensitivity", data=df[df["Genome"] == name], ax=ax)
    # sns.scatterplot(df[df["Genome"] == "Escherichia coli"], "Min", "Max", hue="Sensitivity",
    #                 )
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(s, cbar_ax)
    plt.show()

    plot_fix_min_move_max(df)
    plot_fix_max_move_min(df)
    plot_move_consecutive_blocks(df)

if __name__ == "__main__":
    main(my_env, parsed_args)
