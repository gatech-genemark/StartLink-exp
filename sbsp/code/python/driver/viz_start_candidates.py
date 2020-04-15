# Karl Gemayel
# Georgia Institute of Technology
#
# Created: 3/5/20

import logging
import argparse
import os

import pandas as pd
from typing import *

import seaborn as sns
import matplotlib.pyplot as plt

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import sbsp_log  # runs init in sbsp_log and configures logger

# Custom imports
from sbsp_general import Environment
from sbsp_viz.colormap import ColorMap as CM

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #

parser = argparse.ArgumentParser("Description of driver.")

parser.add_argument('--pf-input', required=True, help="Summary statistics file")

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


def stack_columns_as_rows(df, stack, new_col, labels, label_col="class"):
    # type: (pd.DataFrame, List[str], str, List[str], str) -> pd.DataFrame

    if labels is not None:
        if len(labels) != len(stack):
            raise ValueError("Labels and stack must have same number of elements: {} != {}".format(
                len(labels), len(stack))
            )
    else:
        labels = stack


    list_df = list()

    for c, l in zip(stack, labels):
        if c not in df.columns:
            logger.warning("Cannot find column {} in data frame.".format(c))
            continue

        df_curr = df.copy()
        df_curr.rename(columns={c: new_col}, inplace=True)            # new column name
        df_curr[label_col] = l

        list_df.append(df_curr)

    return pd.concat(list_df, axis=0, sort=False)


def add_percentages(df):
    # type: (pd.DataFrame) -> None

    # fraction predicted by step
    df["Total Candidates"] = df["downstream"] + df["upstream"]



def convert_overlap_consistency_list_to_one_per_row(df, col):
    # type: (pd.DataFrame, str) -> pd.DataFrame

    list_rows = list()

    for index, row in df.iterrows():

        ancestor = row["Ancestor"]
        str_list = row[col]
        if len(str_list.strip("[]").strip()) == 0:
            continue

        l = [float(i) for i in str_list.strip("[]").split(",")]
        for c in l:

            list_rows.append({"Ancestor": ancestor, "Consistency": c})

    return pd.DataFrame(list_rows)


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None

    df = pd.read_csv(args.pf_input, header=0)
    add_percentages(df)
    sns.set_context(context="paper", font_scale=1.5)

    #clean up
    df = df[df["Total Candidates"] < 100]

    colors = ["windows blue", "amber", "faded green", "dusty purple"]
    palette = sns.xkcd_palette(colors)
    sns.palplot(palette)
    plt.show()
    sns.set_palette(palette)

    fig_num = 0

    plt.figure(figsize=(12, 4))
    sns.jointplot(x="gc", y="Total Candidates", data=df)
    # plt.legend(loc='center left', bbox_to_anchor=(1.05, 0.5))
    plt.savefig(os.path.join(env["pd-work"], "{}.pdf".format(fig_num)), bbox_inches='tight')
    plt.show()
    fig_num += 1

    # Average number of candidates per GC
    df_tmp = df.groupby(["gcfid", "ancestor"], as_index=False).agg("mean")
    plt.figure(figsize=(12, 4))
    g = sns.scatterplot(x="gc", y="Total Candidates", data=df_tmp, hue="ancestor", palette=CM.get_map("ancestor"))
    plt.legend(loc='center left', bbox_to_anchor=(1.05, 0.5))
    g.set(ylabel="Average number of candidates")
    g.set(xlabel="GC")
    plt.savefig(os.path.join(env["pd-work"], "{}.pdf".format(fig_num)), bbox_inches='tight')
    plt.show()
    fig_num += 1

    # Average number of candidates per GC
    plt.figure(figsize=(12, 4))
    g = sns.lmplot(x="gc", y="Total Candidates", data=df_tmp, hue="ancestor", aspect=2, legend=False, ci=None,
                   lowess=True, palette=CM.get_map("ancestor"))
    plt.legend(loc='center left', bbox_to_anchor=(1.05, 0.5))
    g.set(ylabel="Average number of candidates")
    g.set(xlabel="GC")
    plt.savefig(os.path.join(env["pd-work"], "{}.pdf".format(fig_num)), bbox_inches='tight')
    plt.show()
    fig_num += 1

    pass


if __name__ == "__main__":
    main(my_env, parsed_args)
