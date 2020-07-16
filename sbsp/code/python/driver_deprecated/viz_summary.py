# Karl Gemayel
# Georgia Institute of Technology
#
# Created: 3/5/20

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

parser = argparse.ArgumentParser("Description of driver.")

parser.add_argument('--pf-summary', required=True, help="Summary statistics file")

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
    steps = ["A", "B", "C", "U"]
    for s in steps:
        df["{}: GMS2=SBSP %".format(s)] = 100 * df["{}: GMS2=SBSP".format(s)] / df["GMS2=SBSP"]

    for s in steps:
        df["{}: GMS2=SBSP=NCBI %".format(s)] = 100 * df["{}: GMS2=SBSP=NCBI".format(s)] / df["{}: GMS2=SBSP".format(s)]
        df["{}: (GMS2=SBSP)!=NCBI".format(s)] = df["{}: GMS2=SBSP".format(s)] - df["{}: GMS2=SBSP=NCBI".format(s)]

    df["GMS2=SBSP %"] = 100 * df["GMS2=SBSP"] / df["GMS2"]
    df["GMS2=SBSP=NCBI %"] = 100 * df["GMS2=SBSP=NCBI"] / df["GMS2=SBSP"]

    df["GMS2=SBSP % From SBSP"] = 100 * df["GMS2=SBSP"] / df["SBSP"]
    df["GMS2=SBSP % From NCBI"] = 100 * df["GMS2=SBSP"] / df["NCBI"]

    df["SBSP % From NCBI"] = 100 * df["SBSP"] / df["NCBI"]

    df["(GMS2=SBSP)!=NCBI"] = df["GMS2=SBSP"] - df["GMS2=SBSP=NCBI"]

    df["(GMS2=SBSP)!=NCBI % From GMS2=SBSP"] = 100 * df["(GMS2=SBSP)!=NCBI"] / df["GMS2=SBSP"]


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

    df = pd.read_csv(args.pf_summary, header=0)
    add_percentages(df)
    steps = ["A", "B", "C", "U"]

    import seaborn as sns
    import matplotlib.pyplot as plt
    sns.set_context(context="paper", font_scale=1.5)

    colors = ["windows blue", "amber", "faded green", "dusty purple"]
    palette = sns.xkcd_palette(colors)
    sns.palplot(palette)
    plt.show()
    sns.set_palette(palette)


    # plot steps
    fig_num = 0
    for col in ["{}: GMS2=SBSP %", "{}: GMS2=SBSP=NCBI %", "{}: (GMS2=SBSP)!=NCBI"]:

        df_stacked = stack_columns_as_rows(
            df,
            [col.format(s) for s in steps],
            col.format("M"),
            steps,
            "Step"
        )

        for kind in ["bar", "box"]: #["point", "bar", "strip", "swarm", "box", "violin", "boxen"]:
        # for kind in ["point", "bar", "strip", "swarm", "box", "violin", "boxen"]:
            # plt.figure(figsize=(12, 4))

            g = sns.catplot(x="Ancestor", y=col.format("M"), hue="Step", kind=kind, data=df_stacked,
                    dodge=True, legend=False, aspect=2)
            g.set(ylabel=col.split(maxsplit=1)[1])

            plt.legend(loc='center right', bbox_to_anchor=(1.125, 0.5))
            plt.savefig(os.path.join(env["pd-work"], "{}.pdf".format(fig_num)), bbox_inches='tight')
            plt.show()

            fig_num += 1

    common={"linewidth": 0}

    plt.figure(figsize=(12,4))
    g = sns.scatterplot(x="GC", y="GMS2=SBSP=NCBI %", hue="Ancestor", data=df, **common)
    g.set(ylim=[50, 100])
    plt.legend(loc='center left',  bbox_to_anchor=(1.05, 0.5))
    plt.savefig(os.path.join(env["pd-work"], "{}.pdf".format(fig_num)), bbox_inches='tight')
    plt.show()
    fig_num += 1


    plt.figure(figsize=(12,4))
    g = sns.scatterplot(x="GC", y="GMS2=SBSP % From SBSP", hue="Ancestor", data=df, **common)
    g.set(ylim=[50, 100])
    plt.legend(loc='center left', bbox_to_anchor=(1.05, 0.5))
    plt.savefig(os.path.join(env["pd-work"], "{}.pdf".format(fig_num)), bbox_inches='tight')
    plt.show()
    fig_num += 1

    plt.figure(figsize=(12,4))
    g = sns.scatterplot(x="GC", y="(GMS2=SBSP)!=NCBI", hue="Ancestor", data=df, **common)
    plt.legend(loc='center left', bbox_to_anchor=(1.05, 0.5))
    plt.savefig(os.path.join(env["pd-work"], "{}.pdf".format(fig_num)), bbox_inches='tight')
    plt.show()
    fig_num += 1

    plt.figure()
    g = sns.catplot(x="Ancestor", y="(GMS2=SBSP)!=NCBI", data=df, kind="box", aspect=1.5)
    # g.set(ylim=[50, 100])
    plt.savefig(os.path.join(env["pd-work"], "{}.pdf".format(fig_num)), bbox_inches='tight')
    plt.show()
    fig_num += 1


    plt.figure(figsize=(12,4))
    g = sns.scatterplot(x="GC", y="(GMS2=SBSP)!=NCBI % From GMS2=SBSP", hue="Ancestor", data=df, **common)
    g.set(ylim=(0,None))
    g.set(ylabel="(GMS2=SBSP)!=NCBI %\nFrom GMS2=SBSP")
    plt.legend(loc='center left', bbox_to_anchor=(1.05, 0.5))
    plt.savefig(os.path.join(env["pd-work"], "{}.pdf".format(fig_num)), bbox_inches='tight')
    plt.show()
    fig_num += 1

    plt.figure()
    g = sns.catplot(x="Ancestor", y="(GMS2=SBSP)!=NCBI % From GMS2=SBSP", data=df, kind="box", aspect=1.5)
    g.set(ylabel="(GMS2=SBSP)!=NCBI %\nFrom GMS2=SBSP")
    # g.set(ylim=[50, 100])
    plt.savefig(os.path.join(env["pd-work"], "{}.pdf".format(fig_num)), bbox_inches='tight')
    plt.show()
    fig_num += 1

    plt.figure(figsize=(12,4))
    g = sns.scatterplot(x="GC", y="GMS2=SBSP %", hue="Ancestor", data=df, **common)
    g.set(ylim=[50, 100])
    plt.legend(loc='center left', bbox_to_anchor=(1.05, 0.5))
    plt.savefig(os.path.join(env["pd-work"], "{}.pdf".format(fig_num)), bbox_inches='tight')
    plt.show()
    fig_num += 1

    plt.figure(figsize=(12, 4))
    g = sns.scatterplot(x="GC", y="GMS2=SBSP % From NCBI", hue="Ancestor", data=df, **common)
    g.set(ylim=[50, 100])
    plt.legend(loc='center left', bbox_to_anchor=(1.05, 0.5))
    plt.savefig(os.path.join(env["pd-work"], "{}.pdf".format(fig_num)), bbox_inches='tight')
    plt.show()
    fig_num += 1

    plt.figure(figsize=(12, 4))
    g = sns.scatterplot(x="GC", y="GMS2=SBSP", hue="Ancestor", data=df, **common)
    plt.legend(loc='center left', bbox_to_anchor=(1.05, 0.5))
    plt.savefig(os.path.join(env["pd-work"], "{}.pdf".format(fig_num)), bbox_inches='tight')
    plt.show()
    fig_num += 1

    plt.figure()
    g = sns.catplot(x="Ancestor", y="GMS2=SBSP %", data=df, kind="box", aspect=1.5)
    # g.set(ylim=[50, 100])
    plt.savefig(os.path.join(env["pd-work"], "{}.pdf".format(fig_num)), bbox_inches='tight')
    plt.show()
    fig_num += 1

    plt.figure()
    g = sns.catplot(x="Ancestor", y="GMS2=SBSP % From SBSP", data=df, kind="box", aspect=1.5)
    g.set(ylim=[50, 100])
    plt.savefig(os.path.join(env["pd-work"], "{}.pdf".format(fig_num)), bbox_inches='tight')
    plt.show()
    fig_num += 1

    g = sns.lmplot(x="GC", y="GMS2=SBSP % From SBSP", hue="Ancestor", data=df, aspect=2, legend=False, ci=None)
    g.set(ylim=[50, 100])
    plt.legend(bbox_to_anchor=(1.05, 0.75), loc=2, borderaxespad=0.)
    plt.savefig(os.path.join(env["pd-work"], "{}.pdf".format(fig_num)), bbox_inches='tight')
    plt.show()
    fig_num += 1


    # number of genes

    plt.figure(figsize=(6, 3))

    g = sns.scatterplot(x="GC", y="GMS2=SBSP", hue="Ancestor", data=df, **common)
    plt.savefig(os.path.join(env["pd-work"], "{}.pdf".format(fig_num)), bbox_inches='tight')
    plt.show()
    fig_num += 1

    # number of genes
    plt.figure(figsize=(6, 3))
    g = sns.catplot(x="Ancestor", y="GMS2=SBSP", kind="box", data=df, legend=False, aspect=1.5)
    plt.savefig(os.path.join(env["pd-work"], "{}.pdf".format(fig_num)), bbox_inches='tight')
    plt.show()
    fig_num += 1


    plt.figure(figsize=(6, 3))
    g = sns.catplot(x="Ancestor", y="SBSP % From NCBI", kind="box", data=df, aspect=1.5)
    g.set(ylim=[50, 100])
    plt.savefig(os.path.join(env["pd-work"], "{}.pdf".format(fig_num)), bbox_inches='tight')
    plt.show()
    fig_num += 1

    plt.figure(figsize=(12, 4))
    g = sns.scatterplot(x="GC", y="SBSP % From NCBI", hue="Ancestor", data=df, **common)
    plt.legend(loc='center left', bbox_to_anchor=(1.05, 0.5))
    g.set(ylim=[50, 100])
    plt.savefig(os.path.join(env["pd-work"], "{}.pdf".format(fig_num)), bbox_inches='tight')
    plt.show()
    fig_num += 1

    plt.figure(figsize=(12, 4))

    g = sns.lmplot(x="GC", y="SBSP % From NCBI", hue="Ancestor", data=df, aspect=2, legend=False, ci=None)
    plt.legend(bbox_to_anchor=(1.05, 0.75), loc=2, borderaxespad=0.)
    g.set(ylim=[50, 100])
    plt.savefig(os.path.join(env["pd-work"], "{}.pdf".format(fig_num)), bbox_inches='tight')
    plt.show()
    fig_num += 1

    # number of genomes
    plt.figure(figsize=(6, 3))
    g = sns.countplot(x="Ancestor", data=df )
    g.set(ylabel="Number of genomes")
    plt.savefig(os.path.join(env["pd-work"], "{}.pdf".format(fig_num)), bbox_inches='tight')
    plt.show()
    fig_num += 1

    # overlap consistency per genome
    plt.figure(figsize=(6, 3))
    g = sns.distplot(df["Overlap Consistency"], bins=50)
    g.set(ylabel="Number of genomes")
    plt.savefig(os.path.join(env["pd-work"], "{}.pdf".format(fig_num)), bbox_inches='tight')
    plt.show()
    fig_num += 1

    # overlap consistency per gene

    pc_columns = ["Overlap Consistency List"] + [x for x in df.columns if "RPC" in x]

    for col in pc_columns:
        df_tmp = convert_overlap_consistency_list_to_one_per_row(df, col)

        plt.figure(figsize=(6, 3))
        g = sns.distplot(df_tmp["Consistency"], bins=50, kde=False)
        g.set(ylabel="Number of genes")
        g.set(title=col)
        plt.savefig(os.path.join(env["pd-work"], "{}.pdf".format(fig_num)), bbox_inches='tight')
        plt.show()
        fig_num += 1
        #
        # # overap per gene per ancestor
        # plt.figure(figsize=(6, 3))
        # unique_vals = sorted(set(df_tmp['Ancestor']))
        # targets = [df_tmp.loc[df_tmp['Ancestor'] == val] for val in unique_vals]
        #
        # for target in targets:
        #     sns.distplot(target[['Consistency']], kde=False)
        #
        # g.set(ylabel="Number of genes")
        # plt.savefig(os.path.join(env["pd-work"], "{}.pdf".format(fig_num)), bbox_inches='tight')
        # plt.show()
        # fig_num += 1
        #
        # # boxplot
        # plt.figure(figsize=(12, 4))
        #
        # g = sns.catplot(x="Ancestor", y="Consistency", data=df_tmp, kind="box", aspect=2)
        # g.set(ylabel="Number of genes")
        # plt.savefig(os.path.join(env["pd-work"], "{}.pdf".format(fig_num)), bbox_inches='tight')
        # plt.show()
        # fig_num += 1


    # boxplot
    plt.figure(figsize=(12, 4))

    g = sns.catplot(x="Ancestor", y="Percent of Overlap Groups", data=df, kind="box", aspect=2)
    g.set(ylabel="Percent of components")
    g.set(title="Percent of components with most common upstream distance < 3nt")
    plt.savefig(os.path.join(env["pd-work"], "{}.pdf".format(fig_num)), bbox_inches='tight')
    plt.show()
    fig_num += 1

    # boxplot
    plt.figure(figsize=(12, 4))

    g = sns.lmplot(x="GC", y="Percent of Overlap Groups", hue="Ancestor", data=df, aspect=2, legend=False, ci=None)
    g.set(ylabel="Percent of components")
    g.set(title="Percent of components with most common upstream distance < 3nt")
    plt.legend(bbox_to_anchor=(1.05, 0.75), loc=2, borderaxespad=0.)

    plt.savefig(os.path.join(env["pd-work"], "{}.pdf".format(fig_num)), bbox_inches='tight')
    plt.show()
    fig_num += 1


    print(df["Ancestor"].value_counts())
    print(df[["Ancestor", "GMS2=SBSP"]].groupby("Ancestor").agg("sum"))





    print("Done")

    pass


if __name__ == "__main__":
    main(my_env, parsed_args)
