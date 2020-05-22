# Karl Gemayel
# Georgia Institute of Technology
#
# Created: 3/24/20
import ast
import collections
import logging
import argparse
import math

import numpy as np
import pandas as pd
from typing import *

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import sbsp_log  # runs init in sbsp_log and configures logger

# Custom imports
from sbsp_general import Environment
from sbsp_general.general import os_join
from sbsp_general.shelf import next_name
from sbsp_io.general import mkdir_p
from sbsp_viz.general import FigureOptions, save_figure
from sbsp_viz import sns
from sbsp_viz.colormap import ColorMap as CM

from plotnine import *

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #

parser = argparse.ArgumentParser("Visualize the stats at a genome level.")

parser.add_argument('--pf-data', required=True, help="Input data generated from analysis_per_query.py")

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

def get_summary_for_gcfid(df):
    # type: (pd.DataFrame) -> Dict[str, Any]
    result = dict()

    # expand


    first_row = df.iloc[0]

    result["GCFID"] = first_row["GCFID"]
    result["Genome GC"] = first_row["Genome GC"]
    result["Name"] = first_row["Name"]
    result["Ancestor"] = first_row["Ancestor"]



    for x in ["SBSP", "GMS2", "NCBI", "Prodigal", "GMS2=SBSP", "GMS2=SBSP=NCBI",
              "GMS2=SBSP=Prodigal", "(GMS2=SBSP)!=(NCBI=Prodigal)", "(GMS2=SBSP)!=NCBI", "(GMS2=SBSP)!=Prodigal"]:
        result[x] = df[x].sum()

    numerator_denominator = [
        ("GMS2=SBSP", "SBSP"),
        ("GMS2=SBSP", "GMS2"),
        ("GMS2=SBSP=NCBI", "GMS2=SBSP"),
        ("GMS2=SBSP=Prodigal", "GMS2=SBSP"),
        ("(GMS2=SBSP)!=NCBI", "GMS2=SBSP"),
        ("(GMS2=SBSP)!=Prodigal", "GMS2=SBSP"),
        ("(GMS2=SBSP)!=(NCBI=Prodigal)", "GMS2=SBSP")
    ]

    for pair in numerator_denominator:
        x, y = pair
        df_subset = df[df[y]]
        result["{} % {}".format(x, y)] = round(100 * df_subset[x].sum() / float(df_subset[y].sum()))

    result["Sen(SBSP,NCBI)"] = 100 * df["SBSP=NCBI"].sum() / float(((df["SBSP"]) & (df["NCBI"])).sum())
    result["Cov(SBSP,NCBI)"] = 100 * df["SBSP"].sum() / float(df.iloc[0]["Total SBSP"])
    result["Cov2(SBSP,NCBI)"] = 100 * ((df["SBSP"]) & (df["NCBI"])).sum() / float(df["NCBI"].sum())

    result["Sen(GMS2,NCBI)"] = 100 * df["GMS2=NCBI"].sum() / float(((df["GMS2"]) & (df["NCBI"])).sum())
    result["Cov(GMS2,NCBI)"] = 100 * df["GMS2"].sum() / float(df.iloc[0]["Total GMS2"])
    result["Cov2(GMS2,NCBI)"] = 100 * ((df["GMS2"]) & (df["NCBI"])).sum() / float(df["NCBI"].sum())

    result["Sen(GMS2=SBSP,NCBI)"] = 100 * df["GMS2=SBSP=NCBI"].sum() / float(((df["GMS2=SBSP"]) & (df["NCBI"])).sum())
    result["Cov(GMS2=SBSP,NCBI)"] = 100 * df["GMS2=SBSP"].sum() / float(df.iloc[0]["Total GMS2=SBSP"])
    result["Cov2(GMS2=SBSP,NCBI)"] = 100 * ((df["GMS2=SBSP"]) & (df["NCBI"])).sum() / float(df["NCBI"].sum())


    return result


def get_summary_per_gcfid(df):
    # type: (pd.DataFrame) -> pd.DataFrame

    result = list()
    for gcfid, df_group in df.groupby("GCFID", as_index=False):
        result.append(get_summary_for_gcfid(df_group))

    return pd.DataFrame(result)


def viz_summary_per_gcfid_per_step(env, df):
    # type: (Environment, pd.DataFrame) -> None

    # gather analysis for steps A, A+B, and A+B+C
    list_df = list()        # type: List[pd.DataFrame]

    # compute total number of predictions per tool, per genome
    for gcfid, df_group in df.groupby("GCFID", as_index=False):
        df.loc[df_group.index, "Total SBSP"] = df.loc[df_group.index, "SBSP"].sum()
        df.loc[df_group.index, "Total GMS2"] = df.loc[df_group.index, "GMS2"].sum()
        df.loc[df_group.index, "Total GMS2=SBSP"] = df.loc[df_group.index, "GMS2=SBSP"].sum()

    # loop over steps A, A+B, and A+B+C and collect stats
    tag = None
    for step in ["A", "B", "C"]:
        if tag is None:
            tag = step
        else:
            tag += "+" + step
        df_summary_per_gcfid = get_summary_per_gcfid(df[df["Predicted-at-step"] <= step])
        df_summary_per_gcfid["SBSP Step"] = tag
        list_df.append(df_summary_per_gcfid)

    df_per_gcfid_per_step = pd.concat(list_df, sort=False)

    import matplotlib.pyplot as plt
    # fig, ax = plt.subplots()
    #
    # sns.lineplot(df_per_gcfid_per_step, "SBSP Step", "SBSP", hue="GCFID", ax=ax,
    #              sns_kwargs={"palette": CM.get_map("verified")},
    #              legend=False
    #              )
    # for l in ax.lines:
    #     l.set_linestyle("--")
    #
    # ax2 = ax.twinx()
    # sns.lineplot(df_per_gcfid_per_step, "SBSP Step", "Sen(SBSP,NCBI)", hue="GCFID", ax=ax2,
    #              sns_kwargs={"palette": CM.get_map("verified")},)
    #
    # fo = FigureOptions(
    #     xlabel="SBSP Step",
    #     ylabel="Percentage",
    #     # ylim=[0, 105],
    #     save_fig=next_name(env["pd-work"])
    # )
    # FigureOptions.set_properties_for_axis(ax, fo)
    # plt.subplots_adjust(bottom=0.2)
    # handles, labels = ax.get_legend_handles_labels()
    # ax.legend(handles=handles[1:], labels=labels[1:],
    #           loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.25))
    #
    # plt.savefig(fo.save_fig)
    # plt.show()


    fig, axes = plt.subplots(3, 2, sharex="all", sharey="row")
    ax = axes[:, 0]

    sns.lineplot(df_per_gcfid_per_step, "SBSP Step", "Sen(SBSP,NCBI)", hue="GCFID", ax=ax[0],
                 sns_kwargs={"palette": CM.get_map("verified")}, legend=False,
                 figure_options=FigureOptions(
            ylabel="Sensitivity",
            ylim=[85,105],

        ))

    sns.lineplot(df_per_gcfid_per_step, "SBSP Step", "Cov(SBSP,NCBI)", hue="GCFID", ax=ax[1],
                 sns_kwargs={"palette": CM.get_map("verified")},
                 legend=False, figure_options=FigureOptions(
            ylabel="Percent of Genes",
            ylim=[0, None]
        )
                 )


    sns.lineplot(df_per_gcfid_per_step, "SBSP Step", "SBSP", hue="GCFID", ax=ax[2],
                 sns_kwargs={"palette": CM.get_map("verified")},
                 legend=False, figure_options=FigureOptions(
            ylabel="Number of Genes",
            ylim=[0, None]
        )
                 )

    fig.align_ylabels(ax)

    # plt.savefig(next_name(env["pd-work"]))
    # plt.show()

    # fig, ax = plt.subplots(3, 1, sharex="all")
    ax = axes[:, 1]
    sns.lineplot(df_per_gcfid_per_step, "SBSP Step", "Sen(GMS2=SBSP,NCBI)", hue="GCFID", ax=ax[0],
                 sns_kwargs={"palette": CM.get_map("verified")}, legend=False,
                 figure_options=FigureOptions(
                     ylabel="Sensitivity",
                     ylim=[85, 105],

                 ))

    sns.lineplot(df_per_gcfid_per_step, "SBSP Step", "Cov(GMS2=SBSP,NCBI)", hue="GCFID", ax=ax[1],
                 sns_kwargs={"palette": CM.get_map("verified")},
                 legend=False, figure_options=FigureOptions(
            ylabel="Percent of Genes",
            ylim=[0, None]
        )
                 )

    sns.lineplot(df_per_gcfid_per_step, "SBSP Step", "GMS2=SBSP", hue="GCFID", ax=ax[2],
                 sns_kwargs={"palette": CM.get_map("verified")},
                figure_options=FigureOptions(
            ylabel="Number of Genes",
            ylim=[0, None]
        )
                 )

    ax[2].get_legend().remove()

    fig.align_ylabels(ax)

    for ax in axes.ravel():
        ax.set_xlabel("Steps")

    axes[0][0].set_title("SBSP")
    axes[0][1].set_title("GMS2=SBSP")

    fig.subplots_adjust(bottom=0.21)

    # handles, labels = ax.get_legend_handles_labels()
    # fig.legend(handles=handles[1:], labels=labels[1:], loc="lower center", ncol=4)#, bbox_to_anchor=(0.5, -0.25))
    handles, labels = ax.get_legend_handles_labels()
    labels[0]="Genome"
    fig.legend(handles=handles, labels=labels, loc="lower center", ncol=3)#, bbox_to_anchor=(0.5, -0.25))
    plt.savefig(next_name(env["pd-work"]))
    plt.show()

    # three plots

    for gcfid, df_group in df.groupby("GCFID", as_index=False):
        df.loc[df_group.index, "Total SBSP"] = ((df_group["SBSP"]) & (df_group["NCBI"])).sum()
        df.loc[df_group.index, "Total GMS2"] = ((df_group["GMS2"]) & (df_group["NCBI"])).sum()
        df.loc[df_group.index, "Total GMS2=SBSP"] = ((df_group["GMS2=SBSP"]) & (df_group["NCBI"])).sum()

    df_all = get_summary_per_gcfid(df)

    # map column names for tables
    columns = ["GCFID", "NCBI", "Sen(SBSP,NCBI)", "Sen(GMS2,NCBI)", "Sen(GMS2=SBSP,NCBI)",
               "Cov2(SBSP,NCBI)", "Cov2(GMS2,NCBI)", "Cov2(GMS2=SBSP,NCBI)"]
    df_sen = df_all.copy()[columns].rename(columns={
        "GCFID": "Genome",
        "NCBI": "Verified",
        "Sen(SBSP,NCBI)": "SBSP",
        "Sen(GMS2,NCBI)": "GMS2",
        "Sen(GMS2=SBSP,NCBI)": "GMS2=SBSP",
    }, inplace=False)
    df_sen[["Genome", "Verified", "SBSP", "GMS2", "GMS2=SBSP"]].to_csv(
        os_join(env["pd-work"], "sensitivity.csv"), index=False)

    # print(df_all[["GCFID", "NCBI", "Cov2(SBSP,NCBI)", "Cov2(GMS2,NCBI)", "Cov2(GMS2=SBSP,NCBI)"]].to_string(index=False))

    df_cov = df_all[columns].rename(columns={
        "GCFID": "Genome",
        "NCBI": "Verified",
        "Cov2(SBSP,NCBI)": "SBSP",
        "Cov2(GMS2,NCBI)": "GMS2",
        "Cov2(GMS2=SBSP,NCBI)": "GMS2=SBSP",
    }, inplace=False)

    df_cov[["Genome", "Verified", "SBSP", "GMS2", "GMS2=SBSP"]].to_csv(
        os_join(env["pd-work"], "coverage.csv"), index=False)


def viz_analysis_per_query(env, df, **kwargs):
    # type: (Environment, pd.DataFrame, Dict[str, Any]) -> None

    pd_work = env["pd-work"]

    df["(GMS2=SBSP)!=NCBI"] = df["GMS2=SBSP"] & df["NCBI"] & ~df["GMS2=SBSP=NCBI"]
    df["(GMS2=SBSP)!=Prodigal"] = df["GMS2=SBSP"] & df["Prodigal"] & ~df["GMS2=SBSP=Prodigal"]

    # print(df[["Ancestor", "GMS2=SBSP"]].groupby("Ancestor", as_index=False).sum().to_string())
    # print(df[["Ancestor", "GCFID"]].groupby("Ancestor")["GCFID"].nunique().to_string())

    viz_summary_per_gcfid_per_step(env, df)


def read_analysis_per_query_to_df(pf_data):
    # type: (str) -> pd.DataFrame
    return pd.read_csv(pf_data)


def fix_names(r):
    # type: (pd.Series) -> str
    return "{}. {}".format(
        r["GCFID"][0], r["GCFID"].split("_")[1]
    )


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None

    df = read_analysis_per_query_to_df(args.pf_data)
    df["GCFID"] = df.apply(fix_names, axis=1)       # Change names from director format to human readable

    viz_analysis_per_query(env, df)


if __name__ == "__main__":
    main(my_env, parsed_args)
