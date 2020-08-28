# Karl Gemayel
# Georgia Institute of Technology
#
# Created: 3/17/20
import logging
import argparse
import pandas as pd
from typing import *
import seaborn
import matplotlib.pyplot as plt

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

parser = argparse.ArgumentParser("Visualization additional statistics collected from log files.")

parser.add_argument('--pf-stats', required=True, help="File containing stats from log files.")

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

def viz_per_genome(env, df):
    # type: (Environment, pd.DataFrame) -> None

    df_grp = df.groupby(["Genome", "Ancestor"], as_index=False).mean()

    sns.catplot(df_grp, "Ancestor", "BLAST", figure_options=FigureOptions(
        save_fig=next_name(env["pd-work"]),
        xlabel="Clade",
        ylabel="Number of BLASTp Hits"
    ),
                sns_kwargs={"palette": CM.get_map("ancestor")}
    )

    # list_grp = list()
    # for _, df_grp in df.groupby("Genome", as_index=False):
    #     indices = df_grp.index
    #
    #     list_grp.append({
    #         "Genome": df.at[indices[0], "Genome"],
    #         "Ancestor": df.at[indices[0], "Ancestor"],
    #         "= 0": len(df_grp[df_grp["BLAST"] == 0]),
    #         **{
    #             f"< {x}": len(df_grp[df_grp["BLAST"] < x]) for x in [5, 10, 20, 50, 100, 500, 1000, 5000, 10000]
    #         },
    #         "> 10000": len(df_grp[df_grp["BLAST"] > 10000])
    #     })
    #
    # df_grp = pd.DataFrame(list_grp)
    # sns.catplot(df_grp, "Ancestor", "= 0")
    # sns.catplot(df_grp, "Ancestor", "< 5")
    # sns.catplot(df_grp, "Ancestor", "< 50")
    # sns.catplot(df_grp, "Ancestor", "< 100")

    # plots
    # 1) x: number of queries with < x targets

    # compute per genome, the % of queries with hits <= 0, 5, 10, 20, 40, 80, 160, ... 240 580 1160, ...
    # plot

    list_entries = list()
    for _, df_grp in df.groupby("Genome", as_index=False):
        indices = df_grp.index
        genome = df.at[indices[0], "Genome"]
        ancestor = df.at[indices[0], "Ancestor"]

        total_queries = len(df_grp)
        curr = 0
        for n in range(40):

            list_entries.append({
                "Genome": genome,
                "Ancestor": ancestor,
                "x": curr,
                "y": 100 * len(df_grp[df_grp["BLAST"] < curr]) / total_queries
            })

            # if list_entries[-1]["y"] == 100:
            #     break



            if curr == 0:
                curr = 5
            else:
                curr *= 1.2


    df_tmp = pd.DataFrame(list_entries)

    sns.lineplot(df_tmp, "x", "y", hue="Ancestor", figure_options=FigureOptions(
        xlabel="Number of BLASTp hits",
        ylabel="Cumulative percentage of queries (per genome)",
        save_fig=next_name(env["pd-work"]),
    ),
                 legend_loc="best",
                 legend_title="",
                 legend_ncol=2,
                 sns_kwargs={"ci": "sd", "palette": CM.get_map("ancestor")})



    sns.lineplot(df_tmp, "y", "x", hue="Ancestor", figure_options=FigureOptions(
        ylabel="Number of BLASTp hits",
        xlabel="Cumulative percentage of queries (per genome)",
        save_fig=next_name(env["pd-work"]),
    ),
                 legend_loc="best",
                 legend_title="",
                 legend_ncol=2,
                 sns_kwargs={"ci": "sd", "palette": CM.get_map("ancestor")})

    fig, axes = plt.subplots(2, 2, sharex="all", sharey="all")

    ancestors = sorted(set(df["Ancestor"]))

    for anc, ax in zip(ancestors, axes.ravel()):

        df_anc = df_tmp[df_tmp["Ancestor"] == anc]
        sns.lineplot(df_anc[df_anc["x"] <= 40], "x", "y", hue="Ancestor",
                     legend=None,
                     ax=ax,
                     sns_kwargs={"ci": "sd", "palette": CM.get_map("ancestor")})
        ax.set_title(anc)
        ax.set_xlabel(None)
        ax.set_ylabel(None)

    figure_options = FigureOptions(
        xlabel="Number of BLASTp hits",
        ylabel="Cumulative percentage of queries (per genome)",
        save_fig=next_name(env["pd-work"]),
    )

    fig.add_subplot(111, frameon=False)
    # # hide tick and tick label of the big axes
    plt.tick_params(top=False, bottom=False, left=False, right=False, which="both",
                    labelbottom=False, labeltop=False, labelleft=False, labelright=False)
    plt.xlabel(figure_options.xlabel, labelpad=20)
    plt.ylabel(figure_options.ylabel, labelpad=30)

    save_figure(figure_options, fig)

    plt.show()

def compute_more(df):
    # type: (pd.DataFrame) -> None

    df["Filtered"] = df["BLAST"] - df["QuickFilter"]
    df["Filtered %"] = 100 * df["Filtered"] / df["BLAST"]


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None
    df = pd.read_csv(args.pf_stats)

    compute_more(df)
    fo = FigureOptions(ylim=[0, 700000])

    viz_per_genome(env, df)

    # sns.distplot(df, "BLAST", sns_kwargs={"kde": False}, figure_options=fo)
    # sns.distplot(df[df["BLAST"] < 10000], "BLAST", sns_kwargs={"kde": False}, figure_options=fo)
    # sns.distplot(df[df["BLAST"] < 2000], "BLAST", sns_kwargs={"kde": False}, figure_options=fo)
    #
    # sns.distplot(df,"QuickFilter", sns_kwargs={"kde": False}, figure_options=fo)
    # sns.distplot(df[df["BLAST"] < 10000], "BLAST", sns_kwargs={"kde": False}, figure_options=fo)
    # sns.distplot(df[df["BLAST"] < 10000], "QuickFilter", sns_kwargs={"kde": False}, figure_options=fo)
    # #
    # # # sns.kdeplot(df, "BLAST", "QuickFilter")
    # # sns.scatterplot(df, "BLAST", "QuickFilter", identity=True,
    # #                 sns_kwargs={"alpha": 0.3, "marker": "o", "s": 0.7})
    # sns.distplot(df[df["BLAST"] > 2000], "Filtered %")
    #
    #
    #
    # sns.kdeplot(df, "BLAST", None, hue="Ancestor")



if __name__ == "__main__":
    main(my_env, parsed_args)
