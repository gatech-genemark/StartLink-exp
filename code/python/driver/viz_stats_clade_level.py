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
from sbsp_viz import sns
from sbsp_viz.colormap import ColorMap as CM

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from sbsp_viz.general import FigureOptions, save_figure

parser = argparse.ArgumentParser("Visualize statistics by clade (for large scale analysis).")

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

TOOL = "StartLink"
TOOLp = TOOL + "+"

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
        ("GMS2=SBSP", "NCBI"),
        ("SBSP", "NCBI"),
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

    return result


def get_summary_per_gcfid(df):
    # type: (pd.DataFrame) -> pd.DataFrame

    result = list()
    for gcfid, df_group in df.groupby("GCFID", as_index=False):
        result.append(get_summary_for_gcfid(df_group))

    return pd.DataFrame(result)


def viz_summary_per_gcfid(env, df, title=None):
    # type: (Environment, pd.DataFrame) -> None
    pd_work = env['pd-work']
    sns.catplot(
        df, "Ancestor", "GMS2=SBSP % SBSP", kind="box",
        figure_options=FigureOptions(
            save_fig=next_name(pd_work),
            ylim=[None, 100],
            title=title,
            ylabel=f"Percentage of {TOOLp} from {TOOL} Predictions"
        ),
        sns_kwargs={"palette": CM.get_map("ancestor")}
    )

    sns.catplot(df, "Ancestor", "(GMS2=SBSP)!=NCBI % GMS2=SBSP", kind="box", figure_options=FigureOptions(
        save_fig=next_name(pd_work),
        ylim=[0, 20],
        # ylabel="1 - Sen(NCBI, GMS2=SBSP)",
        ylabel=f"Diff(PGAP, {TOOLp})",
        xlabel="Clade",
        title=title
    ),
                sns_kwargs={"palette": CM.get_map("ancestor")})

    # per GC
    # sns.scatterplot(df, "Genome GC", "(GMS2=SBSP)!=NCBI % GMS2=SBSP", hue="Ancestor", figure_options=FigureOptions(
    #     save_fig=next_name(pd_work),
    #     ylim=[0, None],
    #     title=title,
    # ),
    #                 legend_loc="best",
    #                 sns_kwargs={"palette": CM.get_map("ancestor")})

    # per GC
    # sns.lmplot(df, "Genome GC", "(GMS2=SBSP)!=NCBI % GMS2=SBSP", hue="Ancestor", figure_options=FigureOptions(
    #     save_fig=next_name(pd_work),
    #     ylim=[0, None],
    #     title=title,
    #     ylabel="Err(NCBI, GMS2=SBSP)"
    # # ylabel="1 - Sen(NCBI, GMS2=SBSP)",
    # ),
    #
    #            sns_kwargs={"palette": CM.get_map("ancestor"), "scatter": False, "lowess": True})

    sns.lmplot(df, "Genome GC", "(GMS2=SBSP)!=NCBI % GMS2=SBSP", hue="Ancestor", figure_options=FigureOptions(
        save_fig=next_name(pd_work),
        ylim=[0, None],
        title=title,
        # ylabel="1 - Sen(NCBI, GMS2=SBSP)",
        ylabel=f"Diff(PGAP, {TOOLp})"
    ),
               legend_loc="best",
               sns_kwargs={"palette": CM.get_map("ancestor"), "scatter": True, "lowess": True,
                           "scatter_kws": {"s": 5}, "aspect": 1.5
                           })

    sns.lmplot(df, "Genome GC", "GMS2=SBSP", hue="Ancestor", figure_options=FigureOptions(
        save_fig=next_name(pd_work),
        ylim=[0, None],
        title=title,
    ),
               sns_kwargs={"palette": CM.get_map("ancestor"), "scatter": True, "lowess": True,
                           "scatter_kws": {"s": 5}
                           })

    sns.lmplot(df, "Genome GC", "GMS2=SBSP % SBSP", hue="Ancestor", figure_options=FigureOptions(
        save_fig=next_name(pd_work),
        ylim=[50, 100],
        title=title,
    ),
               sns_kwargs={"palette": CM.get_map("ancestor"), "scatter": True, "lowess": True,
                           "scatter_kws": {"s": 5}
                           })

    sns.lmplot(df, "Genome GC", "GMS2=SBSP % NCBI", hue="Ancestor", figure_options=FigureOptions(
        save_fig=next_name(pd_work),
        ylim=[50, 100],
        title=title,
        xlabel="GC",
        ylabel=f"Cov({TOOLp}, PGAP)"
    ),
               sns_kwargs={"palette": CM.get_map("ancestor"), "scatter": True, "lowess": True,
                           "scatter_kws": {"s": 5}
                           })

    sns.lmplot(df, "Genome GC", "SBSP % NCBI", hue="Ancestor", figure_options=FigureOptions(
        save_fig=next_name(pd_work),
        ylim=[50, 100],
        title=title,
        xlabel="GC",
        ylabel=f"Cov({TOOL}, PGAP)"
    ),
               sns_kwargs={"palette": CM.get_map("ancestor"), "scatter": True, "lowess": True,
                           "scatter_kws": {"s": 5}
                           })

    sns.catplot(df, "Ancestor", "GMS2=SBSP % NCBI", kind="box", figure_options=FigureOptions(
        save_fig=next_name(pd_work),
        # ylim=[0, 20],
        ylim=[50, 100],
        # ylabel="1 - Sen(NCBI, GMS2=SBSP)",
        ylabel=f"Cov({TOOLp}, PGAP)",
        xlabel="Clade",
        title=title
    ),
                sns_kwargs={"palette": CM.get_map("ancestor")})

    sns.catplot(df, "Ancestor", "SBSP % NCBI", kind="box", figure_options=FigureOptions(
        save_fig=next_name(pd_work),
        # ylim=[0, 20],
        ylim=[50, 100],
        # ylabel="1 - Sen(NCBI, GMS2=SBSP)",
        ylabel=f"Cov({TOOL}, PGAP)",
        xlabel="Clade",
        title=title
    ),
                sns_kwargs={"palette": CM.get_map("ancestor")})

    # sns.lmplot(df, "Genome GC", "GMS2=SBSP % GMS2", hue="Ancestor", figure_options=FigureOptions(
    #     save_fig=next_name(pd_work),
    #     ylim=[50, 100],
    #     title=title,
    # ),
    #            sns_kwargs={"palette": CM.get_map("ancestor"), "scatter": True, "lowess": True,
    #                        "scatter_kws": {"s": 5}
    #                        })

    # sns.scatterplot(df, "NCBI", "(GMS2=SBSP)!=NCBI % GMS2=SBSP", hue="Ancestor", figure_options=FigureOptions(
    #     save_fig=next_name(pd_work),
    #     ylim=[0, None],
    #     title=title,
    # ),
    #            sns_kwargs={"palette": CM.get_map("ancestor"),
    #                        })

    # sns.scatterplot(df, "GMS2=SBSP", "(GMS2=SBSP)!=NCBI % GMS2=SBSP", hue="Ancestor", figure_options=FigureOptions(
    #     save_fig=next_name(pd_work),
    #     ylim=[0, None],
    #     title=title,
    # ),
    #                 sns_kwargs={"palette": CM.get_map("ancestor"),
    #                             })

    # # per GC
    # sns.scatterplot(df, "Genome GC", "(GMS2=SBSP)!=Prodigal % GMS2=SBSP", hue="Ancestor", figure_options=FigureOptions(
    #     save_fig=next_name(pd_work),
    #     ylim=[0, None],
    #     title=title,
    # ),
    #                 sns_kwargs={"palette": CM.get_map("ancestor")})


def contour_kimura_per_ancestor(env, df):
    import seaborn
    import matplotlib.pyplot as plt

    ancestors = sorted(list(set(df["Ancestor"])))
    fig, axes = plt.subplots(2, math.ceil(len(ancestors) / 2), sharex=True, sharey=True, figsize=(6, 6))

    for anc, ax in zip(ancestors, axes.ravel()):
        df_group = df[df["Ancestor"] == anc]
        seaborn.kdeplot(df_group["Min-Kimura"].values, df_group["Max-Kimura"].values, ax=ax)
        ax.set_title(anc)
        # ax.set_ylim([0.45, 0.525])

    # fig.xlabel("Min-Kimura")
    # plt.xlabel("Min-Kimura")
    # plt.ylabel("Max-Kimura")
    # fig.text(0.5, 0.04, 'Min-Kimura', ha='center')
    # fig.text(0.04, 0.5, 'Max-Kimura', va='center', rotation='vertical')
    fig.add_subplot(111, frameon=False)
    # # hide tick and tick label of the big axes
    plt.tick_params(top=False, bottom=False, left=False, right=False, which="both",
                    labelbottom=False, labeltop=False, labelleft=False, labelright=False)
    plt.xlabel("Minimum Kimura distance", labelpad=20)
    plt.ylabel("Maximum Kimura distance", labelpad=30)

    fig.tight_layout()
    save_figure(FigureOptions(
        save_fig=next_name(env["pd-work"])
    ))

    plt.show()


def kimura_dist_plot(env, df):
    import seaborn
    import matplotlib.pyplot as plt

    ancestors = list(set(df["Ancestor"]))
    # fig, axes = plt.subplots(2, math.ceil(len(ancestors)/2), sharex=True, sharey=True)
    #
    # for anc, ax in zip(ancestors, axes.ravel()):
    #
    #     df_group = df[df["Ancestor"] == anc]
    #     seaborn.distplot(df_group["Average-Kimura"], ax=ax, color=CM.get_map("ancestor")[anc],
    #                      hist=False)
    #     ax.set_title(anc)
    # plt.show()

    fig, ax = plt.subplots()  # type: plt.Figure, plt.Axes
    for anc in ancestors:
        df_group = df[df["Ancestor"] == anc]
        seaborn.distplot(df_group["Average-Kimura"], ax=ax, color=CM.get_map("ancestor")[anc],
                         hist=False, label=anc)
        # ax.set_title(anc)

    ax.legend(ancestors)
    ax.set_ylabel("PDF")
    save_figure(FigureOptions(
        save_fig=next_name(env["pd-work"])
    ))
    plt.show()


def heat_map_Kimura_accuracy(env, df_all, x, y, num_steps=20, balance=False,
                             xlabel=None, ylabel=None):
    # type: (Environment, pd.DataFrame, str, str, int) -> None
    import matplotlib.pyplot as plt

    if xlabel is None:
        xlabel = x
    if ylabel is None:
        ylabel = y
    ancestors = sorted(list(set(df_all["Ancestor"])))
    fig, axes = plt.subplots(2, math.ceil(len(ancestors) / 2), sharex=True, sharey=True)
    cbar_ax = fig.add_axes([.91, .3, .03, .4])

    # fig = plt.figure()
    num_rows = 2
    num_cols = math.ceil(len(ancestors) / 2)

    axis_idx = 0
    curr_row = 0
    curr_col = 0
    for ancestor, df in df_all.groupby("Ancestor", as_index=False):
        ax = axes.ravel()[axis_idx]
        # ax = plt.subplot2grid((num_rows, num_cols), (curr_row, curr_col))
        axis_idx += 1
        curr_col += 1
        if curr_col == math.ceil(len(ancestors) / 2):
            curr_row += 1
            curr_col = 0

        min_x = min(df[x])
        max_x = max(df[x]) + 0.000000001

        min_y = min(df[y])
        max_y = max(df[y]) + 0.000000001

        if balance:
            min_x = min_y = min(min_x, min_y)
            max_x = max_y = max(max_x, max_y)

        ss_x = (max_x - min_x) / float(num_steps)
        ss_y = (max_y - min_y) / float(num_steps)

        num_col = num_steps
        num_row = num_steps
        import numpy as np
        gms2_eq_sbsp_and_ncbi = np.zeros([num_row, num_col], dtype=float)
        gms2_eq_sbsp_eq_ncbi = np.zeros([num_row, num_col], dtype=float)

        df_gms2_eq_sbsp_and_ncbi = (df["GMS2=SBSP"]) & (df["NCBI"])
        df_gms2_eq_sbsp_eq_ncbi = (df["GMS2=SBSP=NCBI"])

        for index in df.index:
            x_val = df.at[index, x]
            y_val = df.at[index, y]

            x_pos = int((x_val - min_x) / ss_x)
            y_pos = int((y_val - min_y) / ss_y)

            gms2_eq_sbsp_and_ncbi[x_pos][y_pos] += 1 if df.at[index, "GMS2=SBSP"] and df.at[index, "NCBI"] else 0
            gms2_eq_sbsp_eq_ncbi[x_pos][y_pos] += 1 if df.at[index, "GMS2=SBSP=NCBI"] else 0

        gms2_eq_sbsp_and_ncbi[gms2_eq_sbsp_and_ncbi < 10] = 0
        accuracy = 100 - 100 * np.divide(gms2_eq_sbsp_eq_ncbi, gms2_eq_sbsp_and_ncbi)
        # accuracy = np.flip(accuracy, 0)

        import seaborn
        import matplotlib.pyplot as plt

        xticks = list(range(0, num_steps, int(num_steps / 5)))
        yticks = list(range(0, num_steps, int(num_steps / 5)))

        l_x = np.arange(min_x, max_x, ss_x)
        l_y = np.arange(min_y, max_y, ss_y)
        xticklabels = [round(l_x[i], 2) for i in xticks]
        yticklabels = [round(l_y[i], 2) for i in yticks]
        g = seaborn.heatmap(
            accuracy.transpose(), vmin=0, vmax=40,
            xticklabels=xticklabels, yticklabels=yticklabels, ax=ax,
            # cmap=seaborn.light_palette("green"),
            cmap="Blues",
            cbar=False)
        # cbar_ax=None if axis_idx != 0 else cbar_ax, cbar=axis_idx==0)

        # cbar=g.cbar

        g.invert_yaxis()
        g.set_xticks(xticks)
        g.set_yticks(yticks)
        g.set_xticklabels(xticklabels, rotation=0)

        # g.set_xlabel("Min Kimura")
        # g.set_ylabel("Max Kimura")
        g.set_title(ancestor)
        mappable = ax.collections[0]

    # im = plt.gca().get_children()[0]
    # cax = fig.add_axes([0.8, 0.1, 0.03, 0.8])
    cbar_ax = fig.axes[-1]

    # fig.tight_layout(rect=[0, 0, .9, 1])
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axes
    plt.tick_params(top=False, bottom=False, left=False, right=False, which="both",
                    labelbottom=False, labeltop=False, labelleft=False, labelright=False)
    plt.xlabel(xlabel, labelpad=20)
    plt.ylabel(ylabel, labelpad=30)

    # ax3 = plt.subplot2grid((num_rows, num_cols), (0, num_cols - 1), rowspan=num_rows,
    #                        )
    plt.colorbar(mappable, cax=cbar_ax)
    fig.tight_layout(rect=[0, 0, .9, 1])

    save_figure(FigureOptions(
        save_fig=next_name(env["pd-work"])
    ))

    plt.show()


def bin_data_one_d(env, df_all, feature, num_steps=20):
    # type: (Environment, pd.DataFrame, str, int) -> pd.DataFrame
    min_x = min(df_all[feature])
    max_x = max(df_all[feature]) + 0.000000001
    ss_x = (max_x - min_x) / float(num_steps)

    list_df = list()
    for ancestor, df in df_all.groupby("Ancestor", as_index=False):

        import numpy as np
        gms2_eq_sbsp_and_ncbi = np.zeros(num_steps, dtype=float)
        gms2_eq_sbsp_eq_ncbi = np.zeros(num_steps, dtype=float)

        for index in df.index:
            x_val = df.at[index, feature]

            x_pos = int((x_val - min_x) / ss_x)

            gms2_eq_sbsp_and_ncbi[x_pos] += 1 if df.at[index, "GMS2=SBSP"] and df.at[index, "NCBI"] else 0
            gms2_eq_sbsp_eq_ncbi[x_pos] += 1 if df.at[index, "GMS2=SBSP=NCBI"] else 0

        accuracy = np.divide(gms2_eq_sbsp_eq_ncbi, gms2_eq_sbsp_and_ncbi)

        xticks = list(range(0, num_steps))

        l_x = np.arange(min_x, max_x, ss_x)
        xticklabels = [round(l_x[i], 2) for i in xticks]

        curr_df = pd.DataFrame({
            feature: xticklabels,
            "Accuracy": accuracy,
            "Number-of-queries": gms2_eq_sbsp_and_ncbi
        })
        curr_df["Ancestor"] = ancestor
        list_df.append(curr_df)

    return pd.concat(list_df)  # type: pd.DataFrame


def one_dim_Kimura_accuracy(env, df_all, num_steps=20):
    # type: (Environment, pd.DataFrame, int) -> None
    import matplotlib.pyplot as plt
    pd_work = env["pd-work"]
    ancestors = sorted(list(set(df_all["Ancestor"])))
    # fig, axes = plt.subplots(2, math.ceil(len(ancestors) / 2), sharex=True, sharey=True)

    # min_x = min(df_all["Average-Kimura"])
    # max_x = max(df_all["Average-Kimura"]) + 0.000000001
    # ss_x = (max_x - min_x) / float(num_steps)
    #
    # list_df = list()
    # axis_idx = 0
    # for ancestor, df in df_all.groupby("Ancestor", as_index=False):
    #     # ax = axes.ravel()[axis_idx]
    #     # axis_idx += 1
    #
    #
    #
    #
    #
    #     import numpy as np
    #     gms2_eq_sbsp_and_ncbi = np.zeros(num_steps, dtype=float)
    #     gms2_eq_sbsp_eq_ncbi = np.zeros(num_steps, dtype=float)
    #
    #     df_gms2_eq_sbsp_and_ncbi = (df["GMS2=SBSP"]) & (df["NCBI"])
    #     df_gms2_eq_sbsp_eq_ncbi = (df["GMS2=SBSP=NCBI"])
    #
    #     for index in df.index:
    #
    #         x_val = df.at[index, "Average-Kimura"]
    #
    #         x_pos = int((x_val-min_x) / ss_x)
    #
    #         gms2_eq_sbsp_and_ncbi[x_pos] += 1 if df.at[index, "GMS2=SBSP"] and df.at[index, "NCBI"] else 0
    #         gms2_eq_sbsp_eq_ncbi[x_pos] += 1 if df.at[index, "GMS2=SBSP=NCBI"] else 0
    #
    #     accuracy = np.divide(gms2_eq_sbsp_eq_ncbi, gms2_eq_sbsp_and_ncbi)
    #     # accuracy = np.flip(accuracy, 0)
    #
    #
    #     xticks = list(range(0, num_steps))
    #
    #     l_x = np.arange(min_x, max_x, ss_x)
    #     xticklabels = [round(l_x[i], 2) for i in xticks]
    #     # g = seaborn.heatmap(accuracy.transpose(), vmin=0, vmax=1, xticklabels=xticklabels, yticklabels=yticklabels, ax=ax,
    #     #                     cbar=True)
    #
    #     # g = seaborn.lineplot(xticklabels, accuracy, ax=ax, label=ancestor)
    #
    #     # cbar=g.cbar
    #
    #     # g.set_xticks(xticks)
    #
    #     curr_df = pd.DataFrame({
    #         "Average-Kimura": xticklabels,
    #         "Accuracy": accuracy,
    #         "Number-of-queries": gms2_eq_sbsp_and_ncbi
    #     })
    #     curr_df["Ancestor"] = ancestor
    #     list_df.append(curr_df)
    #
    #     # g.set_xlabel("Min Kimura")
    #     # g.set_ylabel("Max Kimura")
    #     # g.set_title(ancestor)
    #
    # df = pd.concat(list_df)     # type: pd.DataFrame
    df = bin_data_one_d(env, df_all, "Average-Kimura", num_steps)
    sns.lineplot(df, "Average-Kimura", "Accuracy", hue="Ancestor", figure_options=FigureOptions(
        save_fig=next_name(pd_work),
    ), sns_kwargs={"palette": CM.get_map("ancestor")}
                 )
    sns.lineplot(df, "Average-Kimura", "Number-of-queries", hue="Ancestor", figure_options=FigureOptions(
        save_fig=next_name(pd_work),
    ), sns_kwargs={"palette": CM.get_map("ancestor")}
                 )

    total_per_ancestor = {
        ancestor: (df["Ancestor"].isin({ancestor})).sum() for ancestor in ancestors
    }

    df["Percentage-of-queries"] = 0
    df["Cumulative-percentage-of-queries"] = 0
    df.reset_index(inplace=True)
    for ancestor, df_group in df.groupby("Ancestor", as_index=False):  # type: str, pd.DataFrame
        df_group.sort_values("Average-Kimura", inplace=True)
        index = df_group.index

        prev = 0
        total = df_group["Number-of-queries"].sum()
        df.loc[index, "Percentage-of-queries"] = 100 * df.loc[index, "Number-of-queries"] / float(total)

        for i in index:
            df.loc[i, "Cumulative-percentage-of-queries"] = prev + df.loc[i, "Percentage-of-queries"]
            prev = df.loc[i, "Cumulative-percentage-of-queries"]

    fig, ax = plt.subplots(figsize=(8, 4))
    sns.lineplot(df, "Average-Kimura", "Percentage-of-queries", hue="Ancestor", figure_options=FigureOptions(
        save_fig=next_name(pd_work),
        ylabel="Percentage of queries",
        xlabel="Average Kimura distance"
    ),
                 ax=ax,
                 show=True,
                 legend_loc="best",
                 legend_title="",
                 sns_kwargs={"palette": CM.get_map("ancestor")})

    fig, ax = plt.subplots(figsize=(8, 4))
    import seaborn
    seaborn.lineplot("Average-Kimura", "Percentage-of-queries", data=df, hue="Ancestor", ax=ax,
                     palette=CM.get_map("ancestor"))

    ax.set_ylabel("Percentage of queries")
    ax.set_xlabel("Average Kimura distance")
    plt.legend(loc="best", title="", frameon=False)
    # plt.tight_layout()
    # plt.subplots_adjust(bottom=-0.1)
    # plt.savefig(next_name(pd_work))
    fig.savefig(next_name(pd_work), bbox_inches='tight')

    plt.show()

    sns.lineplot(df, "Average-Kimura", "Cumulative-percentage-of-queries", hue="Ancestor", figure_options=FigureOptions(
        save_fig=next_name(pd_work),
    ), sns_kwargs={"palette": CM.get_map("ancestor")})

    # standard dev
    df = bin_data_one_d(env, df_all[df_all["Support"] > 2], "Std-Kimura", num_steps)
    sns.lineplot(df, "Std-Kimura", "Accuracy", hue="Ancestor", figure_options=FigureOptions(
        save_fig=next_name(pd_work),
    ), sns_kwargs={"palette": CM.get_map("ancestor")}
                 )
    sns.lineplot(df, "Std-Kimura", "Number-of-queries", hue="Ancestor", figure_options=FigureOptions(
        save_fig=next_name(pd_work),
    ), sns_kwargs={"palette": CM.get_map("ancestor")}
                 )

    total_per_ancestor = {
        ancestor: (df["Ancestor"].isin({ancestor})).sum() for ancestor in ancestors
    }

    df["Percentage-of-queries"] = 0
    df["Cumulative-percentage-of-queries"] = 0
    df.reset_index(inplace=True)
    for ancestor, df_group in df.groupby("Ancestor", as_index=False):  # type: str, pd.DataFrame
        df_group.sort_values("Std-Kimura", inplace=True)
        index = df_group.index

        prev = 0
        total = df_group["Number-of-queries"].sum()
        df.loc[index, "Percentage-of-queries"] = 100 * df.loc[index, "Number-of-queries"] / float(total)

        for i in index:
            df.loc[i, "Cumulative-percentage-of-queries"] = prev + df.loc[i, "Percentage-of-queries"]
            prev = df.loc[i, "Cumulative-percentage-of-queries"]
    sns.lineplot(df, "Std-Kimura", "Percentage-of-queries", hue="Ancestor", figure_options=FigureOptions(
        save_fig=next_name(pd_work),
    ), sns_kwargs={"palette": CM.get_map("ancestor")})

    sns.lineplot(df, "Std-Kimura", "Cumulative-percentage-of-queries", hue="Ancestor", figure_options=FigureOptions(
        save_fig=next_name(pd_work),
    ), sns_kwargs={"palette": CM.get_map("ancestor")})

    # im = plt.gca().get_children()[0]
    # cax = fig.add_axes([0.8, 0.1, 0.03, 0.8])
    #
    # # fig.tight_layout(rect=[0, 0, .9, 1])
    # fig.add_subplot(111, frameon=False)
    # # hide tick and tick label of the big axes
    # plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    # plt.grid(False)
    # plt.xlabel("Min-Kimura")
    # plt.ylabel("Max-Kimura")
    #
    #
    # save_figure(FigureOptions(
    #     save_fig=next_name(env["pd-work"])
    # ))
    #
    # plt.show()


def analyze_kimura_distances(env, df):
    # type: (Environment, pd.DataFrame) -> None
    pd_work = env["pd-work"]

    df = df[df["Kimura-to-query"] != "[]"].copy()
    df["Kimura-to-query"] = df["Kimura-to-query"].apply(ast.literal_eval)
    df["Average-Kimura"] = df["Kimura-to-query"].apply(np.mean)
    df["Std-Kimura"] = df["Kimura-to-query"].apply(np.std)

    # sns.lmplot(
    #     df, "Genome GC", "Average-Kimura", hue="Ancestor", sns_kwargs={
    #         "scatter": False, "lowess": True, "scatter_kws": {"s": 5},
    #         "palette": CM.get_map("ancestor")
    #     },
    #     figure_options=FigureOptions(save_fig=next_name(pd_work))
    # )
    df_mean = df.groupby(["Ancestor", "GCFID"], as_index=False).mean()

    sns.lmplot(
        df_mean, "Genome GC", "Average-Kimura", hue="Ancestor", sns_kwargs={
            "scatter": True, "lowess": True, "scatter_kws": {"s": 5},
            "palette": CM.get_map("ancestor")
        },
        figure_options=FigureOptions(save_fig=next_name(pd_work))
    )

    # Min/max kimura
    df["Min-Kimura"] = df["Kimura-to-query"].apply(min)
    df["Max-Kimura"] = df["Kimura-to-query"].apply(max)

    # contour_kimura_per_ancestor(env, df)
    one_dim_Kimura_accuracy(env, df)

    kimura_dist_plot(env, df)
    heat_map_Kimura_accuracy(env, df, "Min-Kimura", "Max-Kimura", balance=True,
                             xlabel="Minimum Kimura distance",
                             ylabel="Maximum Kimura distance")
    heat_map_Kimura_accuracy(env, df, "Average-Kimura", "Std-Kimura", balance=False)


def most_frequent(items):
    # type: (Iterable[Any]) -> Any
    occurence_count = collections.Counter(items)
    return occurence_count.most_common(1)[0][0]


def compute_consistency(distances, pivot, flexibility=0):
    # type: (List[int], int, int) -> float

    if flexibility == 0:
        return len([1 for d in distances if d == pivot]) / float(len(distances))
    else:
        numerator = 0
        for d in distances:
            if d is not None and d <= pivot + flexibility and d >= pivot - flexibility:
                numerator += 1

        return float(numerator) / len(distances)


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
        df_curr.rename(columns={c: new_col}, inplace=True)  # new column name
        df_curr[label_col] = l

        df_curr.drop(set(labels).difference({l}), inplace=True, axis=1)

        list_df.append(df_curr)

    return pd.concat(list_df, axis=0, sort=True)


def analyze_upstream_distances(env, df):
    # type: (Environment, pd.DataFrame) -> None
    pd_work = os_join(env["pd-work"], "upstream_distances")
    mkdir_p(pd_work)

    # remove empty lists
    df = df[df["Upstream-distance"] != "[]"].copy()
    df["Upstream-distance"] = df["Upstream-distance"].apply(ast.literal_eval)
    df["Most frequent upstream"] = df["Upstream-distance"].apply(most_frequent)

    # compute consistencies with different flexibilities
    for flexibility in {0, 3}:
        df["DC(x,{})".format(flexibility)] = df[["Most frequent upstream", "Upstream-distance"]].apply(
            lambda r: compute_consistency(r["Upstream-distance"], r["Most frequent upstream"], flexibility),
            axis=1
        )

    df = df[df["Support"] > 10].copy()

    # for mf in range(-20, 50):
    #     df_mf = df[df["Most frequent upstream"] == mf]
    #     if len(df_mf) < 50:
    #         continue
    #
    #     sns.distplot(df_mf, "PC(x,0)", figure_options=FigureOptions(
    #         title="PC({},{})".format(mf, 0),
    #         save_fig=next_name(pd_work),
    #         xlim=(0,1)
    #     ))
    #     sns.distplot(df_mf, "PC(x,3)", figure_options=FigureOptions(
    #         title="PC({},{})".format(mf, 3),
    #         save_fig=next_name(pd_work),
    #         xlim=(0, 1)
    #     ))

    # plot distribution of Average PC
    import seaborn
    import matplotlib.pyplot as plt

    df_tmp = df[(df["Support"] > 10) & (df["Most frequent upstream"] < 100) & (df["Most frequent upstream"] > -50)]
    # NCBI consistency as a func
    df = df[(df["Support"] > 10) & (df["GMS2=SBSP"]) & (df["Most frequent upstream"] < 100) & (
                df["Most frequent upstream"] > -50)]

    df_tmp = stack_columns_as_rows(df_tmp[["Most frequent upstream", "DC(x,0)", "DC(x,3)", "Ancestor"]],
                                   ["DC(x,0)", "DC(x,3)"], "DC(x,f)", None, label_col="Flexibility")
    # seaborn.lmplot("Most frequent upstream", "PC(x,f)", df_tmp,
    #                scatter=False, hue="Flexibility", lowess=True)
    # plt.show()
    #
    # seaborn.lmplot("Most frequent upstream", "PC(x,f)", df_tmp,
    #             hue="Flexibility", lowess=True)
    # plt.show()
    #
    # seaborn.lmplot("Most frequent upstream", "PC(x,f)", df_tmp,
    #                scatter=False, hue="Flexibility")
    # plt.show()

    df_tmp["Most frequent upstream"] -= 1

    sns.lmplot(
        df_tmp, "Most frequent upstream", "DC(x,f)", hue="Flexibility", sns_kwargs={
            "scatter": False, "lowess": True
        },
        legend_loc="best",
        figure_options=FigureOptions(save_fig=next_name(pd_work), xlim=[-10, None], ylim=[0, 1],
                                     ylabel="Distance conservation", xlabel="Most frequent distance to upstream gene")
    )



    # sns.distplot(df, "Most frequent upstream", figure_options=FigureOptions(
    #     save_fig=next_name(pd_work)
    # ),
    #              sns_kwargs={"kde": True})

    df["Most frequent upstream"] -= 1

    fig, ax = plt.subplots()
    import seaborn
    # seaborn.countplot("Most frequent upstream", data=df[(df["Most frequent upstream"] < 10) & (df["Most frequent upstream"] > -10)], hue="Ancestor")
    (df[(df["Most frequent upstream"] < 9) & (df["Most frequent upstream"] >= -10)]
     .groupby("Ancestor")["Most frequent upstream"]
     .value_counts(normalize=True)
     .mul(100)
     .rename('Percentage (by clade)')
     .reset_index()
     .pipe((seaborn.catplot, 'data'), x="Most frequent upstream", y='Percentage (by clade)', hue="Ancestor",
           kind='point', scale=0.5, legend=False,
           order=range(-9, 10),
           palette=CM.get_map("ancestor"), aspect=1.5
           ))

    plt.legend(loc="best", title="Clade")
    figure_options = FigureOptions(
        save_fig=next_name(pd_work),
        xlabel="Most frequent distance to upstream gene",
        ylabel="Percentage of components (by clade)"
    )
    plt.xlabel(figure_options.xlabel)
    plt.ylabel(figure_options.ylabel)
    save_figure(figure_options)

    plt.show()

    fig, ax = plt.subplots()
    # scatter
    (df[(df["Most frequent upstream"] < 9) & (df["Most frequent upstream"] >= -10)]
     .groupby("Ancestor")["Most frequent upstream"]
     .value_counts(normalize=True)
     .mul(100)
     .rename('Percentage (by clade)')
     .reset_index()
     .pipe((seaborn.lineplot, 'data'), x="Most frequent upstream", y='Percentage (by clade)', hue="Ancestor",
           legend=False,
           palette=CM.get_map("ancestor"),
           ))

    plt.legend(loc="best", title="Clade")
    figure_options = FigureOptions(
        save_fig=next_name(pd_work),
        xlabel="Most frequent distance to upstream gene",
        ylabel="Percentage of components (by clade)",
        xlim=[-9.5, 9.5]
    )
    plt.xlabel(figure_options.xlabel)
    plt.ylabel(figure_options.ylabel)
    plt.xlim(figure_options.xlim)
    plt.xticks(range(-8, 9, 2))

    save_figure(figure_options)

    plt.show()

    fig, ax = plt.subplots()

    tmp_df = df[(df["Most frequent upstream"] < 10) & (df["Most frequent upstream"] > -10)].groupby("Ancestor")[
        "Most frequent upstream"].value_counts(normalize=True).mul(100).rename('Percentage (by clade)').reset_index()

    missing_values = set(range(-9, 10)).difference(set(df["Most frequent upstream"]))
    for m in missing_values:
        for a in set(tmp_df["Ancestor"]):
            tmp_df.append({"Most frequent upstream": m, "Percentage (by clade)": 0, "Ancestor": a}, ignore_index=True)

    (tmp_df.pipe((seaborn.barplot, 'data'), x="Most frequent upstream", y='Percentage (by clade)', hue="Ancestor",
                 order=range(-9, 10), palette=CM.get_map("ancestor")))

    plt.legend(loc="best", title="Clade")
    figure_options = FigureOptions(
        save_fig=next_name(pd_work),
        xlabel="Most frequent distance to upstream gene",
        ylabel="Percentage of components (by clade)",
        # xlim=[-9.5, 9.5]
    )
    plt.xlabel(figure_options.xlabel, fontsize=20)
    plt.ylabel(figure_options.ylabel, fontsize=20)
    # plt.xlim(figure_options.xlim)
    # plt.xticks(range(-8,9,2))

    # save_figure(figure_options)
    # plt.subplots_adjust(bottom=-0.1)
    # plt.savefig(next_name(pd_work))
    fig.savefig(next_name(pd_work), bbox_inches='tight')

    plt.show()

    import copy
    import numpy as np
    df_plot = copy.deepcopy(df_tmp)

    df_plot.loc[df_plot["Most frequent upstream"] == 0, "DC(x,f)"] = np.nan

    import matplotlib.gridspec as gridspec
    xlim = [-10, 100]
    fig3 = plt.figure(constrained_layout=True)
    # gs = fig3.add_gridspec(4, 1)
    gs = gridspec.GridSpec(4, 1)
    ax1 = fig3.add_subplot(gs[1:, 0])
    (df_plot.pipe((seaborn.lineplot, 'data'), x="Most frequent upstream", y='DC(x,f)', hue="Flexibility", ax=ax1))
    ax1.set_ylim(0, 1)
    ax1.set_xlim(*xlim)

    ax2 = fig3.add_subplot(gs[0, 0])
    df_t = df_plot[df_plot["Flexibility"] == "DC(x,0)"]
    (df_t[(df_t["Most frequent upstream"] <= xlim[1]) & (df_t["Most frequent upstream"] >= xlim[0])]
     .groupby("Flexibility")["Most frequent upstream"]
     .value_counts(normalize=True)
     .mul(100)
     .rename('Percentage (by clade)')
     .reset_index().pipe((seaborn.lineplot, 'data'), x="Most frequent upstream", y='Percentage (by clade)', ax=ax2))
    # (df_plot[df_plot["Flexibility"] == "DC(x,0)"].pipe((seaborn.barplot, 'data'), x="Most frequent upstream", y='DC(x,f)', ax=ax2, order=range(*xlim)))
    # seaborn.distplot(df_plot[df_plot["Flexibility"] == "DC(x,0)"]["Most frequent upstream"], ax=ax2, norm_hist=False, bins=50)
    ax2.set_xlabel("")
    # ax2.set_xticklabels([])
    # ax2.set_xlim(*xlim)

    plt.show()

    import copy
    import numpy as np
    from matplotlib.ticker import FuncFormatter

    df_plot = copy.deepcopy(df_tmp)

    df_plot.loc[df_plot["Most frequent upstream"] == -5, "DC(x,f)"] = np.nan

    def number_formatter(number, pos=None):
        """Convert a number into a human readable format."""
        magnitude = 0
        while abs(number) >= 1000:
            magnitude += 1
            number /= 1000.0
        return '%.0f%s' % (number, ['', 'K', 'M', 'B', 'T', 'Q'][magnitude])

    import matplotlib.gridspec as gridspec
    xlim = [-10, 100]
    fig3 = plt.figure(constrained_layout=True)
    # gs = fig3.add_gridspec(4, 1)
    gs = gridspec.GridSpec(4, 1, )
    ax1 = fig3.add_subplot(gs[1:, 0])
    (df_plot.pipe((seaborn.lineplot, 'data'), x="Most frequent upstream", y='DC(x,f)', hue="Flexibility", ax=ax1))
    ax1.set_ylim(0, 1)
    ax1.set_xlim(*xlim)
    ax1.set_xlabel("Most frequent distance to upstream gene")

    ax2 = fig3.add_subplot(gs[0, 0])
    df_t = df_plot[df_plot["Flexibility"] == "DC(x,0)"]
    (df_t[(df_t["Most frequent upstream"] <= xlim[1]) & (df_t["Most frequent upstream"] >= xlim[0])]
     .groupby("Flexibility")["Most frequent upstream"]
     .value_counts(normalize=False)
     # .mul(100)
     .rename('Percentage (by clade)')
     .reset_index().pipe((seaborn.lineplot, 'data'), x="Most frequent upstream", y='Percentage (by clade)', ax=ax2))
    # (df_plot[df_plot["Flexibility"] == "DC(x,0)"].pipe((seaborn.barplot, 'data'), x="Most frequent upstream", y='DC(x,f)', ax=ax2, order=range(*xlim)))
    # seaborn.distplot(df_plot[df_plot["Flexibility"] == "DC(x,0)"]["Most frequent upstream"], ax=ax2, norm_hist=False, bins=50)
    ax2.set_xlabel("")
    ax2.set_ylabel("Number of\ncomponents")
    ax2.set_xticklabels([])
    ax2.set_xlim(*xlim)
    ax2.set_ylim(0, None)
    ax2.yaxis.set_major_formatter(FuncFormatter(number_formatter))
    fig3.align_ylabels()
    plt.tight_layout()
    plt.savefig(next_name(env["pd-work"]))
    plt.show()

    xlim = [-10, 100]
    fig, ax = plt.subplots()
    # gs = fig3.add_gridspec(4, 1)
    (df_plot.pipe((seaborn.lineplot, 'data'), x="Most frequent upstream", y='DC(x,f)', hue="Flexibility", ax=ax))
    ax.set_ylim(0, 1)
    ax.set_xlim(*xlim)
    ax.set_xlabel("Most frequent distance to upstream gene")

    plt.tight_layout()
    plt.savefig(next_name(env["pd-work"]))
    plt.show()

    # (df[(df["Most frequent upstream"] < 10) & (df["Most frequent upstream"] > -10)]
    #  .groupby("Ancestor")["Most frequent upstream"]
    #  .value_counts()
    #  .rename('number')
    #  .reset_index()
    #  .pipe((seaborn.catplot, 'data'), x="Most frequent upstream", y='number', hue="Ancestor",
    #        kind='point', scale=0.5, legend=False,
    #        palette=CM.get_map("ancestor"), aspect=1.5
    #        ))
    #
    # plt.legend(loc="best", title="Clade")
    # figure_options = FigureOptions(
    #     save_fig=next_name(pd_work),
    #     xlabel="Most frequent distance to upstream gene",
    #     ylabel="Number of components"
    # )
    # plt.xlabel(figure_options.xlabel)
    # plt.ylabel(figure_options.ylabel)
    # save_figure(figure_options)
    #
    # plt.show()

    # f, ax1 = plt.subplots()
    # ax2 = ax1.twinx()
    # for ancestor, df_group in df.groupby("Ancestor"):
    #     seaborn.distplot(df_group["Most frequent upstream"], kde=False, ax=ax1)
    #
    #     # ax2.set_ylim(0, 3)
    #     ax2.yaxis.set_ticks([])
    #     seaborn.kdeplot(df_group["Most frequent upstream"], ax=ax2)
    #     ax1.set_xlabel('x var')
    #     ax1.set_ylabel('Counts')
    # # g = seaborn.FacetGrid(df, hue="Ancestor")
    # # g = g.map(seaborn.distplot, "Most frequent upstream", hist=True)
    # plt.show()

    # print(df["Most frequent upstream"].value_counts(normalize=True))
    #
    # sns.lmplot(
    #     df, "Most frequent upstream", "PC(x,0)", hue="Ancestor", sns_kwargs={
    #         "scatter": False, "lowess": True,
    #         "palette": CM.get_map("ancestor")
    #     },
    #     figure_options=FigureOptions(save_fig=next_name(pd_work), xlim=[-7, None], ylim=[0,1]),
    # )

    # sns.lmplot(
    #     df, "Most frequent upstream", "PC(x,3)", hue="Ancestor", sns_kwargs={
    #         "scatter": False, "lowess": True,
    #         "palette": CM.get_map("ancestor")
    #     },
    #     figure_options=FigureOptions(save_fig=next_name(pd_work), xlim=[-7, None], ylim=[0,1])
    # )

    # NCBI sensitivity
    # collect:
    # average 5' per ancestor, r,

    # ranges = [(-5,0), (0, 10), (10, 30), (30, 50), (50, 70)]
    # list_collect = list()
    # for r in ranges:
    #
    #     r_filter = (df["Most frequent upstream"] >= r[0]) & (df["Most frequent upstream"] < r[1])
    #
    #     df_summary_per_gcfid = get_summary_per_gcfid(df[r_filter])
    #     # viz_summary_per_gcfid(env, df_summary_per_gcfid, title=str(r))
    #
    #
    #     df_summary_per_gcfid = df_summary_per_gcfid.groupby("Ancestor", as_index=False).mean()
    #     df_summary_per_gcfid["Range"] = str(r)
    #     list_collect.append(df_summary_per_gcfid)
    #
    #
    # df_tmp = pd.concat(list_collect, sort=False)
    #
    # sns.catplot(df_tmp, "Range", "(GMS2=SBSP)!=NCBI % GMS2=SBSP", hue="Ancestor", kind="point",
    #             sns_kwargs={"palette": CM.get_map("ancestor")})
    #
    # sns.catplot(df_tmp, "Range", "GMS2=SBSP", hue="Ancestor", kind="point",
    #             sns_kwargs={"palette": CM.get_map("ancestor")})
    #
    # # do not average per gcfid - average per ancestor
    # list_collect = list()
    #
    # range_avgs = list()
    # range_label = list()
    #
    # for r in ranges:
    #     r_filter = (df["Most frequent upstream"] >= r[0]) & (df["Most frequent upstream"] < r[1])
    #     df_r = df[r_filter]
    #
    #     for ancestor, df_group in df_r.groupby("Ancestor", as_index=False):       # type: str, pd.DataFrame
    #
    #         f_gms2_eq_sbsp_with_ncbi_pred = (df_group["GMS2=SBSP"]) & (df_group["NCBI"])
    #         f_gms2_eq_sbsp_not_eq_ncbi = (f_gms2_eq_sbsp_with_ncbi_pred) & (df_group["(GMS2=SBSP)!=NCBI"])
    #
    #         sensitivity = 100 * f_gms2_eq_sbsp_not_eq_ncbi.sum() / float(f_gms2_eq_sbsp_with_ncbi_pred.sum())
    #         list_collect.append({
    #             "Ancestor": ancestor,
    #             "Range": str(r),
    #             "range_avg": (r[1] + r[0]) / 2.0,
    #             "(GMS2=SBSP)!=NCBI % GMS2=SBSP": sensitivity,
    #             "GMS2=SBSP": f_gms2_eq_sbsp_with_ncbi_pred.sum()
    #         })
    #
    #     range_label.append(r)
    #     range_avgs.append((r[1]+r[0]) / 2.0)
    #
    # df_tmp = pd.DataFrame(list_collect)
    #
    # sns.catplot(df_tmp, "Range", "(GMS2=SBSP)!=NCBI % GMS2=SBSP", hue="Ancestor", kind="point",
    #             sns_kwargs={"palette": CM.get_map("ancestor")})
    #
    # sns.catplot(df_tmp, "Range", "GMS2=SBSP", hue="Ancestor", kind="point",
    #             sns_kwargs={"palette": CM.get_map("ancestor")})
    #
    # ancestors = list(set(df_tmp["Ancestor"]))
    # fig, axes = plt.subplots(len(ancestors), 1, sharex="all",)
    # for ancestor, ax in zip(ancestors, axes.ravel()):       # type: str, plt.Axes
    #     ax2 = ax.twinx()
    #     curr_df = df_tmp[df_tmp["Ancestor"] == ancestor]
    #     seaborn.lineplot("range_avg", "(GMS2=SBSP)!=NCBI % GMS2=SBSP", data=curr_df, ax=ax)
    #     seaborn.lineplot("range_avg", "GMS2=SBSP", data=curr_df,
    #                  color='r', legend=False,
    #                  ax=ax2)
    #     ax.set_ylabel(None)
    #     ax2.set_ylabel(None)
    #     ax.set_xlabel("Range Average")
    #
    #
    # plt.xticks(range_avgs, range_label)
    # plt.show()
    #
    # fig, ax = plt.subplots()
    # ax2 = ax.twinx()
    # seaborn.lineplot("range_avg", "(GMS2=SBSP)!=NCBI % GMS2=SBSP", data=df_tmp, ax=ax, color="b", ci=None, hue="Ancestor")
    # seaborn.lineplot("range_avg", "GMS2=SBSP", data=df_tmp, ci=None,
    #                  color='r', legend=False,
    #                  ax=ax2, hue="Ancestor")
    # # plt.xticks(range_avgs, range_label)
    # ax.set_ylim([0, None])
    # ax2.set_ylim([0, None])
    #
    # ax.set_ylabel("NCBI 5' error rate vs GMS2=SBSP")
    # ax2.set_ylabel("Number of GMS2=SBSP genes")
    # ax.set_xlabel("Range Average")
    #
    # ax.yaxis.label.set_color('b')
    # ax2.yaxis.label.set_color('r')
    # ax.set_xlabel("Distance to upstream gene (nt)")
    # plt.show()


def analyze_support(env, df):
    # type: (Environment, pd.DataFrame) -> None

    sns.lineplot(df, "Genome GC", "Support", hue="Ancestor")
    sns.catplot(df, "Ancestor", "Support",
                sns_kwargs={"palette": CM.get_map("ancestor")})

    df_gcfid = df.groupby(["GCFID", "Ancestor"], as_index=False).mean()
    sns.catplot(
        df_gcfid, "Ancestor", "Support",
        sns_kwargs={"palette": CM.get_map("ancestor")},
        figure_options=FigureOptions(
            xlabel="Clade",
            ylabel="Average number of targets\nper query, per species",
            ylim=[0,None],
            save_fig=next_name(env['pd-work'])
        )
    )
    pass


def viz_summary_per_gcfid_per_step(env, df):
    # type: (Environment, pd.DataFrame) -> None
    pd_work = env['pd-work']

    list_df = list()
    for step in ["A", "B", "C"]:
        df_summary_per_gcfid = get_summary_per_gcfid(df[df["Predicted-at-step"] == step])
        df_summary_per_gcfid["SBSP Step"] = step
        list_df.append(df_summary_per_gcfid)

    df_per_gcfid_per_step = pd.concat(list_df, sort=False)

    sns.catplot(
        df_per_gcfid_per_step, "Ancestor", "(GMS2=SBSP)!=NCBI % GMS2=SBSP", hue="SBSP Step", kind="box",
        legend_loc="best",
        figure_options=FigureOptions(
            save_fig=next_name(pd_work),
            xlabel="Clade",
            ylabel=f"Diff(PGAP,{TOOLp})",
            # label_fontsize="x-large",
            # title_fontsize="x-large",

        ),

    )


def viz_analysis_per_query(env, df, **kwargs):
    # type: (Environment, pd.DataFrame, Dict[str, Any]) -> None

    pd_work = env["pd-work"]
    df.reset_index(inplace=True)

    df["(GMS2=SBSP)!=NCBI"] = df["GMS2=SBSP"] & df["NCBI"] & ~df["GMS2=SBSP=NCBI"]
    df["(GMS2=SBSP)!=Prodigal"] = df["GMS2=SBSP"] & df["Prodigal"] & ~df["GMS2=SBSP=Prodigal"]

    print(df[["Ancestor", "SBSP", "GMS2=SBSP"]].groupby("Ancestor", as_index=False).sum().to_string())
    print(df[["Ancestor", "GCFID"]].groupby("Ancestor")["GCFID"].nunique().to_string())

    df.drop(df[((df["Ancestor"] == "FCB group") & (df["Genome GC"] > 0.6))].index, axis=0, inplace=True)
    # remove all step C
    # df.drop(df.index[df["Predicted-at-step"] == "C"], inplace=True)
    # df.loc[df["Predicted-at-step"] == "B", "Predicted-at-step"] = "C"
    # df.loc[df["Predicted-at-step"] == "U", "Predicted-at-step"] = "B"
    # df.drop(df.index[df["Support"] < 5], inplace=True)

    viz_summary_per_gcfid_per_step(env, df)


    df_summary_per_gcfid = get_summary_per_gcfid(df)
    viz_summary_per_gcfid(env, df_summary_per_gcfid)

    analyze_upstream_distances(env, df[~df["Upstream-distance"].isnull()].copy())
    analyze_kimura_distances(env, df[~df["Kimura-to-query"].isnull()].copy())
    analyze_support(env, df)


def read_analysis_per_query_to_df(pf_data):
    # type: (str) -> pd.DataFrame
    generic = lambda x: ast.literal_eval(x)
    conv = {
        # 'Kimura-to-query': generic,
        # 'Upstream-distance': generic
    }
    df = pd.read_csv(pf_data, converters=conv)

    return df


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None

    df = read_analysis_per_query_to_df(args.pf_data)
    df = df[df["Ancestor"] != "Alphaproteobacteria"].copy()
    # df = df.sample(100).copy()
    viz_analysis_per_query(env, df)


if __name__ == "__main__":
    main(my_env, parsed_args)
