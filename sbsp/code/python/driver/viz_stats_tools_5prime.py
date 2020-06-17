# Karl Gemayel
# Georgia Institute of Technology
#
# Created: 3/17/20
import logging
import argparse
import math

import pandas as pd
from typing import *
import matplotlib.pyplot as plt

# noinspection All
from skmisc.loess._loess import loess

# noinspection PyUnresolvedReferences
import sbsp_log  # runs init in sbsp_log and configures logger

# Custom imports
from sbsp_general import Environment
from sbsp_general.shelf import next_name
# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from sbsp_viz.general import FigureOptions, save_figure
from sbsp_viz.shelf import loess_with_stde, create_mappable_for_colorbar

parser = argparse.ArgumentParser("Description of driver.")

parser.add_argument('--pf-stats', required=True)
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

def get_tags_for_5prime(df):
    # type: (pd.DataFrame) -> List[str]
    return [x.strip().split(":")[1] for x in df.columns if x.startswith("5p:")]


def compute_percentages(df):
    # type: (pd.DataFrame) -> None

    # get all tool combinations
    list_tags = get_tags_for_5prime(df)

    for tag in list_tags:
        df[f"M:{tag}"] = 100 - 100 * df[f"5p:{tag}"] / df[f"3p:{tag}"]



# def loess_with_stde(x, y, ax, label):
#     from skmisc.loess import loess
#
#     l = loess(x, y)
#     l.fit()
#     pred = l.predict(x, stderror=True)
#     conf = pred.confidence()
#     stderr = pred.stderr
#     lowess = pred.values
#     ll = conf.lower
#     ul = conf.upper
#
#     ul = lowess + 1.96 * stderr
#     ll = lowess - 1.96 * stderr
#
#     # ax.plot(x, y, '+')
#     ax.plot(x, lowess, label=label)
#     ax.fill_between(x, ll, ul, alpha=.33)
#




import numpy as np


def lowess(x, y, f=1./3.):
    """
    Basic LOWESS smoother with uncertainty.
    Note:
        - Not robust (so no iteration) and
             only normally distributed errors.
        - No higher order polynomials d=1
            so linear smoother.
    """
    # get some paras
    xwidth = f*(x.max()-x.min()) # effective width after reduction factor
    N = len(x) # number of obs
    # Don't assume the data is sorted
    order = np.argsort(x)
    # storage
    y_sm = np.zeros_like(y)
    y_stderr = np.zeros_like(y)
    # define the weigthing function -- clipping too!
    tricube = lambda d : np.clip((1- np.abs(d)**3)**3, 0, 1)
    # run the regression for each observation i
    for i in range(N):
        dist = np.abs((x[order][i]-x[order]))/xwidth
        w = tricube(dist)
        # form linear system with the weights
        A = np.stack([w, x[order]*w]).T
        b = w * y[order]
        ATA = A.T.dot(A)
        ATb = A.T.dot(b)
        # solve the syste
        sol = np.linalg.solve(ATA, ATb)
        # predict for the observation only
        yest = A[i].dot(sol)# equiv of A.dot(yest) just for k
        place = order[i]
        y_sm[place]=yest
        sigma2 = (np.sum((A.dot(sol) -y [order])**2)/N )
        # Calculate the standard error
        y_stderr[place] = np.sqrt(sigma2 *
                                A[i].dot(np.linalg.inv(ATA)
                                                    ).dot(A[i]))
    return y_sm, y_stderr

def plot_per_tool_by_genome_type(env, df):
    # type: (Environment, pd.DataFrame) -> None

    list_tags = get_tags_for_5prime(df)

    num_tags = len(list_tags)

    fig, ax = plt.subplots(2, math.ceil(num_tags/2), sharey="all", sharex="all")
    fig.add_axes([.91, .3, .03, .4])
    cbar_ax = fig.axes[-1]
    #
    # save_figure(FigureOptions(
    #     save_fig=next_name(env["pd-work"])
    #         ), fig)
    #
    # plt.show()
    # return

    import numpy as np
    kws = {
        # "levels": np.arange(0, 1, 0.2),
        # "vmin": 0, "vmax": 0.55,
        # "norm": True
        "xlim": [0.2, 0.8],
        "ylim": [0, 35],
        "cbar_max": 1,
        "num_steps": 35,
    }

    cbar_enable = {        "cbar_ax": cbar_ax, "cbar": True,}

    counter = 0
    for tag, c, a in zip(list_tags, ["b", "g", "r", "o"], ax.ravel()):
        x, y, y_l, y_u = loess_with_stde(df, "GC", f"M:{tag}", a, tag.replace("=", ","),
                                         **kws, **cbar_enable if counter == 0 else dict())

        a.set_title(tag.replace("=", ",").replace("NCBI", "PGAP").replace("GMS2", "GeneMarkS-2"))
        # a.set_ylim([65,100])
        # a.set_ylim([0, 35])
        # eps_x = [z for z in a.get_ylim()]
        # eps_x[0] -= 0.01
        # eps_x[1] += 0.01
        #
        # a.set_xlim(eps_x)
        # if counter % 2 == 0:
        #     a.set_ylabel("Percentage of gene-start differences")
        # if counter >= math.ceil(num_tags/2):
        #     a.set_xlabel("GC")
        counter += 1

        mappable = a.collections[0]


    # plt.legend(loc="best")
    figure_options = FigureOptions(
        save_fig=next_name(env["pd-work"])
            )
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axes
    plt.tick_params(top=False, bottom=False, left=False, right=False, which="both",
                    labelbottom=False, labeltop=False, labelleft=False, labelright=False)
    plt.xlabel("GC", labelpad=30)
    plt.ylabel("Percentage of gene-start differences", labelpad=30)
    # plt.xlabel("GC")
    # plt.ylabel("Percent 5' Match")

    # mappable=create_mappable_for_colorbar(np.arange(0, 0.4, 0.05), "Reds")
    # plt.colorbar(mappable, cax=cbar_ax, cmap="Reds")
    fig.tight_layout(rect=[-0.02, -0.02, .9, 1])

    # plt.tight_layout()
    # FigureOptions.set_properties_for_axis(ax, figure_options)

    save_figure(figure_options, fig)
    plt.show()
    #
    # for tag in list_tags:
    #     sns.jointplot(df, "GC", f"M:{tag}")
    #
    #
    # x = df["GC"].values
    # y = df[f"M:{list_tags[0]}"].values
    # order = np.argsort(x)
    # # run it
    # y_sm, y_std = lowess(x, y, f=1. / 5.)
    # # plot it
    # plt.plot(x[order], y_sm[order], color='tomato', label='LOWESS')
    # plt.fill_between(x[order], y_sm[order] - 1.96 * y_std[order],
    #                  y_sm[order] + 1.96 * y_std[order], alpha=0.3, label='LOWESS uncertainty')
    # # plt.plot(x, y, 'k.', label='Observations')
    # # plt.legend(loc='best')
    # # run it
    # y_sm, y_std = lowess(x, y, f=1. / 5.)
    # # plot it
    # plt.plot(x[order], y_sm[order], color='tomato', label='LOWESS')
    # plt.fill_between(x[order], y_sm[order] - y_std[order],
    #                  y_sm[order] + y_std[order], alpha=0.3, label='LOWESS uncertainty')
    # # plt.plot(x, y, 'k.', label='Observations')
    # plt.legend(loc='best')
    # plt.show()

    # calculate a 60 day rolling mean and plot
    # calculate a 60 day rolling mean and plot

    # df_stacked = stack_columns_as_rows(
    #     df, [f"M:{tag}" for tag in list_tags], "Percent 5p Match", [f"M:{tag}" for tag in list_tags], "Tools"
    # )
    #
    #
    # sns.lmplot(
    #     df_stacked, "GC", "Percent 5p Match", hue="Tools",
    #     figure_options=FigureOptions(
    #         xlabel="Genome GC",
    #         ylim=[70, 100]
    #     ),
    #     legend_loc="best",
    #     sns_kwargs={"scatter_kws": {"s": 5, "alpha": 0.3}, "lowess": False, "scatter": False, "aspect": 1.5}
    # )
    # # sns.tsplot(df_stacked, "GC", "Percent 5p Match", hue="Tools", sns_kwargs={"ci":"sd"})
    # fig, ax = plt.subplots(1, 1)
    # seaborn.lineplot(df["GC"], df[f"M:{list_tags[0]}"])
    # # seaborn.tsplot(df, "GC", f"M:{list_tags[0]}" , ci="sd")
    # plt.show()



    plt.show()
def main(env, args):
    # type: (Environment, argparse.Namespace) -> None
    df = pd.read_csv(args.pf_stats)
    # df = df.sample(50)
    compute_percentages(df)

    # cleanup
    df = df[(df["M:GMS2=NCBI"] < 50) & (df["M:GMS2=Prodigal"] < 50)].copy()     # type: pd.DataFrame
    # df = df.sample(1000)
    df.sort_values("GC", inplace=True)

    plot_per_tool_by_genome_type(env, df)

if __name__ == "__main__":
    main(my_env, parsed_args)
