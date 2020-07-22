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
from sbsp_general.shelf import next_name
from sbsp_viz.colormap import ColorMap as CM
import sbsp_viz.sns as sns

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from sbsp_viz.general import FigureOptions

parser = argparse.ArgumentParser("Description of driver.")

parser.add_argument('--pf-data', required=True)

parser.add_argument('--with-mgm', required=False,  default=False, action="store_true")

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
    df["chunk-size"] /= 1000

    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()



    sns.lineplot(df[df["Tool"] == "SBSP"], "chunk-size", "percentage-common-3prime-and-5prime-from-common-3prime",
                 hue="Genome",
                 sns_kwargs={"palette": CM.get_map("verified"), "linestyle": "dashed"},
                 ax=ax,
                 legend=False,
                 figure_options=FigureOptions(
                     xlabel="Chunk size (mb)",
                     ylabel="Accuracy",
                     ylim=[74, 101],
                     save_fig=next_name(env["pd-work"])
                 ))

    for l in ax.lines:
        l.set_linestyle("--")

    sns.lineplot(df[df["Tool"] == "GMS2"], "chunk-size", "percentage-common-3prime-and-5prime-from-common-3prime",
                 hue="Genome",
                 sns_kwargs={"palette": CM.get_map("verified")},
                 legend_loc="best",
                 legend_ncol=2,
                 ax=ax)




    if args.with_mgm:
        y_max = ax.get_ylim()[1]
        ax.axvline(50, 0, y_max, color="grey", linestyle="dashed")
        ax.axhline(74, 5, 49, color="grey", linestyle="dashed")
        ax.annotate("MGM", (5, 72))

    if "MGM" in set(df["Tool"]):
        sns.lineplot(df[df["Tool"] == "MGM"], "chunk-size", "percentage-common-3prime-and-5prime-from-common-3prime",
                     hue="Genome",
                     sns_kwargs={"palette": CM.get_map("verified"), "linestyle": "-."},
                     ax=ax,
                     legend=False)

    for l in ax.lines[len(ax.lines)-5:]:
        l.set_linestyle(":")

    fo = FigureOptions(
                     xlabel="Chunk size (mb)",
                     ylabel="Accuracy",
                     ylim=[74,101],
                     save_fig=next_name(env["pd-work"])
                 )
    FigureOptions.set_properties_for_axis(ax, fo)
    plt.savefig(fo.save_fig)
    plt.show()


if __name__ == "__main__":
    main(my_env, parsed_args)
