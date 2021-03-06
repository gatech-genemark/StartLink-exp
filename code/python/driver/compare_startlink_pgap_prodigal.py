# Karl Gemayel
# Georgia Institute of Technology
#
# Created: 08/30/2020
import copy
import logging
import argparse
import os

import pandas as pd
from typing import *

# noinspection All
import seaborn

import pathmagic

# noinspection PyUnresolvedReferences
import sbsp_log  # runs init in sbsp_log and configures logger

# Custom imports
from sbsp_container.genome_list import GenomeInfoList, GenomeInfo
from sbsp_general import Environment
from sbsp_general.general import get_value, os_join, run_shell_cmd
from sbsp_general.labels import Labels
from sbsp_general.labels_comparison_detailed import LabelsComparisonDetailed
from sbsp_general.shelf import compute_gc_from_file, next_name
from sbsp_io.general import mkdir_p
from sbsp_io.labels import read_labels_from_file
from sbsp_io.objects import save_obj, load_obj
from sbsp_options.parallelization import ParallelizationOptions
from sbsp_parallelization.pbs import PBS
from sbsp_pbs_data.mergers import merge_identity
from sbsp_pbs_data.splitters import split_genome_info_list
import sbsp_argparse.parallelization

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from sbsp_viz.colormap import ColorMap
from sbsp_viz.general import FigureOptions
from sbsp_viz.sns import lmplot

parser = argparse.ArgumentParser("Run external prediction tools on genome list.")

parser.add_argument('--pf-genome-list', required=True, help="List of genomes")
parser.add_argument('--pf-checkpoint')
parser.add_argument('--pf-parallelization-options')

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


def get_differences_for_gi(env, gi, **kwargs):
    # type: (Environment, GenomeInfo, Dict[str, Any]) -> pd.DataFrame

    pf_sequences = os_join(env["pd-data"], gi.name, "sequence.fasta")
    gc = 100 * compute_gc_from_file(pf_sequences)

    clade = gi.attributes.get("ancestor")

    try:
        pf_startlink = os_join(env["pd-runs"], gi.name, "sbsp", "sbsp.gff")
        pf_ncbi = os_join(env["pd-data"], gi.name, "ncbi.gff")
        pf_gms2 = os_join(env["pd-runs"], gi.name, "gms2", "gms2.gff")
        pf_prodigal = os_join(env["pd-runs"], gi.name, "prodigal", "prodigal.gff")

        labels_startlink = read_labels_from_file(pf_startlink)
        labels_ncbi = read_labels_from_file(pf_ncbi)
        labels_gms2 = read_labels_from_file(pf_gms2)
        labels_prodigal = read_labels_from_file(pf_prodigal)

        labels_startlink_plus = LabelsComparisonDetailed(labels_startlink, labels_gms2).match_3p_5p("a")
    except FileNotFoundError:
        return pd.DataFrame()

    def _helper_stats(ref, pred, ref_name, ref_pred):
        # type: (Labels, Labels, str, str) -> Dict[str, Any]
        lcd = LabelsComparisonDetailed(ref, pred)

        return {
            "Difference": len(lcd.match_3p_not_5p("a")) / float(len(lcd.match_3p("a"))),
            "Tool": f"{ref_name},{ref_pred}",
            "GC": gc,
            "Clade": clade,
            "GCFID": gi.name,
            "Name": gi.attributes.get("name")
        }

    list_entries = [
        _helper_stats(labels_startlink_plus, labels_ncbi, "StartLink+", "PGAP"),
        _helper_stats(labels_startlink_plus, labels_prodigal, "StartLink+", "Prodigal"),
    ]

    return pd.DataFrame(list_entries)


def get_differences(env, gil, **kwargs):
    # type: (Environment, GenomeInfoList, Dict[str, Any]) -> pd.DataFrame()

    list_df = list()
    for gi in gil:
        list_df.append(
            get_differences_for_gi(env, gi, **kwargs)
        )

    return pd.concat(list_df, ignore_index=True, sort=False)


def plot_difference_to_startlink(env, df):
    # type: (Environment, pd.DataFrame) -> None

    df_combined = copy.deepcopy(df)
    df_combined["Clade"] = "Total"

    df = pd.concat([df, df_combined], ignore_index=True, sort=False)
    import matplotlib.pyplot as plt

    cmap = ColorMap.get_map("ancestor")
    cmap["Total"] = "black"
    df["Difference"] *= 100

    fig, axes = plt.subplots(1, 2, sharey="all", sharex="all", figsize=(12,4))

    for i, (tool, ax) in enumerate(zip(df["Tool"].unique(), axes)):

        df_tool = df[df["Tool"] == tool]
        for hue in df_tool["Clade"].unique():
            if hue == "Total":
                continue

            seaborn.regplot(
            "GC", "Difference", data=df_tool[df_tool["Clade"] == hue], lowess=True,
            scatter_kws={"s": 2, "alpha": 0.3}, color=cmap[hue],
            ax=ax, label=hue
            )

        seaborn.regplot(
            "GC", "Difference", data=df_tool[df_tool["Clade"] == "Total"], lowess=True,
            scatter_kws={"s": 0}, line_kws={"linestyle": "dashed"}, color=cmap["Total"],
            ax=ax, label="Average"
        )

        handles, labels = ax.get_legend_handles_labels()

        if i == 1:
            ax.set_ylabel("")
        else:
            ax.set_ylabel("Percentage of gene-start differences")
        ax.set_ylim(0, 20)
        ax.set_title(tool)

        print(df_tool.sort_values("Difference", ascending=False).to_csv())

    fig.subplots_adjust(0, 0.1, 1, 1)
    leg = fig.legend(handles=handles, labels=labels, loc="lower center", ncol=5, frameon=False)#, bbox_to_anchor=(0.5, -0.25))
    for lh in leg.legendHandles:
        lh.set_alpha(1)
        lh._sizes = [50]

    fig.tight_layout(rect=[0,0.1,1,1])
    plt.savefig(next_name(env["pd-work"]))
    plt.show()

    # lmplot(df, "GC", "Difference", hue="Clade",
    #        sns_kwargs={"col": "Tool", "lowess": True},
    #        figure_options=FigureOptions(
    #     save_fig=next_name(env["pd-work"]),
    #     xlabel="GC",
    #     ylabel="% difference in gene starts"
    # ))


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None

    gil = GenomeInfoList.init_from_file(args.pf_genome_list)

    pf_checkpoint = args.pf_checkpoint

    if not pf_checkpoint or not os.path.isfile(pf_checkpoint):
        prl_options = ParallelizationOptions.init_from_dict(env, vars(args))

        if prl_options is not None:
            pbs = PBS(env, prl_options,
                      splitter=split_genome_info_list,
                      merger=merge_identity
                      )

            output = pbs.run(
                data={"gil": gil},
                func=get_differences,
                func_kwargs={
                    "env": env,
                }
            )
            df = pd.concat(output, ignore_index=True, sort=False)
        else:
            df = get_differences(env, gil)

        if pf_checkpoint:
            save_obj(df, pf_checkpoint)
    else:
        df = load_obj(pf_checkpoint)

    plot_difference_to_startlink(env, df)


if __name__ == "__main__":
    main(my_env, parsed_args)
