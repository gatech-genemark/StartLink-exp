# Karl Gemayel
# Georgia Institute of Technology
#
# Created: 3/17/20
import logging
import argparse
import pandas as pd
import numpy as np
from typing import *

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import sbsp_argparse.parallelization
import sbsp_log  # runs init in sbsp_log and configures logger

# Custom imports
from sbsp_container.genome_list import GenomeInfoList, GenomeInfo
from sbsp_container.gms2_mod import GMS2Mod
from sbsp_general import Environment
import sbsp_viz.sns as sns

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from sbsp_general.GMS2Noncoding import GMS2Noncoding
from sbsp_general.MotifModel import MotifModel
from sbsp_general.general import os_join, get_value, run_shell_cmd
from sbsp_general.labels import Labels
from sbsp_general.labels_comparison_detailed import LabelsComparisonDetailed
from sbsp_general.shelf import relative_entropy, run_gms2_prediction_with_model, train_and_create_models, \
    add_toolp_rbs_to_gms2_model, next_name
from sbsp_io.general import remove_p, convert_multi_fasta_to_single, read_lst, write_lst, read_gff, mkdir_p
from sbsp_io.labels import write_labels_to_file, read_labels_from_file
from sbsp_options.parallelization import ParallelizationOptions
from sbsp_parallelization.pbs import PBS
from sbsp_pbs_data.mergers import merge_identity
from sbsp_pbs_data.splitters import split_genome_info_list
from sbsp_viz.general import FigureOptions

parser = argparse.ArgumentParser("Description of driver.")

parser.add_argument('--pf-genome-list', required=True, help="Genome List")
sbsp_argparse.parallelization.add_parallelization_options(parser)

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

def train_gms2_model(env, pf_new_seq, pf_new_labels, **kwargs):
    group = get_value(kwargs, "group", "A", default_if_none=True)
    clean = get_value(kwargs, "clean", True)
    pf_mod = get_value(kwargs, "pf_mod", os_join(env["pd-work"], "a.mod"), default_if_none=True)

    cmd = f"cd {env['pd-work']}; "
    cmd += f"/storage4/karl/sbsp/biogem/sbsp/bin_external/gms2/biogem gms2-training -s {pf_new_seq} -l {pf_new_labels} -m {pf_mod} --order-coding 5 --order-noncoding 2 --only-train-on-native 1 --genetic-code 11 --order-start-context 2 --fgio-dist-thr 25 --genome-group {group} --ga-upstr-len-rbs 20 --align right --ga-width-rbs 6"
    run_shell_cmd(
        cmd
    )
    mod = GMS2Mod.init_from_file(pf_mod)
    if not clean:
        remove_p(pf_mod)

    return mod

def train_with_fraction_of_genes(env, gi, percent):
    # type: (Environment, GenomeInfo, float) -> [str, str]
    pf_gms2 = os_join(env["pd-runs"], gi.name, "gms2", "gms2.gff")
    pf_sequence = os_join(env["pd-data"], gi.name, "sequence.fasta")


def randomly_select_labels(pf_labels, pf_labels_percent, percent):
    # type: (str, str, float) -> None
    labels = [l for l in read_labels_from_file(pf_labels)]



    # for seqname in labels.keys():
    #     total = len(labels[seqname])
    #     tmp = sorted(np.random.choice(labels[seqname], size=int(total * percent / float(100)), replace=False),
    #                  key=lambda x: x["left"])
    #     new_labels[seqname] = tmp


    new_labels = Labels(
        sorted(np.random.choice(labels, size=int(len(labels) * percent / float(100)), replace=False),
                            key=lambda l: l.left()),
    )

    write_labels_to_file(new_labels, pf_labels_percent)

def logo_rbs_from_gms2_mod_file(pd_figures, pf_mod, title=""):
    # type: (str, str, str) -> None

    mod = GMS2Mod.init_from_file(pf_mod)
    mm = MotifModel(mod.items["RBS_MAT"], mod.items["RBS_POS_DISTR"])
    non = GMS2Noncoding(mod.items["NON_MAT"])
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1,2)
    import logomaker as lm
    lm.Logo(
        lm.transform_matrix(mm.pwm_to_df(), from_type="probability", to_type="information", background=non.pwm_to_array(0)),
        ax=axes[0])
    axes[0].set_title(title)
    axes[0].set_ylim(0,2)

    df_spacer = pd.DataFrame({"Distance from start": range(len(mm._spacer)), "Probability": mm._spacer})
    sns.lineplot(df_spacer, "Distance from start", "Probability", ax=axes[1], figure_options=FigureOptions(
        ylim=[0,0.4]
    ))
    plt.tight_layout()
    plt.savefig(next_name(pd_figures))

    plt.show()




def relative_entropy_analysis_for_gi_for_percent(env, pf_sequence, pf_labels, pf_mod, pf_verified, group, percent, pd_figures):
    # type: (Environment, str, str, str, str, str, float, str) -> Dict[str, Any]


    # 1)  randomly select percent of labels
    pf_labels_percent = os_join(env["pd-work"], "labels_percent.lst")
    pf_mod_percent = os_join(env["pd-work"], "model_percent.mod")
    pf_labels_predict = os_join(env["pd-work"], "labels_predict.lst")


    randomly_select_labels(pf_labels, pf_labels_percent, percent)

    # train new model
    mod_percent = train_and_create_models(env, pf_sequences=pf_sequence, pf_labels=pf_labels_percent, group=group, clean=False, pf_mod=pf_mod_percent)

    # add RBSB to GMS2 model
    add_toolp_rbs_to_gms2_model(env, pf_sequence, pf_labels_percent, pf_mod, pf_mod_percent)

    logo_rbs_from_gms2_mod_file(pd_figures, pf_mod_percent, str(percent))

    # run prediction with new model
    run_gms2_prediction_with_model(pf_sequence, pf_mod_percent, pf_labels_predict)

    # compare predictions
    lcd = LabelsComparisonDetailed(read_labels_from_file(pf_labels_predict), read_labels_from_file(pf_verified))

    mm = MotifModel(mod_percent.items["RBS_MAT"], mod_percent.items["RBS_POS_DISTR"])
    non = GMS2Noncoding(mod_percent.items["NON_MAT"])
    return {
        "RE": relative_entropy(mm, non),
        "RE Motif": relative_entropy(mm, non, component="motif"),
        "RE Spacer": relative_entropy(mm, non, component="spacer"),
        "Error": 100 - 100 * len(lcd.match_3p_5p('a')) / len(lcd.match_3p('a'))
    }


def set_up_labels_and_sequence_for_genome(env, gi):
    # type: (Environment, GenomeInfo) -> Dict[str, Any]
    pf_gms2 = os_join(env["pd-runs"], gi.name, "gms2", "gms2.gff")
    pf_sequence = os_join(env["pd-data"], gi.name, "sequence.fasta")
    pf_mod = os_join(env["pd-runs"], gi.name, "gms2", "GMS2.mod")

    mod = GMS2Mod.init_from_file(pf_mod)
    group = mod.items["GENOME_TYPE"].split("-")[1].upper()

    return {
        "pf_labels": pf_gms2,
        "pf_sequence": pf_sequence,
        "pf_verified": os_join(env["pd-data"], gi.name, "verified.gff"),
        "pf_mod": pf_mod,
        "group": group
    }



def relative_entropy_analysis_for_gi(env, gi, prl_options):
    # type: (Environment, GenomeInfo, ParallelizationOptions) -> pd.DataFrame
    # Steps:

    list_entries = list()

    # set up labels (as lst) and sequence for genome
    setup_info = set_up_labels_and_sequence_for_genome(env, gi)

    if prl_options["use-pbs"]:
        pd_figures = os_join(prl_options["pbs-pd-head"], gi.name)
    else:
        pd_figures = os_join(env["pd-work"], gi.name)
    mkdir_p(pd_figures)


    for percent in range(10, 101, 5):
        for trial in range(10):
            info = relative_entropy_analysis_for_gi_for_percent(
                env, pf_sequence=setup_info["pf_sequence"],
                pf_labels=setup_info["pf_labels"],
                group=setup_info["group"],
                pf_mod=setup_info["pf_mod"],
                pf_verified=setup_info["pf_verified"],
                percent=percent,
                pd_figures=pd_figures
            )

            list_entries.append({
                "Genome": gi.name,
                "Percent": percent,
                "Trial": trial,
                **info
            })

    return pd.DataFrame(list_entries)

def fix_names(r):
    # type: (pd.Series) -> str
    return "{}. {}".format(
        r["Genome"][0], r["Genome"].split("_")[1]
    )

def relative_entropy_analysis(env, gil, prl_options):
    # type: (Environment, GenomeInfoList, ParallelizationOptions) -> None

    list_df = list()
    if not prl_options["use-pbs"]:
        for gi in gil:
            list_df.append(relative_entropy_analysis_for_gi(env, gi, prl_options))
    else:
        pbs = PBS(env, prl_options, splitter=split_genome_info_list, merger=merge_identity)
        list_df = pbs.run(
            data={"gil": gil},
            func=relative_entropy_analysis_for_gi,
            func_kwargs={"env": env, "prl_options": prl_options}
        )


    df = pd.concat(list_df, ignore_index=True, sort=False)
    df["Genome"] = df.apply(fix_names, axis=1)
    # TODO: plot what's needed
    df.to_csv(os_join(env["pd-work"], "summary.csv"), index=False)


    pd_figures = os_join(env["pd-work"], "summary_figures")
    mkdir_p(pd_figures)

    sns.scatterplot(df, "Percent", "Error", figure_options=FigureOptions(ylim=[0,20], save_fig=next_name(pd_figures)))
    sns.lineplot(df, "RE", "Error", hue="Genome", figure_options=FigureOptions(ylim=[0,20], save_fig=next_name(pd_figures)))
    sns.lineplot(df, "RE Motif", "Error",hue="Genome", figure_options=FigureOptions(ylim=[0,20], save_fig=next_name(pd_figures)))
    sns.lineplot(df, "RE Spacer", "Error", hue="Genome",figure_options=FigureOptions(ylim=[0,20], save_fig=next_name(pd_figures)))
    sns.scatterplot(df, "RE Motif", "RE Spacer", hue="Genome", identity=True, figure_options=FigureOptions(save_fig=next_name(pd_figures)))

    sns.lmplot(df, "Percent", "Error", hue="Genome", figure_options=FigureOptions(ylim=[0, 20], save_fig=next_name(pd_figures)))
    sns.lmplot(df, "RE", "Error", hue="Genome", figure_options=FigureOptions(ylim=[0, 20], save_fig=next_name(pd_figures)))
    sns.lmplot(df, "RE Motif", "Error", hue="Genome", figure_options=FigureOptions(ylim=[0, 20], save_fig=next_name(pd_figures)))
    sns.lmplot(df, "RE Spacer", "Error", hue="Genome", figure_options=FigureOptions(ylim=[0, 20], save_fig=next_name(pd_figures)))
    sns.lmplot(df, "Percent", "RE", hue="Genome",
               figure_options=FigureOptions(save_fig=next_name(pd_figures)))





def main(env, args):
    # type: (Environment, argparse.Namespace) -> None
    gil = GenomeInfoList.init_from_file(args.pf_genome_list)
    prl_options = ParallelizationOptions.init_from_dict(env, vars(args))



    relative_entropy_analysis(env, gil, prl_options)


if __name__ == "__main__":
    main(my_env, parsed_args)
