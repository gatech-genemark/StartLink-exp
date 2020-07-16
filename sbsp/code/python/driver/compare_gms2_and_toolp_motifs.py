# Karl Gemayel
# Georgia Institute of Technology
#
# Created: 3/17/20
import logging
import argparse
from tqdm import tqdm
import pandas as pd

import matplotlib.pyplot as plt
import logomaker as lm


# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import sbsp_log  # runs init in sbsp_log and configures logger

# Custom imports
from sbsp_container.genome_list import GenomeInfoList, GenomeInfo
from sbsp_container.gms2_mod import GMS2Mod
from sbsp_general import Environment
from sbsp_general.GMS2Noncoding import GMS2Noncoding
from sbsp_general.MotifModel import MotifModel

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from sbsp_general.general import os_join, run_shell_cmd, get_value
from sbsp_general.labels_comparison_detailed import LabelsComparisonDetailed
from sbsp_general.shelf import next_name, compute_gc_from_file, relative_entropy, run_gms2_prediction_with_model, \
    train_and_create_models, add_toolp_rbs_to_gms2_model
from sbsp_io.general import mkdir_p
from sbsp_io.labels import read_labels_from_file
from sbsp_viz.general import FigureOptions

parser = argparse.ArgumentParser("Description of driver.")

parser.add_argument('--pf-genome-list', required=True)
parser.add_argument('--verified', required=False, action="store_true", default=False)

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


def get_identital_labels(pf_gms2, pf_sbsp, pf_toolp, **kwargs):

    pf_lst = get_value(kwargs, "pf_lst", None)
    if pf_lst is not None:
        run_shell_cmd()

    else:
        run_shell_cmd(
            f"compp -a {pf_gms2} -b {pf_sbsp} -I -q -n | grep -v \"#\" > {pf_toolp}"
        )


def compare_gms2_and_toolp_motifs_for_gi(env, gi):
    # type: (Environment, GenomeInfo) -> [GMS2Mod, GMS2Mod]

    pf_gms2 = os_join(env["pd-runs"], gi.name, "gms2", "gms2.gff")
    pf_sbsp = os_join(env["pd-runs"], gi.name, "sbsp_submission/accuracy", f"{gi.name}.gff")
    pf_sequence = os_join(env["pd-data"], gi.name, "sequence.fasta")
    pf_toolp = os_join(env["pd-work"], "toolp.gff")

    # get toolp predictions
    get_identital_labels(
        pf_gms2, pf_sbsp, pf_toolp
    )

    # get gms2 model
    mod_gms2 = train_and_create_models(
        env,
        pf_labels = pf_gms2,
        pf_sequences = pf_sequence
    )

    mod_toolp = train_and_create_models(
        env,
        pf_labels=pf_toolp,
        pf_sequences=pf_sequence
    )

    return mod_gms2, mod_toolp


def compare_gms2_start_predictions_with_motif_from_toolp(env, gi):
    # type: (Environment, GenomeInfo) -> float

    pf_gms2 = os_join(env["pd-runs"], gi.name, "gms2", "gms2.gff")
    pf_gms2_mod = os_join(env["pd-runs"], gi.name, "gms2", "GMS2.mod")
    pf_sbsp = os_join(env["pd-runs"], gi.name, "sbsp_submission/accuracy", f"{gi.name}.gff")
    pf_sequence = os_join(env["pd-data"], gi.name, "sequence.fasta")
    pf_toolp = os_join(env["pd-work"], "toolp.gff")

    # get toolp predictions
    get_identital_labels(
        pf_gms2, pf_sbsp, pf_toolp
    )

    # create new motif model with toolp and add it to new model file
    pf_new_mod = os_join(env["pd-work"], "toolp.mod")
    add_toolp_rbs_to_gms2_model(env, pf_sequence, pf_toolp, pf_gms2_mod, pf_new_mod)

    # run prediction with new model
    pf_new_pred = os_join(env["pd-work"], "new_pred.gff")
    run_gms2_prediction_with_model(pf_sequence, pf_new_mod, pf_new_pred)

    # compare predictions
    lcd = LabelsComparisonDetailed(read_labels_from_file(pf_gms2), read_labels_from_file(pf_new_pred))

    return 100 * len(lcd.match_3p_5p('a')) / len(lcd.match_3p('a'))



def compare_gms2_start_predictions_with_motif_from_toolp_verified(env, gi, **kwargs):
    # type: (Environment, GenomeInfo) -> [float, float]

    group = get_value(kwargs, "group", None)

    pf_gms2 = os_join(env["pd-runs"], gi.name, "gms2", "gms2.gff")
    pf_gms2_mod = os_join(env["pd-runs"], gi.name, "gms2", "GMS2.mod")
    pf_sbsp = os_join(env["pd-runs"], gi.name, "sbsp_submission/accuracy", f"{gi.name}.gff")
    pf_sequence = os_join(env["pd-data"], gi.name, "sequence.fasta")
    pf_toolp = os_join(env["pd-work"], "toolp.gff")
    pf_verified = os_join(env["pd-data"], gi.name, "verified.gff")

    # get toolp predictions
    get_identital_labels(
        pf_gms2, pf_sbsp, pf_toolp
    )

    # create new motif model with toolp and add it to new model file
    pf_new_mod = os_join(env["pd-work"], "toolp.mod")
    add_toolp_rbs_to_gms2_model(env, pf_sequence, pf_toolp, pf_gms2_mod, pf_new_mod, group=group)

    # run prediction with new model
    pf_new_pred = os_join(env["pd-work"], "new_pred.gff")
    run_gms2_prediction_with_model(pf_sequence, pf_new_mod, pf_new_pred)

    # compare predictions
    lcd1 = LabelsComparisonDetailed(read_labels_from_file(pf_gms2), read_labels_from_file(pf_verified))
    lcd2 = LabelsComparisonDetailed(read_labels_from_file(pf_new_pred), read_labels_from_file(pf_verified))

    return [100 * len(lcd.match_3p_5p('a')) / len(lcd.match_3p('a')) for lcd in [lcd1, lcd2]]

def fix_names(genome):
    # type: (str) -> None
    return "{}. {}".format(
        genome[0], genome.split("_")[1]
    )

def main(env, args):
    # type: (Environment, argparse.Namespace) -> None

    gil = GenomeInfoList.init_from_file(args.pf_genome_list)

    pd_figures = os_join(env["pd-work"], "figures")
    mkdir_p(pd_figures)


    list_run_info = list()

    for gi in tqdm(gil, total=len(gil)):
        # get gms2 and toolp models
        mod_gms2, mod_toolp = compare_gms2_and_toolp_motifs_for_gi(env, gi)

        group = mod_gms2.items["GENOME_TYPE"].split("-")[1].upper()


        mm_gms2 = MotifModel(mod_gms2.items["RBS_MAT"], None)
        mm_toolp = MotifModel(mod_toolp.items["RBS_MAT"], None)
        non_gms2 = GMS2Noncoding(mod_gms2.items["NON_MAT"])

        df_gms2 = mm_gms2.pwm_to_df()
        df_toolp = mm_toolp.pwm_to_df()

        fig, axes = plt.subplots(1, 2, sharex="all", sharey="all", figsize=(8, 4))

        # relative
        rel_mat = lm.transform_matrix(df_gms2, from_type="probability", to_type="information")
        lm.Logo(rel_mat, color_scheme="classic", ax=axes[0])
        axes[0].set_ylim(*[0, 2])
        axes[0].set_title("GeneMarkS-2")

        # shannon
        sha_mat = lm.transform_matrix(df_toolp, from_type="probability", to_type="information")
        lm.Logo(sha_mat, color_scheme="classic", ax=axes[1])
        axes[1].set_ylim(*[0, 2])
        axes[1].set_title("StartLink+")
        plt.tight_layout()
        plt.savefig(next_name(pd_figures))
        plt.show()

        rel_gms2 = relative_entropy(mm_gms2, non_gms2)
        rel_toolp = relative_entropy(mm_toolp, non_gms2)
        gc = 100 * compute_gc_from_file(os_join(env["pd-data"], gi.name, "sequence.fasta"))

        if not args.verified:
            list_run_info.append({
                "GC": gc,
                "Accuracy": 100 - compare_gms2_start_predictions_with_motif_from_toolp(env, gi),
                "RE GMS2": rel_gms2,
                "RE toolp": rel_toolp
            })
        else:
            # verified
            comp = compare_gms2_start_predictions_with_motif_from_toolp_verified(env, gi, group=group)
            list_run_info.append({
                "Genome": fix_names(gi.name),
                "Error": 100 - comp[0],
                "Tool": "GMS2",
                "RE": rel_gms2,
                "GC": gc
                })
            list_run_info.append({
                "Genome": fix_names(gi.name),
                "Error": 100 - comp[1],
                "Tool": "GMS2 with SL",
                "RE": rel_toolp,
                "GC": gc
                })

            print(list_run_info[-2:])

    import sbsp_viz.sns as sns
    if args.verified:
        df = pd.DataFrame(list_run_info)
        df.to_csv(next_name(env["pd-work"], ext="csv"))

        sns.lineplot(df, "Genome", "Error", hue="Tool", figure_options=FigureOptions(
            save_fig=next_name(env["pd-work"]),
            xlabel="Genome",
            ylabel="Error"))

        sns.lineplot(df, "Genome", "RE", hue="Tool",
                        figure_options=FigureOptions(
                            save_fig=next_name(env["pd-work"]),
                            xlabel="Genome",
                            ylabel="Relative entropy",
                        ))


    else:

        df = pd.DataFrame(list_run_info)
        sns.scatterplot(df, "GC", "Accuracy",
                    figure_options=FigureOptions(
                        save_fig=next_name(env["pd-work"]),
                        xlabel="GC",
                        ylabel="Percentage of different 5' ends",
                        ylim=[0,10],
                    ))

        df.to_csv(next_name(env["pd-work"], ext="csv"))

        sns.scatterplot(df, "RE GMS2", "RE toolp", identity=True, figure_options=FigureOptions(
            save_fig=next_name(env["pd-work"])
        ))


        print("Average Error: {}".format(df["Accuracy"].mean()))

        df = pd.DataFrame(list_run_info)
        df = df[df["Accuracy"] < 2].copy()
        sns.scatterplot(df, "GC", "Accuracy",
                    figure_options=FigureOptions(
                        save_fig=next_name(env["pd-work"]),
                        xlabel="GC",
                        ylabel="Percentage of different 5' ends",
                        ylim=[0,10],
                    ))

        sns.scatterplot(df, "RE GMS2", "RE toolp", identity=True, figure_options=FigureOptions(
            save_fig=next_name(env["pd-work"])
        ))

        print("Average Error: {}".format(df["Accuracy"].mean()))

        df.to_csv(next_name(env["pd-work"], ext="csv"))





if __name__ == "__main__":
    main(my_env, parsed_args)
