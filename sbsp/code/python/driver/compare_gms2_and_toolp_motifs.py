# Karl Gemayel
# Georgia Institute of Technology
#
# Created: 3/17/20
import logging
import os
import argparse
from tqdm import tqdm
import pandas as pd
from typing import *
import re
import Bio.Seq
from Bio import SeqIO

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
from sbsp_general.shelf import next_name, compute_gc_from_file, relative_entropy
from sbsp_io.general import remove_p, mkdir_p
from sbsp_io.labels import read_labels_from_file
from sbsp_io.sequences import write_fasta_hash_to_file
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


def read_fasta_into_hash(fnsequences, stop_at_first_space=True):
    # type: (str) -> dict[str, Bio.Seq.Seq]
    # read all sequences
    sequences = {}  # will contain all protein sequences from file
    try:
        for record in SeqIO.parse(fnsequences, "fasta"):
            if stop_at_first_space:
                sequences[record.name.strip().split()[0]] = record.seq
            else:
                sequences[record.description.strip()] = record.seq

    except Exception:
        pass

    return sequences

def create_attribute_dict(attribute_string, delimiter=";"):
    # type: (str, str) -> Dict[str, Any]

    attributes = dict()

    for current in attribute_string.strip().split(sep=delimiter):
        if len(current.strip().split(maxsplit=1)) == 2:
            k, v = current.strip().split(maxsplit=1)
            attributes[k] = v

    return attributes

def read_gff(pf_labels, shift=-1):
    # type: (str, int) -> Dict[str, List[Dict[str, Any]]]


    labels = dict() # type: Dict[str, List[Dict[str, Any]]]

    pattern = re.compile(r"([^\t]+)\t([^\t]+)\t(CDS)\t(\d+)\t(\d+)\t([^\t]+)\t([+-])\t([^\t]+)\t([^\t]+)")

    with open(pf_labels, "r") as f:

        for line in f:

            line = line.strip()

            m = pattern.match(line)
            if m:

                attributes = create_attribute_dict(m.group(9))
                seqname = m.group(1).strip()

                label = {
                    "left" : int(m.group(4)) + shift,
                    "right" : int(m.group(5)) + shift,
                    "strand" : m.group(7),
                    "seqname" : m.group(1),
                    "attributes": attributes
                }

                if seqname not in labels.keys():
                    labels[seqname] = list()

                labels[seqname].append(label)
    return labels

def convert_multi_fasta_into_single_fasta(sequences, labels, new_seqname):
    # type: (Dict[str, str], Dict[str, List[Dict[str, Any]]], str) -> [Dict[str, str], Dict[str, List[Dict[str, Any]]]]

    sorted_seqnames = sorted(sequences.keys())

    joined_seq = {new_seqname: ""}
    joined_labels = {new_seqname: list()}

    curr_offset = 0

    for seqname in sorted_seqnames:

        if seqname not in labels.keys():
            continue

        seqname_sequence = sequences[seqname]
        seqname_labels = labels[seqname]

        for l in seqname_labels:
            new_l = l.copy()
            new_l["left"] += curr_offset
            new_l["right"] += curr_offset
            new_l["seqname"] = new_seqname

            joined_labels[new_seqname].append(new_l)


        joined_seq[new_seqname] += seqname_sequence
        curr_offset += len(seqname_sequence)

    return [joined_seq, joined_labels]

def write_lst(labels, pf_output, shift=+1):
    # type: (Dict[str, List[Dict[str, Any]]], str, int) -> None

    # seqname, source, feature, left, right, score, strand, frame, attributes
    with open(pf_output, "w") as f:
        f.write("# GeneMark.hmm-2 LST format\n"
"# GeneMark.hmm-2 prokaryotic version: 1.14\n"
"# File with sequence: tmpseq.fna\n"
"# File with native parameters: itr_1.mod\n"
"# Native species name and build: gms2-training\n"
"# File with MetaGeneMark parameters: /storage4/karl/sbsp/biogem/sbsp/bin_external/gms2/mgm_11.mod\n"
"# translation table: 11\n"
"# output date start: Mon Jun  8 09:26:44 2020\n\n")
        for seqname, seqname_labels in labels.items():
            f.write(f"SequenceID: {seqname}\n")
            counter = 1

            for counter, l in enumerate(seqname_labels):

                out = str(counter)
                out += " " + str(l["strand"])
                out += " " + str(l["left"] + shift)
                out += " " + str(l["right"] + shift)
                out += " " + str(l["right"] - l["left"] + 1)
                out += " " "nativebac" + " AGGAGG 6 1"
                out += " " + l["attributes"]["gene_type"]


                f.write(out + "\n")

def convert_multi_fasta_to_single(env, pf_sequences, pf_labels):
    # pf_sequence = sys.argv[1]
    # pf_labels = sys.argv[2]
    # pd_work = sys.argv[3]

    org_seq = read_fasta_into_hash(pf_sequences)
    org_labels = read_gff(pf_labels, shift=0)

    new_seq, new_labels = convert_multi_fasta_into_single_fasta(org_seq, org_labels, "anydef")
    pd_work = env["pd-work"]

    pf_new_seq = os_join(pd_work, "sequence_joined")
    pf_new_labels = os_join(pd_work, "labels_joined_lst")
    import os
    # write_gff(new_labels, os.path.join(pd_work, "labels_joined"), shift=0)
    write_lst(new_labels, pf_new_labels, shift=0)
    write_fasta_hash_to_file(new_seq, pf_new_seq)

    return pf_new_seq, pf_new_labels


def get_identital_labels(pf_gms2, pf_sbsp, pf_toolp, **kwargs):

    pf_lst = get_value(kwargs, "pf_lst", None)
    if pf_lst is not None:
        run_shell_cmd()

    else:
        run_shell_cmd(
            f"compp -a {pf_gms2} -b {pf_sbsp} -I -q -n | grep -v \"#\" > {pf_toolp}"
        )


def train_gms2_model(env, pf_new_seq, pf_new_labels, **kwargs):
    group = get_value(kwargs, "group", "A", default_if_none=True)
    pf_mod = os_join(env["pd-work"], "a.mod")
    cmd = f"cd {env['pd-work']}; "
    cmd += f"/storage4/karl/sbsp/biogem/sbsp/bin_external/gms2/biogem gms2-training -s {pf_new_seq} -l {pf_new_labels} -m {pf_mod} --order-coding 5 --order-noncoding 2 --only-train-on-native 1 --genetic-code 11 --order-start-context 2 --fgio-dist-thr 25 --genome-group {group} --ga-upstr-len-rbs 20 --align right --ga-width-rbs 6"
    run_shell_cmd(
        cmd
    )
    mod = GMS2Mod.init_from_file(pf_mod)
    remove_p(pf_mod)

    return mod

def train_and_create_models(env, pf_labels, pf_sequences, **kwargs):
    # type: (Environment, str, str) -> GMS2Mod
    pf_new_seq, pf_new_labels = convert_multi_fasta_to_single(env, pf_sequences, pf_labels)

    mod = train_gms2_model(env, pf_new_seq, pf_new_labels, **kwargs)
    remove_p(pf_new_labels)
    remove_p(pf_new_seq)

    return mod


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

def append_to_file(a_str, pf_output):
    # type: (str, str) -> None

    with open(pf_output, "a") as f:
        f.write(a_str)


def add_toolp_rbs_to_gms2_model(env, pf_sequence, pf_toolp, pf_gms2_mod, pf_new_mod, **kwargs):
    # type: (Environment, str, str, str, str) -> None

    group = get_value(kwargs, "group", None)

    # run toolp and create model file
    rbs_toolp = train_and_create_models(
        env,
        pf_labels=pf_toolp,
        pf_sequences=pf_sequence,
        group=group
    ).items["RBS_MAT"]      # type: Dict[str, List[float]]

    cmd = ""

    # remove RBS_MAT from new model
    cmd += " awk '{if ($1 == \"$RBS_MAT\") NR += 4 ; else print }' " + "{} > {}".format(pf_gms2_mod, pf_new_mod)

    run_shell_cmd(cmd)

    # write toolp RBS_MAT to new model file
    rbs_as_str = "\n\n$RBS_MAT\n"
    for i in sorted(rbs_toolp.keys()):
        rbs_as_str += str(i) + " " + " ".join([str(x) for x in rbs_toolp[i]]) + "\n"
    rbs_as_str += "\n\n"
    append_to_file(
        rbs_as_str, pf_new_mod
    )

    return


def run_gms2_prediction_with_model(pf_sequence, pf_new_mod, pf_new_pred):
    # type: (str, str, str) -> None

    prog = "/storage4/karl/sbsp/biogem/sbsp/bin_external/gms2/gmhmmp2"
    mgm_mod = "/home/karl/gms2-install/gms2_linux_64/mgm_11.mod"
    cmd = f"{prog} -m {pf_new_mod} -M {mgm_mod} -s {pf_sequence} -o {pf_new_pred} --format gff"
    run_shell_cmd(cmd)


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
