import copy
import os
import logging
from typing import *
import numpy as np
import sbsp_general
from sbsp_general.dataframe import df_print_labels
from sbsp_general.general import except_if_not_in_set
from sbsp_general.labels import create_key_3prime_from_label
from sbsp_general.labels_comparison_detailed import LabelsComparisonDetailed
from sbsp_options.pipeline_sbsp import PipelineSBSPOptions
from sbsp_pbs_data.splitters import *
from sbsp_viz.labels_comparison_detailed import LabelsComparisonDetailedViz

log = logging.getLogger(__name__)


def df_add_is_true_start(df, true_labels, source, suffix_is_true_column="is-true",
                         suffix_3prime="3prime",
                         suffix_5prime_3prime="5prime-3prime",
                         coordinates_suffix=None):

    # type: (pd.DataFrame, sbsp_general.labels.Labels, str, str) -> None

    except_if_not_in_set(source, ["q-", "t"])

    key_3prime_to_gene_key = dict()
    for l in true_labels:
        key_3prime = create_key_3prime_from_label(l)
        key_3prime_to_gene_key[key_3prime] = sbsp_general.labels.create_gene_key_from_label(l)

    sbsp_general.general.df_add_5prime_3prime_key(df, source, suffix_5prime_3prime, coordinates_suffix=coordinates_suffix)


    column_is_true = "{}{}".format(source, suffix_is_true_column)
    column_5prime_3prime = "{}{}".format(source, suffix_5prime_3prime)

    df[column_is_true] = -1         # -1 means unknown

    for index, row in df.iterrows():

        # construct gene key
        row_5prime_3prime = row[column_5prime_3prime]

        [_, mv_accession, mv_left, mv_right, mv_strand] = row_5prime_3prime.strip().split(";")

        row_3prime_key = sbsp_general.general.create_3prime_key(None, mv_accession, mv_left, mv_right,
                                                                mv_strand)

        row_5prime_3prime_key = sbsp_general.general.create_gene_key(None, mv_accession, mv_left, mv_right,
                                                                mv_strand)

        if row_3prime_key in key_3prime_to_gene_key:

            # correct mv start
            if row_5prime_3prime_key == key_3prime_to_gene_key[row_3prime_key]:

                df.at[index, column_is_true] = 1
            # wrong mv start
            else:
                df.at[index, column_is_true] = 0


def df_add_distance_between_predicted_and_true(df, true_labels, source, suffix_distance_to_true="distance-to-true",
                         suffix_3prime="3prime",
                         suffix_5prime_3prime="5prime-3prime",
                         coordinates_suffix=None):

    # type: (pd.DataFrame, sbsp_general.labels.Labels, str, str) -> None

    except_if_not_in_set(source, ["q-", "t"])

    key_3prime_to_gene_key = dict()
    for l in true_labels:
        key_3prime = create_key_3prime_from_label(l)
        key_3prime_to_gene_key[key_3prime] = sbsp_general.labels.create_gene_key_from_label(l)

    sbsp_general.general.df_add_5prime_3prime_key(df, source, suffix_5prime_3prime, coordinates_suffix=coordinates_suffix)


    column_distance = "{}{}".format(source, suffix_distance_to_true)
    column_5prime_3prime = "{}{}".format(source, suffix_5prime_3prime)

    df[column_distance] = np.nan         # -1 means unknown

    for index, row in df.iterrows():

        # construct gene key
        row_5prime_3prime = row[column_5prime_3prime]

        [_, mv_accession, mv_left, mv_right, mv_strand] = row_5prime_3prime.strip().split(";")

        row_3prime_key = sbsp_general.general.create_3prime_key(None, mv_accession, mv_left, mv_right,
                                                                mv_strand)

        row_5prime_3prime_key = sbsp_general.general.create_gene_key(None, mv_accession, mv_left, mv_right,
                                                                mv_strand)

        if row_3prime_key in key_3prime_to_gene_key:

            if mv_strand == "+":

                [_,_,true_left,_,_] = key_3prime_to_gene_key[row_3prime_key].strip().split(";")
                predicted_left = mv_left

                df.at[index, column_distance] = float(predicted_left) - float(true_left)
            else:
                [_, _, _, true_right, _] = key_3prime_to_gene_key[row_3prime_key].strip().split(";")
                predicted_right = mv_right

                df.at[index, column_distance] = float(true_right) - float(predicted_right)







def pipeline_step_compute_accuracy(env, df, pipeline_options):
    # type: (Environment, pd.DataFrame, PipelineSBSPOptions) -> pd.DataFrame

    log.info("Pipeline Step: Compute accuracy...")

    from sbsp_io.labels import read_labels_from_file

    for genome in set(df["q-genome"]):
        pf_q_labels_true = os.path.join(env["pd-data"], genome, pipeline_options["fn-q-labels-true"])

        labels = read_labels_from_file(pf_q_labels_true, shift=0)

        df_add_is_true_start(df, labels, "q-", "is-true",
                                                    coordinates_suffix="-sbsp")
        df_add_distance_between_predicted_and_true(
            df, labels, "q-", "distance-to-true",
            coordinates_suffix="-sbsp")


    # get labels
    genome_to_pf_labels = df_print_labels(env, df, "q", suffix_coordinates="sbsp",
                                           suffix_fname="")

    # print accuracies
    from sbsp_general.labels_comparison import LabelsComparison
    genome_to_comparison = dict()

    for genome in genome_to_pf_labels:
        pf_q_labels_true = os.path.join(env["pd-data"], genome, pipeline_options["fn-q-labels-true"])

        genome_to_comparison[genome] = LabelsComparison(env, pf_q_labels_true, genome_to_pf_labels[genome])

        labels_a = read_labels_from_file(pf_q_labels_true)
        labels_b = read_labels_from_file(genome_to_pf_labels[genome])

        lcd = LabelsComparisonDetailed(labels_a, labels_b,
                                       name_a="Verified",
                                       name_b="SBSP",
                                       tag=genome,
                                       split_on_attributes=["predicted-at-step"])

        # LabelsComparisonDetailedViz(lcd).run(env["pd-work"])

    accuracy = LabelsComparison.stringify_genome_accuracies(genome_to_comparison, ",")
    import sbsp_io.general
    pf_accuracy = os.path.join(env["pd-work"], pipeline_options["fn-accuracy"])
    sbsp_io.general.write_string_to_file(accuracy, pf_accuracy)

    return df



def separate_msa_outputs_by_stats(env, df, dn_msa_output):
    # type: (Environment, pd.DataFrame, str) -> None

    if dn_msa_output is None:
        dn_msa_output = "msa_output"

    pd_msa_output_true = os.path.join(env['pd-work'], "{}_true".format(dn_msa_output))
    pd_msa_output_false = os.path.join(env['pd-work'], "{}_false".format(dn_msa_output))


    if not os.path.exists(pd_msa_output_true):
        os.makedirs(pd_msa_output_true)
    if not os.path.exists(pd_msa_output_false):
        os.makedirs(pd_msa_output_false)

    df_group = df.groupby(["q-is-true", "pf-msa-output"], as_index=False).agg("first")

    df_group_true = df_group[df_group["q-is-true"] == 1]
    df_group_false = df_group[df_group["q-is-true"] == 0]

    from shutil import copyfile

    for index, row in df_group_true.iterrows():
        fn_msa = os.path.basename(row["pf-msa-output"])
        copyfile(row["pf-msa-output"], os.path.join(pd_msa_output_true, fn_msa))
        try:
            copyfile("{}_nt".format(row["pf-msa-output"]), os.path.join(pd_msa_output_true, "{}_nt".format(fn_msa)))
        except IOError:
            pass

    for index, row in df_group_false.iterrows():

        fn_msa = os.path.basename(row["pf-msa-output"])
        copyfile(row["pf-msa-output"], os.path.join(pd_msa_output_false, fn_msa))

        try:
            copyfile("{}_nt".format(row["pf-msa-output"]), os.path.join(pd_msa_output_false, "{}_nt".format(fn_msa)))
        except IOError:
            pass
