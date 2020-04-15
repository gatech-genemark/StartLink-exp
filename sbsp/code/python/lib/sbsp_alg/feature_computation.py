import os
import logging
import numpy as np
import pandas as pd
from typing import *

from Bio.SubsMat import MatrixInfo as matlist

import sbsp_general
import sbsp_general.data
import sbsp_io.sequences
from sbsp_alg.phylogeny import global_alignment_aa_with_gap, k2p_distance
from sbsp_general.general import get_value, except_if_not_in_set
from sbsp_io.general import mkdir_p
from sbsp_general import Environment
from sbsp_alg.gene_distances import *


log = logging.getLogger(__name__)



def compute_distance(distance_type, q_align_aa, t_align_aa, q_align_nt, t_align_nt, on_fail=100.0, **kwargs):
    # type: (str, str, str, str, str, float, Dict[str, Any]) -> float

    try:
        if distance_type == "kimura":
            return k2p_distance(q_align_nt, t_align_nt, kimura_on_3rd=False)
        if distance_type == "kimura-on-3rd" or distance_type == "kimura3":
            return k2p_distance(q_align_nt, t_align_nt, kimura_on_3rd=True)

        if distance_type == "ds":
            return compute_distance_ds(q_align_nt, t_align_nt, **kwargs)

        if distance_type == "dn":
            return compute_distance_dn(q_align_nt, t_align_nt, **kwargs)

        if distance_type == "mismatch-aa":
            return compute_distance_mismatch_aa(q_align_aa, t_align_aa, **kwargs)

        if distance_type == "syn-fraction":
            return compute_synonymous_fraction(q_align_aa, t_align_aa, q_align_nt, t_align_nt)

        if distance_type == "non-syn-fraction":
            return compute_non_synonymous_fraction(q_align_aa, t_align_aa, q_align_nt, t_align_nt)

        if distance_type == "syn-poisson":
            return compute_synonymous_poisson(q_align_aa, t_align_aa, q_align_nt, t_align_nt)

        if distance_type == "non-syn-poisson":
            return compute_non_synonymous_poisson(q_align_aa, t_align_aa, q_align_nt, t_align_nt)

    except ValueError:
        return on_fail

    raise ValueError("Unknown distance type: {}".format(distance_type))



def count_aa_mismatches(seq_a, seq_b):
    if len(seq_a) != len(seq_b):
        raise ValueError("Sequence sizes are not the same: {} != {}".format(len(seq_a), len(seq_b)))

    from math import log, sqrt

    alignment_length = len(seq_a)

    matches = 0

    ungapped_length = 0

    for i in range(alignment_length):

        # ignore gaps
        if seq_a[i] == "-" or seq_b[i] == "-":
            continue

        ungapped_length += 1

        if seq_a[i] == seq_b[i]:
            matches += 1

    if ungapped_length == 0:
        return 0.0

    return matches / float(ungapped_length)


# FIXME: found in another file - put somewhere
def add_gaps_to_nt_based_on_aa(seq_nt, seq_aa_with_gaps):
    # type: (str, str) -> str

    # make sure number of nt is 3 times number of amino acids (exclude gaps)
    num_aa = len(seq_aa_with_gaps) - seq_aa_with_gaps.count('-')
    if len(seq_nt) != 3 * num_aa:
        raise ValueError("Number of nucleotides ({}) should be 3 times the number of amino acids ({})".format(
            len(seq_nt), num_aa))

    seq_nt_with_gaps = ""

    pos_in_nt = 0
    pos_in_aa_with_gaps = 0

    while pos_in_aa_with_gaps < len(seq_aa_with_gaps):

        curr_aa = seq_aa_with_gaps[pos_in_aa_with_gaps]

        # if gap
        if curr_aa == "-":
            seq_nt_with_gaps += "---"       # 3 nt gaps = 1 aa gap
        else:
            seq_nt_with_gaps += seq_nt[pos_in_nt:pos_in_nt+3]       # add next 3 nucleotides
            pos_in_nt += 3

        pos_in_aa_with_gaps += 1

    return seq_nt_with_gaps

def df_add_labeled_sequences(env, df, **kwargs):
    # type: (Dict, pd.DataFrame, Dict[str, Any]) -> pd.DataFrame

    """
    Add labeled sequences to data frame

    :param df:
    :param env: environment variable
    :param kwargs: Variables to help define dataframe columns
    :return:

    """

    source = get_value(kwargs, "source", "both")
    suffix_coordinates = get_value(kwargs, "suffix_coordinates", None)
    suffix_gene_sequence = get_value(kwargs, "suffix_gene_sequence", "gene-sequence")

    upstream_length_nt = get_value(kwargs, "upstream_length_nt", None)
    downstream_length_nt = get_value(kwargs, "downstream_length_nt", None)
    limit_upstream_to_first_candidate = get_value(kwargs, "limit_upstream_to_first_candidate", False)

    except_if_not_in_set(source, ["both", "q", "t"])

    # Code sketch:
    # 1) Find all needed genome sequences
    # 2) Map genome to indices in dataframe
    # 3) Read genomes one by one and add to dataframe

    q_labels_per_genome = sbsp_general.data.get_labels_per_genome(df, "query", coordinates_suffix=suffix_coordinates)
    t_labels_per_genome = sbsp_general.data.get_labels_per_genome(df, "target", coordinates_suffix=suffix_coordinates)

    # get nucleotide and protein labels
    q_sequences_per_genome = sbsp_io.sequences.read_sequences_from_genome_labels_pairs(
        env, q_labels_per_genome,
        leave_out_gene_stop=True,
        **kwargs
    )

    t_sequences_per_genome = sbsp_io.sequences.read_sequences_from_genome_labels_pairs(
        env, t_labels_per_genome,
        leave_out_gene_stop=True,
        **kwargs
    )

    sequences_per_genome = {
        "q": q_sequences_per_genome,
        "t": t_sequences_per_genome
    }

    sources = ["q", "t"]
    if source != "both":
        sources = [source]
    df["remove"] = False

    # add to dataframe
    for index, row in df.iterrows():

        # for each source
        for s in sources:

            # key for gene
            key = sbsp_general.general.create_gene_key(
                row["{}-{}".format(s, "genome")],
                row["{}-{}".format(s, "accession")],
                row["{}-{}".format(s, "left")],
                row["{}-{}".format(s, "right")],
                row["{}-{}".format(s, "strand")]
            )

            if key not in sequences_per_genome[s].keys():
                df.at[index, "remove"] = True
                # raise ValueError("Couldn't find sequence of key ({})".format(key))
                continue

            # for types of sequence
            for sequence_type in ["prot", "nucl"]:           # FIXME: change prot nucl to aa/nt

                if sequence_type in sequences_per_genome[s][key].keys():
                    df.at[index, "{}-{}-{}".format(s, sequence_type, suffix_gene_sequence)] = \
                        sequences_per_genome[s][key][sequence_type]
                    df.at[index, "{}-{}-pos-5prime-in-frag-{}".format(s, sequence_type, suffix_gene_sequence)] = \
                        sequences_per_genome[s][key]["{}-pos-5prime-in-frag".format(sequence_type)]

                    df.at[index, "{}-{}-position-of-5prime-in-msa-fragment-no-gaps".format(s, sequence_type)] = \
                        sequences_per_genome[s][key]["{}-pos-5prime-in-frag".format(sequence_type)]

    df = df[df["remove"] == False]

    return df


def compute_feature_helper(env, pf_data, **kwargs):
    # type: (Environment, str, Dict[str, Any]) -> pd.DataFrame
    # assumes sequences extracted

    df = pd.read_csv(pf_data, header=0)

    suffix_gene_sequence = get_value(kwargs, "suffix_gene_sequence", "gene-sequence")
    column_output = get_value(kwargs, "column_output", "k2p-distance")
    kimura_on_3rd = get_value(kwargs, "kimura_on_3rd", False)
    distance_types = get_value(kwargs, "distance_types", {"kimura"})

    df = df_add_labeled_sequences(env, df,
                                  source="both",
                                  suffix_gene_sequence=suffix_gene_sequence)

    matrix = matlist.blosum62
    import sbsp_alg.phylogeny
    sbsp_alg.phylogeny.add_stop_codon_to_blosum(matrix)

    df[column_output] = np.nan
    df["aa-match-fraction"] = np.nan

    for index, row in df.iterrows():
        # perform alignment of sequences
        q_sequence = row["q-prot-{}".format(suffix_gene_sequence)]
        t_sequence = row["t-prot-{}".format(suffix_gene_sequence)]

        q_sequence_nt = row["q-nucl-{}".format(suffix_gene_sequence)]
        t_sequence_nt = row["t-nucl-{}".format(suffix_gene_sequence)]

        [q_align, t_align, _, _, _] = \
            global_alignment_aa_with_gap(q_sequence, t_sequence, matrix)

        q_align_nt = add_gaps_to_nt_based_on_aa(q_sequence_nt, q_align)
        t_align_nt = add_gaps_to_nt_based_on_aa(t_sequence_nt, t_align)

        match_aa_frac = count_aa_mismatches(q_align, t_align)

        df.at[index, "aa-match-fraction"] = match_aa_frac

        try:
            df.at[index, column_output] = k2p_distance(q_align_nt, t_align_nt, kimura_on_3rd=kimura_on_3rd)
        except ValueError:
            df.at[index, column_output] = 100

        try:
            df.at[index, "kimura3"] = k2p_distance(q_align_nt, t_align_nt, kimura_on_3rd=True)
        except ValueError:
            df.at[index, "kimura3"] = 100

        for distance_type in distance_types:
            try:
                df.at[index, distance_type] = compute_distance(
                    distance_type,
                    q_align, t_align,
                    q_align_nt, t_align_nt,
                    on_fail=100,
                    **kwargs
                )
            except ValueError:
                df.at[index, distance_type] = 100

    return df

def compute_features(env, pf_data, pf_output, **kwargs):
    # type: (Environment, str, str, Dict[str, Any]) -> str

    pd_work = env["pd-work"]

    mkdir_p(pd_work)

    df = compute_feature_helper(env, pf_data)

    # clean up
    df.drop("q-nucl-gene-sequence", axis=1, inplace=True)
    df.drop("q-prot-gene-sequence", axis=1, inplace=True)
    df.drop("t-nucl-gene-sequence", axis=1, inplace=True)
    df.drop("t-prot-gene-sequence", axis=1, inplace=True)

    df.to_csv(pf_output, index=False)
    return pf_output
