import os

import numpy as np
import pandas as pd

import sbsp_general.data
import sbsp_general.labels
import sbsp_io.general
import sbsp_io.sequences
from sbsp_alg.phylogeny import k2p_distance, global_alignment_aa_with_gap
from sbsp_alg.gene_distances import *
from sbsp_general import Environment
from sbsp_general.general import get_value, except_if_not_in_set
from sbsp_options.msa import MSAOptions
from typing import *
# from memory_profiler import profile

from Bio.SubsMat import MatrixInfo as matlist

from sbsp_general.labels import Labels


def find_indices_for_column_value(df, column_name):
    # type: (pd.DataFrame, str) -> Dict[str, List[int]]

    """
    Retrieves the row indices for each unique value in the dataframe column
    :param df: a data frame
    :param column_name: the column
    :return: dictionary of list of indices
    """

    if column_name not in df.columns.values:
        raise ValueError("Invalid column name ({})".format(column_name))

    indices = dict()

    for index, row in df.iterrows():

        if row[column_name] not in indices.keys():
            indices[column_name] = list()

        indices[column_name].append(index)

    return indices


def _get_genome_to_source_to_indices(df, source):
    # type: (pd.DataFrame, str) -> Dict[str, Dict[str, List[int]]]

    def merge_and_add_source(merge_into, other, asource):
        # type: (Dict[str, Dict[str, List]], Dict[str, List], str) -> None

        if merge_into is None:
            raise ValueError("Can't merge into None")

        if other is None:
            return

        for k in other.keys():

            if k not in merge_into.keys():
                merge_into[k] = dict()

            if asource not in merge_into[k].keys():
                merge_into[k][asource] = list()

            merge_into[k][asource].append(other[k])

    # Structure:
    # d[genome_name][source] = list of indices
    genome_to_source_to_indices = dict()

    if source == "both":
        source = ["q", "t"]
    else:
        source = [source]

    for s in source:
        column_genome = "{}-{}".format(s, "genome")
        curr_dict = find_indices_for_column_value(df, column_genome)

        merge_and_add_source(genome_to_source_to_indices, curr_dict, s)

    return genome_to_source_to_indices

# @profile
def df_add_labeled_sequences(env, df, **kwargs):
    # type: (Environment, pd.DataFrame, Dict[str, Any]) -> pd.DataFrame

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
                raise ValueError("Couldn't find sequence of key ({})".format(key))

            # for types of sequence
            for sequence_type in ["prot", "nucl"]:           # FIXME: change prot nucl to aa/nt

                if sequence_type in sequences_per_genome[s][key].keys():
                    df.at[index, "{}-{}-{}".format(s, sequence_type, suffix_gene_sequence)] = \
                        sequences_per_genome[s][key][sequence_type]
                    df.at[index, "{}-{}-pos-5prime-in-frag-{}".format(s, sequence_type, suffix_gene_sequence)] = \
                        sequences_per_genome[s][key]["{}-pos-5prime-in-frag".format(sequence_type)]

                    df.at[index, "{}-{}-position-of-5prime-in-msa-fragment-no-gaps".format(s, sequence_type)] = \
                        sequences_per_genome[s][key]["{}-pos-5prime-in-frag".format(sequence_type)]

    return df


def df_collect_values_per_key(df, column_key, columns_values):
    # type: (pd.DataFrame, str, List[str]) -> Dict[str, Dict[str, List]]

    # make sure all column names exist
    for column_name in [column_key] + columns_values:
        if column_name not in df.columns.values:
            raise ValueError("Unknown column name ({})".format(column_name))

    result = dict()         # type: Dict[str, Dict[str, List]]

    for index, row in df.iterrows():

        key = row[column_key]

        if key not in result:
            result[key] = {c: list() for c in columns_values}
            result[key]["index"] = list()

        for c in columns_values:
            result[key][c].append(row[c])

        result[key]["index"].append(index)

    return result


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





def df_compute_kimura_helper(df, **kwargs):
    # type: (pd.DataFrame, **str) -> pd.DataFrame
    # assumes sequences extracted

    suffix_gene_sequence = get_value(kwargs, "suffix_gene_sequence", "gene-sequence")
    column_output = get_value(kwargs, "column_output", "k2p-distance")
    kimura_on_3rd = get_value(kwargs, "kimura_on_3rd", False)
    distance_types = get_value(kwargs, "distance_types", {"kimura"})

    matrix = matlist.blosum62
    import sbsp_alg.phylogeny
    from sbsp_alg.msa import add_gaps_to_nt_based_on_aa
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



def setup_simple_parallelization(computations_df, num_processes, func, *args, **kwargs):
    num_computations_needed = len(computations_df)
    computations_per_job = max(int(num_computations_needed / float(num_processes)), 1)

    job_vector = list()

    import multiprocessing as mp
    pool = mp.Pool(processes=num_processes)

    blocks = list()
    done = False
    block_number = 0

    while not done:
        start_block_idx = block_number * computations_per_job
        end_block_idx = start_block_idx + computations_per_job
        if end_block_idx > num_computations_needed:
            end_block_idx = num_computations_needed

        block = computations_df.iloc[start_block_idx:end_block_idx]
        blocks.append(block)

        block_number += 1

        # check if done
        if end_block_idx == num_computations_needed:
            done = True

    # take chunks of computations needed and do them separately
    counter = 0
    for block in blocks:

        args_with_block = (block,) + args

        job_vector.append(pool.apply_async(
            func,
            args=args_with_block,
            kwds=kwargs
        ))

        counter += 1

    results_vector = list()
    for p in job_vector:
        results_vector.append(p.get())

    df_result = pd.concat(results_vector)
    return df_result



def df_compute_kimura(env, df, **kwargs):

    suffix_coordinates = get_value(kwargs, "suffix_coordinates", None)
    suffix_gene_sequence = get_value(kwargs, "suffix_gene_sequence", "gene-sequence")
    column_k2p_distance = get_value(kwargs, "k2p_distance", "k2p-distance")
    num_processors = get_value(kwargs, "num_processors", 1)

    df = df_add_labeled_sequences(env, df,
                                  source="both",
                                  suffix_coordinates=suffix_coordinates,
                                  suffix_gene_sequence=suffix_gene_sequence)

    if num_processors == 1:

        df = df_compute_kimura_helper(df, suffix_gene_sequence=suffix_gene_sequence,
                                      column_output=column_k2p_distance, pd_work=env["pd-work"],
                                      **kwargs)

    else:
        df = setup_simple_parallelization(df, num_processors,
                                          df_compute_kimura_helper, suffix_gene_sequence=suffix_gene_sequence,
                                      column_output=column_k2p_distance, pd_work=env["pd-work"], **kwargs)

    # clean up
    df.drop("q-nucl-gene-sequence", axis=1, inplace=True)
    df.drop("q-prot-gene-sequence", axis=1, inplace=True)
    df.drop("t-nucl-gene-sequence", axis=1, inplace=True)
    df.drop("t-prot-gene-sequence", axis=1, inplace=True)

    return df


def df_filter(env, df, **kwargs):

    msa_options = get_value(kwargs, "msa_options", MSAOptions(env))     # type: (MSAOptions)
    pf_filter_stats = get_value(kwargs, "pf_filter_stats", None)
    filter_non_group_only = get_value(kwargs, "filter_non_group_only", True)


    suffix_coordinates = get_value(kwargs, "suffix_coordinates", None)
    tag_msa = get_value(kwargs, "tag_msa", "msa")

    upstream_length_nt = get_value(kwargs, "upstream_length_nt", None)
    downstream_length_nt = get_value(kwargs, "downstream_length_nt", None)

    from sbsp_alg.msa import filter_df, print_filter_stats_to_file, print_filter_stats

    column_distance = "k2p-distance"
    if msa_options.safe_get("column-distance"):
        column_distance = msa_options.safe_get("column-distance")

    filter_stats = {}
    df = filter_df(df, msa_options, filter_stats=filter_stats, filter_non_group_only=filter_non_group_only,
                   column_distance=column_distance)

    # df = df_add_labeled_sequences(env, df,
    #                               source="both",
    #                               suffix_coordinates=suffix_coordinates,
    #                               suffix_gene_sequence=tag_msa,
    #                               upstream_length_nt=upstream_length_nt,
    #                               downstream_length_nt=downstream_length_nt)

    if pf_filter_stats:
        print_filter_stats_to_file(filter_stats, pf_filter_stats)
    else:
        print_filter_stats(filter_stats)


    return df





def df_get_label_from_row(df, row_number, source, suffix_coordinates=None):
    # type: (pd.DataFrame, int, str, str) -> sbsp_general.labels.Label

    except_if_not_in_set(source, ["q", "t"])

    suffix = ""
    if suffix_coordinates is not None:
        suffix = "-{}".format(suffix_coordinates)

    return sbsp_general.labels.Label.from_fields({
        "left": df.at[row_number, "{}-{}{}".format(source, "left", suffix)] - 1,
        "right": df.at[row_number, "{}-{}{}".format(source, "right", suffix)] - 1,
        "strand": df.at[row_number, "{}-{}{}".format(source, "strand", suffix)],
        "seqname": df.at[row_number, "{}-{}".format(source, "accession")]
    })


def df_coordinates_to_row(df, row_number, label, source, suffix_coordinates=None):
    # type: (pd.DataFrame, int, sbsp_general.labels.Label, str, str) -> None

    except_if_not_in_set(source, ["q", "t"])

    suffix = ""
    if suffix_coordinates is not None:
        suffix = "-{}".format(suffix_coordinates)

    df.at[row_number, "{}-{}{}".format(source, "left", suffix)] = int(label.coordinates().left + 1)
    df.at[row_number, "{}-{}{}".format(source, "right", suffix)] = int(label.coordinates().right + 1)
    df.at[row_number, "{}-{}{}".format(source, "strand", suffix)] = label.coordinates().strand


def df_get_labels_per_genome(df, source, **kwargs):
    # type: (pd.DataFrame, str, **str) -> Dict[str, sbsp_general.labels.Labels]

    suffix_coordinates = get_value(kwargs, "suffix_coordinates", None)

    created_keys = set()            # to handle duplication

    labels_per_genome = dict()

    for index, row in df.iterrows():

        label = df_get_label_from_row(df, index, source, suffix_coordinates)
        genome = df.at[index, "{}-{}".format(source, "genome")]

        curr_key = sbsp_general.labels.create_gene_key_from_label(label)

        if curr_key not in created_keys:
            created_keys.add(curr_key)

            # add label
            if genome not in labels_per_genome.keys():
                labels_per_genome[genome] = sbsp_general.labels.Labels()

            labels_per_genome[genome].add(label)

    return labels_per_genome


def df_print_labels(env, df, source, **kwargs):
    # type: (Environment, pd.DataFrame, str, **str) -> Dict[str, str]

    labels_per_genome = df_get_labels_per_genome(df, source, **kwargs)
    suffix_fname = get_value(kwargs, "suffix_fname", "")
    if suffix_fname is None:
        suffix_fname = ""
    if suffix_fname != "":
        suffix_fname = "-" + suffix_fname

    genome_to_file = dict()

    for genome in labels_per_genome.keys():

        pf_curr = os.path.join(env['pd-work'], "{}.gff".format(genome))

        genome_to_file[genome] = pf_curr

        sbsp_io.general.write_string_to_file(labels_per_genome[genome].to_string(shift_coordinates_by=1), pf_curr)

    return genome_to_file


# @profile
def df_get_index_of_label_per_genome(env, df, source, **kwargs):
    # type: (dict, pd.DataFrame, str, **str) -> (Dict[str, Dict[str, Union[int]]], Dict)

    if source == "q" or source == "t":
        sources = [source]
    else:
        sources = ["q", "t"]

    index_of_label_per_genome = dict()  # type: Dict[str, Dict[str, int]]
    sorted_labels_per_genome = dict()  # type: Dict[str, List[sbsp_general.labels.Label]]

    for source in sources:

        genomes = set(df["{}-{}".format(source, "genome")])

        # read all labels for each genome
        import sbsp_io.labels
        labels_per_genome = dict()
        for g in genomes:
            labels_per_genome[g] = sbsp_io.labels.read_labels_from_file(
                os.path.join(env["pd-data"], g, "ncbi.gff")
            )

        for genome in labels_per_genome.keys():

            sorted_labels_per_genome[genome] = list()

            index_of_label_per_genome[genome] = dict()

            for label in labels_per_genome[genome]:
                sorted_labels_per_genome[genome].append(label)

            # sort labels list by left
            sorted_labels_per_genome[genome].sort(key=lambda r: r.coordinates().left)

            # get index for each label
            for i in range(len(sorted_labels_per_genome[genome])):
                label = sorted_labels_per_genome[genome][i]
                # label_key = sbsp_general.labels.create_gene_key_from_label(label)
                label_key = sbsp_general.general.create_3prime_key(accession=label.seqname(),
                                                                   left=label.coordinates().left,
                                                                   right=label.coordinates().right,
                                                                   strand=label.strand())
                index_of_label_per_genome[genome][label_key] = i

    return [index_of_label_per_genome, sorted_labels_per_genome]


# def df_add_distance_to_upstream_gene(env, df, source,
#                                      index_of_label_per_genome,
#                                      sorted_labels_per_genome,
#                                      **kwargs):
#     # type: (Dict[str, Any], pd.DataFrame, str, Dict[str, Any], Dict[str, Any], Dict[str, Any]) -> None
#
#     suffix_coordinates = get_value(kwargs, "suffix_corrdinates", None)
#     suffix_upstream_distance = get_value(kwargs, "suffix_upstream_distance", "upstream-distance")
#
#     # for each dataframe row, find the label and check upstream
#     inf = 10000000
#     for index, row in df.iterrows():
#
#         label = df_get_label_from_row(df, index, source, suffix_coordinates)
#         # label_key = sbsp_general.labels.create_gene_key_from_label(label)
#         label_key = sbsp_general.general.create_3prime_key(accession=label.seqname(),
#                                                            left=label.coordinates().left,
#                                                            right=label.coordinates().right,
#                                                            strand=label.strand())
#         genome = row["{}-genome".format(source)]
#         pos = index_of_label_per_genome[genome][label_key]
#
#         distance = inf
#
#         # positive strand, check before
#         if label.strand() == "+":
#             if pos != 0:
#                 upstream_label = sorted_labels_per_genome[genome][pos - 1]
#
#                 # compute distance to upstream
#                 distance = label.coordinates().left - upstream_label.coordinates().right
#
#         # negative strand, check after
#         elif label.strand() == "-":
#
#             if pos < len(sorted_labels_per_genome[genome]) - 1:
#                 upstream_label = sorted_labels_per_genome[genome][pos + 1]
#
#                 # compute distance to upstream
#                 distance = upstream_label.coordinates().left - label.coordinates().right
#
#         df.at[index, "{}-{}".format(source, suffix_upstream_distance)] = distance




def df_add_position_of_upstream_gene_in_msa_no_gaps(df):
    # type: (pd.DataFrame) -> None

    for s in ["q", "t"]:

        df["{}-nucl-position-of-upstream-gene-in-msa-fragment-no-gaps".format(s)] = int(0)

        for index, row in df.iterrows():

            position_of_5prime_in_msa = row["{}-nucl-position-of-5prime-in-msa-fragment-no-gaps".format(s)]
            distance_nt_to_upstream_gene = row["{}-nucl-distance-to-upstream-gene".format(s)]

            pos = position_of_5prime_in_msa - distance_nt_to_upstream_gene
            df.at[index, "{}-nucl-position-of-upstream-gene-in-msa-fragment-no-gaps".format(s)] = pos


def get_reference_labels_per_genome(env, genome_names):
    # type: (Dict[str, Any], Iterable[str]) -> Dict[str, Labels]

    from sbsp_io.labels import read_labels_from_file

    labels_per_genome = {
        g: read_labels_from_file(os.path.join(env["pd-data"], g, "ncbi.gff")) for g in genome_names
    }

    return labels_per_genome

def get_labels_per_genome_from_df(df, source):
    # type: (pd.DataFrame, str) -> Dict[str, Labels]

    labels_per_genome = dict()

    for index, row in df.iterrows():
        genome_name = row["{}-genome".format(source)]

        if genome_name not in labels_per_genome:
            labels_per_genome[genome_name] = Labels(name=genome_name)

        labels_per_genome[genome_name].add(
            df_get_label_from_row(df, index, source)
        )

    return labels_per_genome

def get_distance_to_upstream_gene(reference_labels, query_labels):
    # type: (Labels, Labels) -> Dict[str, int]

    merged_labels = reference_labels.update(query_labels)

    return merged_labels.compute_distance_to_upstream_gene()


def get_distance_to_upstream_gene_per_genome(reference_labels_per_genome, labels_per_genome_from_df):
    # type: (Dict[str, Labels], Dict[str, Labels]) -> Dict[str, Dict[str, int]]

    genome_names = reference_labels_per_genome.keys()

    distance_to_upstream_per_genome = dict()

    for g in genome_names:
        distance_to_upstream_per_genome[g] = \
            get_distance_to_upstream_gene(reference_labels_per_genome[g], labels_per_genome_from_df[g])

    return distance_to_upstream_per_genome

def df_add_distance_to_upstream_gene(env, df, source, **kwargs):
    # type: (Dict[str, Any], pd.DataFrame, str, Dict[str, Any]) -> None

    suffix_coordinates = get_value(kwargs, "suffix_corrdinates", None)
    suffix_upstream_distance = get_value(kwargs, "suffix_upstream_distance", "upstream-distance")

    # get all reference labels per genome for given source
    genome_names = set(df["{}-{}".format(source, "genome")])
    reference_labels_per_genome = get_reference_labels_per_genome(env, genome_names)

    # add labels from df to this set, and overwrite those with same 3prime end
    labels_per_genome_from_df = get_labels_per_genome_from_df(df, source)

    # get distance to upstream gene for each label, indexed by key
    distance_to_upstream_gene_per_genome = get_distance_to_upstream_gene_per_genome(
        reference_labels_per_genome, labels_per_genome_from_df
    )

    # update df
    for index, row in df.iterrows():
        label = df_get_label_from_row(df, index, source, suffix_coordinates)
        genome = row["{}-genome".format(source)]

        label_key = sbsp_general.labels.create_gene_key_from_label(label)

        distance = distance_to_upstream_gene_per_genome[genome][label_key]

        if distance is None:
            distance = 1000000

        df.at[index, "{}-{}".format(source, suffix_upstream_distance)] = distance
        df.at[index, "{}-nucl-distance-to-upstream-gene".format(source)] = distance


def df_add_distance_to_upstream_gene_DEPRECATED(env, df, source, **kwargs):
    # type: (dict, pd.DataFrame, str, **str) -> None

    suffix_coordinates = get_value(kwargs, "suffix_corrdinates", None)
    suffix_upstream_distance = get_value(kwargs, "suffix_upstream_distance", "upstream-distance")

    labels_info = get_value(kwargs, "labels_info", None)
    if labels_info is None:

        genomes = set(df["{}-{}".format(source, "genome")])

        # read all labels for each genome
        import sbsp_io.labels
        labels_per_genome = dict()
        for g in genomes:
            labels_per_genome[g] = sbsp_io.labels.read_labels_from_file(
                os.path.join(env["pd-data"], g, "ncbi.gff")
            )

            # labels_per_genome[g].sort_by("left", in_place=True)

        index_of_label_per_genome = dict()          # type: Dict[str, Dict[str, int]]
        sorted_labels_per_genome = dict()           # type: Dict[str, List[sbsp_general.labels.Label]]

        for genome in labels_per_genome.keys():

            sorted_labels_per_genome[genome] = list()

            index_of_label_per_genome[genome] = dict()

            for label in labels_per_genome[genome]:
                sorted_labels_per_genome[genome].append(label)

            # sort labels list by left
            sorted_labels_per_genome[genome].sort(key=lambda r: r.coordinates().left)

            # get index for each label
            for i in range(len(sorted_labels_per_genome[genome])):
                label = sorted_labels_per_genome[genome][i]
                # label_key = sbsp_general.labels.create_gene_key_from_label(label)
                label_key = sbsp_general.general.create_3prime_key(accession=label.seqname(),
                                                                   left=label.coordinates().left,
                                                                   right=label.coordinates().right,
                                                                   strand=label.strand())
                index_of_label_per_genome[genome][label_key] = i
    else:
        [index_of_label_per_genome, sorted_labels_per_genome] = labels_info


    # for each dataframe row, find the label and check upstream
    inf = 10000000
    for index, row in df.iterrows():

        label = df_get_label_from_row(df, index, source, suffix_coordinates)
        # label_key = sbsp_general.labels.create_gene_key_from_label(label)
        label_key = sbsp_general.general.create_3prime_key(accession=label.seqname(),
                                                           left=label.coordinates().left,
                                                           right=label.coordinates().right,
                                                           strand=label.strand())
        genome = row["{}-genome".format(source)]
        pos = index_of_label_per_genome[genome][label_key]

        distance = inf

        # positive strand, check before
        if label.strand() == "+":
            if pos != 0:
                upstream_label = sorted_labels_per_genome[genome][pos-1]

                # compute distance to upstream
                distance = label.coordinates().left - upstream_label.coordinates().right

        # negative strand, check after
        elif label.strand() == "-":

            if pos < len(sorted_labels_per_genome[genome])-1:

                upstream_label = sorted_labels_per_genome[genome][pos+1]

                # compute distance to upstream
                distance = upstream_label.coordinates().left - label.coordinates().right

        df.at[index, "{}-{}".format(source, suffix_upstream_distance)] = distance
        df.at[index, "{}-nucl-distance-to-upstream-gene".format(source)] = distance
