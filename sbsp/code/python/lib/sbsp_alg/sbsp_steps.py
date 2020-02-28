import copy
import os
import logging
import time
import timeit
from multiprocessing import Process, Manager, Lock
from typing import *

from Bio.Align import MultipleSeqAlignment
from Bio.Align.Applications import ClustalwCommandline, ClustalOmegaCommandline
from Bio.Blast import NCBIXML, Record
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SubsMat import MatrixInfo as matlist
from numpy import isclose

import sbsp_io
import sbsp_ml
import sbsp_ml.msa_features_2
from sbsp_alg.feature_computation import compute_features
from sbsp_alg.filtering import filter_orthologs
from sbsp_alg.msa import run_sbsp_msa, get_files_per_key, run_sbsp_msa_from_multiple, \
    run_sbsp_msa_from_multiple_for_multiple_queries, perform_msa_on_df, move_files_using_scp, should_count_in_neighbor, \
    filter_by_pairwise_kimura_from_msa
from sbsp_general.blast import run_blast
from sbsp_general.labels import Label, Coordinates
from sbsp_general.msa_2 import MSAType, MSASinglePointMarker
from sbsp_io.general import mkdir_p, remove_p, write_string_to_file
from sbsp_general.general import get_value, except_if_not_in_set
from sbsp_alg.ortholog_finder import get_orthologs_from_files, extract_labeled_sequences_for_genomes, \
    unpack_fasta_header, select_representative_hsp, create_info_for_query_target_pair, \
    compute_distance_based_on_local_alignment, compute_distance_based_on_global_alignment_from_sequences, \
    run_blast_on_sequence_file, is_valid_start
from sbsp_alg.sbsp_compute_accuracy import pipeline_step_compute_accuracy, separate_msa_outputs_by_stats
from sbsp_general import Environment
from sbsp_io.general import read_rows_to_list
from sbsp_io.msa_2 import add_true_starts_to_msa_output
from sbsp_io.sequences import read_fasta_into_hash, write_fasta_hash_to_file
from sbsp_ml.msa_features_2 import ScoringMatrix
from sbsp_options.sbsp import SBSPOptions
from sbsp_options.pbs import PBSOptions
from sbsp_options.pipeline_sbsp import PipelineSBSPOptions
from sbsp_parallelization.pbs import PBS
from sbsp_pbs_data.mergers import merge_identity
from sbsp_pbs_data.splitters import *

logger = logging.getLogger(__name__)


def duplicate_pbs_options_with_updated_paths(env, pbs_options, **kwargs):
    # type: (Environment, PBSOptions, Dict[str, Any]) -> PBSOptions
    keep_on_head = get_value(kwargs, "keep_on_head", False, default_if_none=True)

    pbs_options = copy.deepcopy(pbs_options)
    pbs_options["pd-head"] = os.path.abspath(env["pd-work"])

    if keep_on_head:
        pbs_options["pd-root-compute"] = os.path.abspath(env["pd-work"])
    elif pbs_options["pd-root-compute"] is None:
        pbs_options["pd-root-compute"] = os.path.abspath(env["pd-work"])

    return pbs_options


def run_step_generic(env, pipeline_options, step_name, splitter, merger, data, func, func_kwargs, **kwargs):
    # type: (Environment, PipelineSBSPOptions, str, Callable, Callable, Dict[str, Any], Callable, Dict[str, Any], Dict[str, Any]) -> Dict[str, Any]

    output = {
        "pf-list-output": os.path.join(env["pd-work"], "pbs-summary.txt")
    }

    if pipeline_options.use_pbs():
        pbs_options = duplicate_pbs_options_with_updated_paths(env, pipeline_options["pbs-options"])

        if pbs_options.safe_get("pd-data-compute"):
            env = env.duplicate({"pd-data": pbs_options["pd-data-compute"]})

        pbs = PBS(env, pbs_options,
                  splitter=splitter,
                  merger=merger
                  )

        if pipeline_options.perform_step(step_name):

            output = pbs.run(
                data=data,
                func=func,
                func_kwargs=func_kwargs
            )
        else:
            # read data from file
            list_pf_output_packages = read_rows_to_list(os.path.join(env["pd-work"], "pbs-summary.txt"))
            output = pbs.merge_output_package_files(list_pf_output_packages)

    return output




def sbsp_step_get_orthologs(env, pipeline_options):
    # type: (Environment, PipelineSBSPOptions) -> Dict[str, Any]
    """
    Given a list of query and target genomes, find the set of related genes
    for each query
    """

    logger.debug("Running: sbsp_step_get_orthologs")

    output = {
        "pf-list-output": os.path.join(env["pd-work"], "pbs-summary.txt")
    }

    if pipeline_options.use_pbs():

        pbs_options = duplicate_pbs_options_with_updated_paths(env, pipeline_options["pbs-options"])

        if pbs_options.safe_get("pd-data-compute"):
            env = env.duplicate({"pd-data": pbs_options["pd-data-compute"]})

        pbs = PBS(env, pbs_options,
                  splitter=split_query_genomes_target_genomes_one_vs_group,
                  merger=merge_identity
        )

        if pipeline_options.perform_step("get-orthologs"):
            output = pbs.run(
                data={"pf_q_list": pipeline_options["pf-q-list"],
                      "pf_output_template": os.path.join(pbs_options["pd-head"], pipeline_options["fn-orthologs"] + "_{}")},
                func=get_orthologs_from_files,
                func_kwargs={
                    "env": env,
                    "pf_t_db": pipeline_options["pf-t-db"],
                    "fn_q_labels": pipeline_options["fn-q-labels"], "fn_t_labels": pipeline_options["fn-t-labels"],
                    "max_evalue": pipeline_options.safe_get("max-evalue"),
                    "clean": True,
                    "sbsp_options": pipeline_options.safe_get("msa-options")
                }
            )
        else:
            # read data from file
            list_pf_output_packages = read_rows_to_list(os.path.join(env["pd-work"], "pbs-summary.txt"))
            output = pbs.merge_output_package_files(list_pf_output_packages)


    return output


def sbsp_step_compute_features(env, pipeline_options, list_pf_previous):
    # type: (Environment, PipelineSBSPOptions, List[str]) -> Dict[str, Any]
    """
    Given a list of query and target genomes, find the set of related genes
    for each query
    """

    logger.debug("Running: sbsp_step_compute_features")

    output = {
        "pf-list-output": os.path.join(env["pd-work"], "pbs-summary.txt")
    }

    if pipeline_options.use_pbs():

        pbs_options = duplicate_pbs_options_with_updated_paths(env, pipeline_options["pbs-options"])

        if pbs_options.safe_get("pd-data-compute"):
            env = env.duplicate({"pd-data": pbs_options["pd-data-compute"]})

        pbs = PBS(env, pbs_options,
                  splitter=split_list,
                  merger=merge_identity
                  )

        if pipeline_options.perform_step("compute-features"):

            output = pbs.run(
                data={"list_pf_data": list_pf_previous,
                      "pf_output_template": os.path.join(pbs_options["pd-head"],
                                                         pipeline_options["fn-compute-features"] + "_{}")},
                func=compute_features,
                func_kwargs={
                    "env": env,
                    "clean": True
                }
            )
        else:
            # read data from file
            list_pf_output_packages = read_rows_to_list(os.path.join(env["pd-work"], "pbs-summary.txt"))
            output = pbs.merge_output_package_files(list_pf_output_packages)

    return output

def sbsp_step_filter(env, pipeline_options, list_pf_previous):
    # type: (Environment, PipelineSBSPOptions, List[str]) -> Dict[str, Any]
    """
    Given a list of query and target genomes, find the set of related genes
    for each query
    """

    logger.debug("Running: sbsp_step_filter")

    output = {
        "pf-list-output": os.path.join(env["pd-work"], "pbs-summary.txt")
    }

    if pipeline_options.use_pbs():
        pbs_options = duplicate_pbs_options_with_updated_paths(env, pipeline_options["pbs-options"])

        if pbs_options.safe_get("pd-data-compute"):
            env = env.duplicate({"pd-data": pbs_options["pd-data-compute"]})

        pbs = PBS(env, pbs_options,
                  splitter=split_list,
                  merger=merge_identity
                  )


        if pipeline_options.perform_step("filter"):

            output = pbs.run(
                data={"list_pf_data": list_pf_previous,
                      "pf_output_template": os.path.join(pbs_options["pd-head"],
                                                         pipeline_options["fn-filter"] + "_{}")},
                func=filter_orthologs,
                func_kwargs={
                    "env": env,
                    "msa_options": pipeline_options["msa-options"],
                    "clean": True
                }
            )
        else:
            # read data from file
            list_pf_output_packages = read_rows_to_list(os.path.join(env["pd-work"], "pbs-summary.txt"))
            output = pbs.merge_output_package_files(list_pf_output_packages)

    return output


def sbsp_step_msa(env, pipeline_options, list_pf_previous):
    # type: (Environment, PipelineSBSPOptions, List[str]) -> Dict[str, Any]
    """
    Given a list of query and target genomes, find the set of related genes
    for each query
    """

    logger.debug("Running: sbsp_step_msa")

    output = {
        "pf-list-output": os.path.join(env["pd-work"], "pbs-summary.txt")
    }

    if pipeline_options.use_pbs():
        pbs_options = duplicate_pbs_options_with_updated_paths(env, pipeline_options["pbs-options"], keep_on_head=False)

        if pbs_options.safe_get("pd-data-compute"):
            env = env.duplicate({"pd-data": pbs_options["pd-data-compute"]})

        # pbs = PBS(env, pbs_options,
        #           splitter=split_list_and_remerge_by_key,
        #           merger=merge_identity
        #           )

        pbs = PBS(env, pbs_options,
                  splitter=split_q3prime_to_list_of_data_files,
                  merger=merge_identity
                  )

        if pipeline_options.perform_step("build-msa"):

            # get files per 3prime key
            q3prime_to_list_pf = get_files_per_key(list_pf_previous)
            pd_msa = os.path.join(pbs_options["pd-head"], "msa")
            mkdir_p(pd_msa)


            output = pbs.run(
                data={"q3prime_to_list_pf": q3prime_to_list_pf,
                      "pf_output_template": os.path.join(pbs_options["pd-head"],
                                                         pipeline_options["fn-msa"] + "_{}")},
                func=run_sbsp_msa_from_multiple_for_multiple_queries,
                func_kwargs={
                    "env": env,
                    "msa_options": pipeline_options["msa-options"],
                    "clean": True,
                    "pd_msa_final": pd_msa
                }
            )


            # output = pbs.run(
            #     data={"list_pf_data": list_pf_previous, "group_key": "q-3prime",
            #           "pf_output_template": os.path.join(pbs_options["pd-head"],
            #                                              pipeline_options["fn-msa"] + "_{}")},
            #     func=run_sbsp_msa,
            #     func_kwargs={
            #         "env": env,
            #         "msa_options": pipeline_options["msa-options"]
            #     }
            # )
        else:
            # read data from file
            list_pf_output_packages = read_rows_to_list(os.path.join(env["pd-work"], "pbs-summary.txt"))
            output = pbs.merge_output_package_files(list_pf_output_packages)

    return output


def run_blast_on_sequences(env, q_sequences, pf_t_db, pf_blast_output, sbsp_options, **kwargs):
    # type: (Environment, Dict[str, Seq], str, str, SBSPOptions, Dict[str, Any]) -> None

    fn_tmp_prefix = get_value(kwargs, "fn_tmp_prefix", None)

    # write sequences to a file
    pf_q_sequences = os.path.join(env["pd-work"], "{}query_sequences.fasta".format(fn_tmp_prefix))
    write_fasta_hash_to_file(q_sequences, pf_q_sequences)

    # start clean
    remove_p(pf_blast_output)

    try:
        run_blast_on_sequence_file(env, pf_q_sequences, pf_t_db, pf_blast_output, sbsp_options=sbsp_options)
    except ValueError:
        raise ValueError("Couldn't run blast")


def create_data_frame_for_msa_search_from_blast_results(r, sbsp_options, **kwargs):
    # type: (Record, SBSPOptions, Dict[str, Any]) -> pd.DataFrame

    df = pd.DataFrame()
    query_info = unpack_fasta_header(r.query)

    # if query_info["left"] != 178270:
    #     logger.debug("{} != 178270".format(query_info["left"]))
    #     return pd.DataFrame()

    distance_min = sbsp_options.safe_get("distance-min")
    distance_max = sbsp_options.safe_get("distance-max")
    max_targets = 50

    filter_orthologs_with_equal_kimura_to_query = sbsp_options.safe_get("filter-orthologs-with-equal-kimura")
    set_target_kimuras = set()

    # for each alignment to a target protein for the current query
    list_entries = list()
    logger.debug("Reading {} targets from blast".format(len(r.alignments)))
    for alignment in r.alignments:
        if len(list_entries) > max_targets:
            logger.debug("Reached limit on number of targets: {} from {}".format(max_targets, len(r.alignments)))
            break

        target_info = unpack_fasta_header(alignment.title)
        hsp = select_representative_hsp(alignment, "")  # get reference hit for target

        # get nucleotide sequence that corresponds to proteins
        original_q_nt = Seq(query_info["lorf_nt"][query_info["offset"]:])
        original_t_nt = Seq(target_info["lorf_nt"][target_info["offset"]:])

        distance = compute_distance_based_on_local_alignment(query_info, target_info, hsp,
                                                             original_q_nt=original_q_nt,
                                                             original_t_nt=original_t_nt,
                                                             **kwargs)


        if distance > distance_min and distance < distance_max:
            if filter_orthologs_with_equal_kimura_to_query is not None:
                rounded = round(distance, filter_orthologs_with_equal_kimura_to_query)
                # if another target with equal (to decimal place) exists, skip this one.
                if rounded in set_target_kimuras:
                    logger.debug("Filtered due to equal Kimura")
                    continue
                else:
                    set_target_kimuras.add(rounded)


            output_info = create_info_for_query_target_pair(
                query_info, target_info, hsp,
                distance_blast=distance,
                distance=distance,
                local_distance=distance,
                local_length=hsp.align_length,
                local_length_without_gaps=sum([
                    1 for i in range(len(hsp.query)) if hsp.query[i] != "-" and hsp.sbjct[i] != "-"
                ])
            )

            list_entries.append(output_info)

    if len(list_entries) > 0:
        df = pd.DataFrame(list_entries)

    return df

def extract_sequences_from_df_for_msa(df):
    # type: (pd.DataFrame) -> List[Seq]

    list_sequences = list()

    if len(df) > 0:

        list_sequences.append(Seq(df.iloc[0]["q-lorf_nt"]).translate())

        for _, row in df.iterrows():
            list_sequences.append(Seq(row["t-lorf_nt"]).translate())

    return list_sequences


def write_sequence_list_to_fasta_file(sequences, pf_sequences):
    # type: (List[Seq], str) -> None

    data = ""
    for i in range(len(sequences)):
        data += ">{}\n{}\n".format(i, sequences[i])

    write_string_to_file(data, pf_sequences)


def run_msa_on_sequence_file(pf_fasta, sbsp_options, pf_msa, **kwargs):
    # type: (str, SBSPOptions, str, Dict[str, Any]) -> None

    gapopen = sbsp_options.safe_get("msa-gapopen")
    gapext = sbsp_options.safe_get("msa-gapext")

    num_processors = get_value(kwargs, "num_processors", None)

    # clustalw_cline = ClustalwCommandline(
    #     "clustalw2", infile=pf_fasta, outfile=pf_msa,
    #     gapopen=gapopen,
    #     gapext=gapext,
    #     outorder="input"
    # )

    logger.debug("Number of processors for MSA: {}".format(num_processors))
    clustalw_cline = ClustalOmegaCommandline(
        "clustalo", infile=pf_fasta, outfile=pf_msa,
        #gapopen=gapopen,
        #gapext=gapext,
        outputorder="input-order",
        force=True,
        outfmt="clustal",
        threads=num_processors
    )

    clustalw_cline()


def run_msa_on_sequences(env, sequences, sbsp_options, **kwargs):
    # type: (Environment, List[Seq], SBSPOptions, Dict[str, Any]) -> MSAType

    pd_work = env["pd-work"]
    fn_tmp_prefix = get_value(kwargs, "fn_tmp_prefix", "", default_if_none=True)

    # write sequences to file
    pf_fasta = os.path.join(pd_work, "{}tmp_sequences.fasta".format(fn_tmp_prefix))
    remove_p(pf_fasta)
    write_sequence_list_to_fasta_file(sequences, pf_fasta)

    # run msa
    pf_msa = os.path.join(pd_work, "{}tmp_msa.txt".format(fn_tmp_prefix))
    run_msa_on_sequence_file(pf_fasta, sbsp_options, pf_msa, **kwargs)

    msa_t = MSAType.init_from_file(pf_msa)

    remove_p(pf_msa, pf_fasta)

    return msa_t


def convert_gapped_aa_to_gapped_nt(seq_aa, seq_nt_no_gaps):
    # type: (Seq, str) -> Seq

    seq_nt_with_gaps = ""
    pos_no_gaps = 0
    for i in range(len(seq_aa)):
        if seq_aa[i] == "-":
            seq_nt_with_gaps += "---"
        else:
            seq_nt_with_gaps += seq_nt_no_gaps[pos_no_gaps:pos_no_gaps+3]
            pos_no_gaps += 3

    return Seq(seq_nt_with_gaps)


def convert_msa_aa_to_nt(msa_t_aa, df):
    # type: (MSAType, pd.DataFrame) -> MSAType

    seq_record_list = list()
    # query sequence
    seq_record_list.append(SeqRecord(convert_gapped_aa_to_gapped_nt(msa_t_aa[0].seq, df.iloc[0]["q-lorf_nt"])))

    # targets
    for i in range(1, msa_t_aa.number_of_sequences()):
        row = df.iloc[i-1]          # -1 to account for query as first sequence
        seq_record_list.append(SeqRecord(convert_gapped_aa_to_gapped_nt(msa_t_aa[i].seq, row["t-lorf_nt"])))

    return MSAType(MultipleSeqAlignment(seq_record_list))


def lower_case_non_5prime_in_msa(msa_t_aa, msa_t_nt):
    # type: (MSAType, MSAType) -> MSAType

    seq_record_list = list()

    for i in range(msa_t_aa.number_of_sequences()):

        new_seq_aa = ""
        for j_aa in range(msa_t_aa.alignment_length()):
            j_nt = j_aa * 3

            if is_valid_start(msa_t_nt[i].seq._data[j_nt:j_nt+3], "+"):
                new_seq_aa += msa_t_aa[i].seq._data[j_aa].upper()
            else:
                new_seq_aa += msa_t_aa[i].seq._data[j_aa].lower()

        seq_record_list.append(SeqRecord(Seq(new_seq_aa), id=msa_t_aa[i].id))

    return MSAType(MultipleSeqAlignment(seq_record_list))




def construct_msa_from_df(env, df, sbsp_options, **kwargs):
    # type: (Environment, pd.DataFrame, SBSPOptions, Dict[str, Any]) -> Tuple[Union[None, MSAType], Union[None, MSAType]]

    # extract sequences
    sequences = extract_sequences_from_df_for_msa(df)
    if len(sequences) == 0:
        return None, None

    msa_t_aa = run_msa_on_sequences(env, sequences, sbsp_options, **kwargs)

    # get nucleotide version of msa
    msa_t_nt = convert_msa_aa_to_nt(msa_t_aa, df)

    msa_t_aa = lower_case_non_5prime_in_msa(msa_t_aa, msa_t_nt)

    return msa_t_aa, msa_t_nt

def number_of_sequences_with_gap_in_position(msa_t, pos):
    # type: (MSAType, int) -> int

    if msa_t.number_of_sequences() == 0:
        return 0

    return sum(1 for i in range(msa_t.number_of_sequences()) if msa_t[i][pos] == "-")

def get_position_from_which_to_start_gap_filtering(msa_t):
    # type: (MSAType) -> int

    # get to first start codon in the query and make sure most targets have reached that point
    first_start_codon_position = None
    passed_go = False
    for i in range(msa_t.alignment_length()):
        if msa_t[0][i].isupper():
            first_start_codon_position = i

        if number_of_sequences_with_gap_in_position(msa_t, i) / float(msa_t.number_of_sequences()) < 0.3:
            passed_go = True

        if passed_go and first_start_codon_position is not None:
            break

    if first_start_codon_position is None:
        first_start_codon_position = 0

    return first_start_codon_position


def filter_df_based_on_msa(df, msa_t, msa_t_nt, sbsp_options, inplace=False, multiple_filterings=False):
    # type: (pd.DataFrame, MSAType, MSAType, SBSPOptions, bool, bool) -> pd.DataFrame

    if not inplace:
        df = df.copy()

    row_numbers_to_drop = set()

    # pairwise Kimura
    if sbsp_options.safe_get("filter-by-pairwise-kimura-from-msa"):
        _, indices_to_keep = filter_by_pairwise_kimura_from_msa(
            [msa_t_nt[i].seq._data for i in range(msa_t_nt.number_of_sequences())], sbsp_options
        )

        indices_in_msa_to_remove = set(range(msa_t_nt.number_of_sequences())).difference(indices_to_keep)
        for i in indices_in_msa_to_remove:
            row_numbers_to_drop.add(i-1)


    params = sbsp_options.safe_get("filter-remove-sequences-that-introduce-gaps")
    gap_width = params[0]
    seq_frac = params[1]

    alignment_length = msa_t.alignment_length()
    num_sequences_aligned = msa_t.number_of_sequences()

    # get to first start codon in the query and make sure most targets have reached that point
    first_start_codon_position = get_position_from_which_to_start_gap_filtering(msa_t)



    end_search = int(alignment_length / 2.0) - gap_width
    for i in range(first_start_codon_position, end_search):

        # if chunk of gaps detected in query
        block_detected = msa_t[0][i:i + gap_width].seq._data.count("-") == gap_width

        sequences_that_contribute_to_block = list()

        if block_detected:

            # find sequences that have no gaps in that region
            for j in range(1, num_sequences_aligned):
                if msa_t[j][i:i + gap_width].seq._data.count("-") == 0:
                    sequences_that_contribute_to_block.append(j-1)

            # compute fraction of these sequences
            num_sequences_that_contribute_to_block = len(sequences_that_contribute_to_block)
            fraction = float(num_sequences_that_contribute_to_block) / num_sequences_aligned

            # remove those sequences if they are very few
            if fraction != 0 and fraction < seq_frac:
                for s in sequences_that_contribute_to_block:
                    row_numbers_to_drop.add(s)

                if not multiple_filterings:
                    break

    # remove any rows necessary
    if len(row_numbers_to_drop) > 0:
        df.drop(df.index[list(row_numbers_to_drop)], inplace=True)

    return df

def get_next_position_in_msa(l_curr_pos, l_msa_t, l_direction, l_skip_gaps_in_query=False):
    # type: (int, MSAType, str, bool) -> Union[int, None]

    if l_direction == "upstream":

        new_pos = l_curr_pos - 1
        while True:
            if new_pos < 0:
                return None
            if l_skip_gaps_in_query and l_msa_t[0][new_pos] == "-":
                new_pos = new_pos - 1
                continue
            break
    else:
        new_pos = l_curr_pos + 1
        while True:
            if new_pos >= l_msa_t.alignment_length():
                return None
            if l_skip_gaps_in_query and l_msa_t[0][new_pos] == "-":
                new_pos = new_pos + 1
                continue
            break
    return new_pos

def compute_conservation_in_region(msa_t, start, end, scorer, **kwargs):
    # type: (MSAType, int, int, ScoringMatrix, Dict[str, Any]) -> float
    """

    :param msa_t:
    :param start:
    :param end: (exclusive)
    :param scorer:
    :param kwargs:
    :return:
    """

    direction = get_value(kwargs, "direction", choices=["upstream", "downstream"], default="downstream")
    score_on_all_pairs = get_value(kwargs, "score_on_all_pairs", False)
    skip_gaps_in_query = get_value(kwargs, "skip_gaps_in_query", False)

    if start < 0 or start >= msa_t.alignment_length():
        raise ValueError("Start of region out of bounds: {} not in [{},{}]".format(start, 0, msa_t.alignment_length()))
    if end < start or end > msa_t.alignment_length():
        raise ValueError("Start of region out of bounds: {} not in [{},{})".format(end, start, msa_t.alignment_length()))



    num_positions_to_analyze = end - start

    curr_pos = start if direction == "downstream" else end - 1
    total_number_of_computations = 0
    pos_score = 0
    num_sequences = msa_t.number_of_sequences()

    for n in range(num_positions_to_analyze):

        if curr_pos == None:
            raise ValueError("Not enough region to compute score")


        if not score_on_all_pairs:
            for i in range(num_sequences):
                element_i = msa_t[i][curr_pos]
                for j in range(i+1, num_sequences):
                    element_j = msa_t[j][curr_pos]

                    pos_score += scorer.score(element_i, element_j)
                    total_number_of_computations += 1
        else:
            for i in range(1, num_sequences):
                pos_score += scorer.score(msa_t[0][curr_pos], msa_t[i][curr_pos])
                total_number_of_computations += 1

        curr_pos = get_next_position_in_msa(curr_pos, msa_t, direction, skip_gaps_in_query)


    return pos_score / float(total_number_of_computations)







def get_all_candidates_before_conserved_block(msa_t, sbsp_options, at_least_until=0):
    # type: (MSAType, SBSPOptions, int) -> List[int]

    logger.debug("Func: get-candidates-without-upstream-conservation")

    region_length = sbsp_options["block-region-length-aa"]     # get_value(kwargs, "region_length", 10)
    threshold = 0.5         # FIXME
    score_on_all_pairs = False      # get_value(kwargs, "score_on_all_pairs", False)

    if at_least_until is None:
        at_least_until = 0

    scorer = ScoringMatrix("identity")          # FIXME: get from sbsp options

    # find the positions of the first two candidate starts in the query
    start = 0
    end = msa_t.alignment_length()

    candidates = list()             # type: List[int]
    for i in range(start, end):

        if msa_t[0][i].isupper():

            if i < region_length or i <= at_least_until:
                candidates.append(i)
            else:
                # compute conservation of block upstream of candidate
                try:
                    conservation = compute_conservation_in_region(
                        msa_t, i - region_length, i,
                        scorer=scorer,
                        direction="upstream",
                        skip_gaps_in_query=True,
                        score_on_all_pairs=score_on_all_pairs,
                    )

                    # if block not conserved, add candidate
                    if conservation < threshold:
                        candidates.append(i)
                    else:
                        break
                except ValueError:
                    candidates.append(i)

    return candidates


def step_a_check_if_at_lorf(candidate_positions):
    # type: (List[int]) -> Union[int, None]
    logger.debug("Func: select-start-positions-from-msa-for-lorf")

    if len(candidate_positions) == 1:
        logger.debug("Only one candidate. Selecting position {}".format(candidate_positions[0]))
        return candidate_positions[0]

    return None




def count_number_of_5prime_candidates_at_position(msa_t, curr_pos, sbsp_options):
    # type: (MSAType, int, SBSPOptions) -> int

    i = curr_pos
    num_upper = 0
    q_curr_type = msa_t[0][i]

    for j in range(msa_t.number_of_sequences()):

        should_count = False

        letter_at_i_j = msa_t[j][i]

        if msa_t[j][i].isupper():

            should_count = True

            t_curr_type = msa_t[j][i]

            if sbsp_options["search-ignore-m-to-l-mutation"]:
                if q_curr_type == "M" and t_curr_type == "L":
                    should_count = False

            if sbsp_options["search-ignore-l-to-m-mutation"]:
                if q_curr_type == "L" and t_curr_type == "M":
                    should_count = False

        # if current position isn't upper, check neighbors
        else:
            if sbsp_options["search-neighbor"]:
                should_count = should_count_in_neighbor(i, msa_t[j].seq._data, sbsp_options, q_curr_type)

        if should_count:
            num_upper += 1

        # penalize
        if sbsp_options.safe_get("search-penalize-standard-aa") is not None:
            if letter_at_i_j in {"v", "l", "i"}:
                num_upper -= sbsp_options.safe_get("search-penalize-standard-aa")

        if sbsp_options.safe_get("search-penalize-no-sequence") is not None:
            if letter_at_i_j == "-":
                if msa_t[j][0:i].count("-") == i:
                    num_upper -= sbsp_options.safe_get("search-penalize-no-sequence")

    return num_upper



def find_first_5prime_that_satisfies_5prime_threshold(msa_t, sbsp_options, begin, radius_aa, direction, skip_gaps_in_query):
    # type: (MSAType, SBSPOptions, int, int, str, bool) -> Union[int, None]

    start_position_in_msa = None
    threshold = 0.5     # FIXME
    num_sequences = msa_t.number_of_sequences()

    curr_pos = begin
    for i in range(radius_aa):
        if curr_pos is None or curr_pos < 0 or curr_pos >= msa_t.alignment_length():
            break

        if not msa_t[0][curr_pos].isupper():
            curr_pos = get_next_position_in_msa(curr_pos, msa_t, direction, skip_gaps_in_query)
            continue

        number_5prime = count_number_of_5prime_candidates_at_position(msa_t, curr_pos, sbsp_options)

        # if threshold
        if number_5prime / float(num_sequences) > threshold:
            start_position_in_msa = curr_pos
            break

        curr_pos = get_next_position_in_msa(curr_pos, msa_t, direction, skip_gaps_in_query)

    return start_position_in_msa

def select_by_upstream_1_4_rule(msa_t, sbsp_options, pos_of_upstream_in_msa):
    # type: (MSAType, SBSPOptions, int) -> Union[int, None]

    start_position_in_msa = None
    radius_aa = 2

    if pos_of_upstream_in_msa is not None and pos_of_upstream_in_msa >= 0:
        # check upstream of position
        start_position_in_msa = find_first_5prime_that_satisfies_5prime_threshold(
            msa_t, sbsp_options, pos_of_upstream_in_msa, radius_aa+1, "upstream", True
        )

        # if not found, try downstream of position
        if start_position_in_msa is None and pos_of_upstream_in_msa < msa_t.alignment_length()-1:
            start_position_in_msa = find_first_5prime_that_satisfies_5prime_threshold(
                msa_t, sbsp_options, pos_of_upstream_in_msa+1, radius_aa, "downstream", True
            )

    return start_position_in_msa


def region_between_two_positions_is_conserved(msa_t, sbsp_options, pos_a, pos_b):
    # type: (MSAType, SBSPOptions, int, int) -> bool
    # compute conservation of block upstream of candidate
    scorer = ScoringMatrix()
    score_on_all_pairs = sbsp_options.safe_get("score-on-all-pairs")
    threshold = 0.5     # FIXME

    conservation = compute_conservation_in_region(
        msa_t, pos_a, pos_b+1,
        scorer=scorer,
        direction="downstream",
        skip_gaps_in_query=False,
        score_on_all_pairs=score_on_all_pairs,
    )

    return conservation > threshold


def candidate_b_has_better_support(msa_t, sbsp_options, pos_a, pos_b, by_at_least=0):
    # type: (MSAType, SBSPOptions, int, int, float) -> bool

    num_5prime_a = count_number_of_5prime_candidates_at_position(msa_t, pos_a, sbsp_options)
    num_5prime_b = count_number_of_5prime_candidates_at_position(msa_t, pos_b, sbsp_options)

    support_a = by_at_least + num_5prime_a / float(msa_t.number_of_sequences())
    support_b = num_5prime_b / float(msa_t.number_of_sequences())

    if isclose(support_a, support_b):

        # run more stringent count
        copy_sbsp_options = copy.deepcopy(sbsp_options)
        copy_sbsp_options["search-neighbor"] = 0

        num_5prime_a = count_number_of_5prime_candidates_at_position(msa_t, pos_a, sbsp_options)
        num_5prime_b = count_number_of_5prime_candidates_at_position(msa_t, pos_b, sbsp_options)

        support_a = by_at_least + num_5prime_a / float(msa_t.number_of_sequences())
        support_b = num_5prime_b / float(msa_t.number_of_sequences())

        if isclose(support_a, support_b):
            return False
        else:
            return support_b > support_a
    # if not close, return better one
    else:
        return support_b > support_a





def select_from_two_neighboring_candidates(msa_t, sbsp_options, pos_a,
                                           pos_b):
    # type: (MSAType, SBSPOptions, int, int) -> int

    if pos_a > pos_b:
        t = pos_a
        pos_a = pos_b
        pos_b = t

    selected = None

    if region_between_two_positions_is_conserved(msa_t, sbsp_options, pos_a, pos_b):
        if candidate_b_has_better_support(msa_t, sbsp_options, pos_a, pos_b):
            selected = pos_b
        else:
            selected = pos_a
    else:

        if sbsp_options.safe_get("search-favor-m"):
            element_a = msa_t[0][pos_a]
            element_b = msa_t[0][pos_b]

            if element_a != "M" and element_b == "M":
                selected = pos_b
            elif element_a == "M" and element_b != "M":
                selected = pos_a

        if selected is None:
            if candidate_b_has_better_support(msa_t, sbsp_options, pos_a, pos_b):
                selected = pos_b
            else:
                selected = pos_a


    return selected



def step_b_find_first_candidate_with_strong_5prime_score(msa_t, candidate_positions, sbsp_options,
                                                         pos_of_upstream_in_msa):
    # type: (MSAType, List[int], SBSPOptions, int) -> Union[int, None]

    threshold = 0.5         # FIXME

    # skip over
    begin = 0
    if pos_of_upstream_in_msa is not None and pos_of_upstream_in_msa >= 0:
        begin = pos_of_upstream_in_msa

    idx_first_valid_candidate = 0
    num_sequences = msa_t.number_of_sequences()
    start_position_in_msa = None

    # skip over all candidates up until 'begin' position
    while idx_first_valid_candidate < len(candidate_positions) and candidate_positions[idx_first_valid_candidate] < begin:
        idx_first_valid_candidate += 1

    for i in range(idx_first_valid_candidate, len(candidate_positions)):

        curr_pos = candidate_positions[i]
        number_5prime = count_number_of_5prime_candidates_at_position(msa_t, curr_pos, sbsp_options)

        # if threshold
        if number_5prime / float(num_sequences) > threshold:
            start_position_in_msa = curr_pos
            break

    # check for nearby downstream
    if start_position_in_msa is not None:
        downstream_start_position_in_msa = find_first_5prime_that_satisfies_5prime_threshold(
            msa_t, sbsp_options, start_position_in_msa+1, sbsp_options["search-better-downstream-aa"],
            "downstream", True
        )

        if downstream_start_position_in_msa is not None:
            start_position_in_msa = select_from_two_neighboring_candidates(
                msa_t, sbsp_options, start_position_in_msa, downstream_start_position_in_msa,
            )

    return start_position_in_msa


def step_c_find_rightmost_by_standard_aa_score(msa_t, candidate_positions, sbsp_options, pos_of_upstream_in_msa):
    # type: (MSAType, List[int], SBSPOptions, int) -> Union[int, None]
    threshold = sbsp_options["search-skip-by-standard-aa-score"]        # type: float

    start_position_in_msa = None

    logger.debug("Func: find-rightmost-by-standard-aa-score")
    for i in reversed(candidate_positions):

        penalized_start_score = sbsp_ml.msa_features_2.compute_simple_saas(msa_t, i)

        logger.debug("Candidate {}, SAAS = {}".format(i, penalized_start_score))

        if penalized_start_score < threshold:
            logger.debug("SAAS < threshold ({}). Select it".format(penalized_start_score))
            start_position_in_msa = i

    logger.debug("No candidate found with low SAAS")

    return start_position_in_msa


def compute_position_of_upstream_in_lorf_nt(series, s):
    # type: (pd.Series, str) -> Union[int, None]

    except_if_not_in_set(s, {"q", "t"})

    # no upstream label
    if series["{}-upstream_left".format(s)] == -1 or series["{}-upstream_right".format(s)] == -1:
        return None

    s_strand = series["{}-strand".format(s)]

    # Compute distance to upstream label (+ means no overlap, - means overlap
    # if the label's strand is "+", then use right of upstream label
    if s_strand == "+":
        distance_of_current_to_upstream = series["{}-left".format(s)] - series["{}-upstream_right".format(s)]
    # otherwise, use left of upstream label (since it's on reverse strand)
    else:
        distance_of_current_to_upstream = series["{}-upstream_left".format(s)] - series["{}-right".format(s)]

    offset_upstream_nt = series["{}-offset".format(s)] - distance_of_current_to_upstream

    return offset_upstream_nt


def convert_ungapped_position_to_gapped(ungapped_position, seq):
    # type: (int, Seq) -> Union[int, None]

    if ungapped_position is None or ungapped_position < 0:
        return None

    seq_length = len(seq)

    curr_pos = 0

    # skip gaps until first none gap
    while curr_pos < seq_length and seq[curr_pos] == "-":
        curr_pos += 1

    for i in range(ungapped_position):

        # state: at non-gap, go over it
        curr_pos += 1

        # if at gap, skip all gaps
        while curr_pos < seq_length and seq[curr_pos] == "-":
            curr_pos += 1

    if curr_pos >= seq_length:
        return None

    return curr_pos

def convert_gapped_position_to_ungapped(gapped_position, seq):
    # type: (int, Seq) -> Union[int, None]

    if gapped_position is None or gapped_position < 0:
        return None

    number_of_gaps_until_position = sum(1 for i in range(gapped_position) if seq[i] == "-")

    return gapped_position - number_of_gaps_until_position


def compute_position_of_upstream_in_msa_for_query(df, msa_t):
    # type: (pd.DataFrame, MSAType) -> Union[int, None]

    pos_of_upstream_in_lorf_nt = compute_position_of_upstream_in_lorf_nt(df.iloc[0], "q")
    if pos_of_upstream_in_lorf_nt is None or pos_of_upstream_in_lorf_nt < 0:
        return None

    pos_of_upstream_in_lorf_aa = int(pos_of_upstream_in_lorf_nt / 3)        # approximate to nearest AA

    return convert_ungapped_position_to_gapped(pos_of_upstream_in_lorf_aa, msa_t[0].seq)


def get_label_from_start_position_in_msa(series, msa_t, start_position_in_msa, s="q"):
    # type: (pd.Series, MSAType, int, str) -> Label

    # convert gapped position to ungapped
    ungapped_offset_aa = convert_gapped_position_to_ungapped(start_position_in_msa, msa_t[0].seq)

    ungapped_offset_nt = ungapped_offset_aa * 3

    s_strand = series["{}-strand".format(s)]

    left = series["{}-left".format(s)]
    right = series["{}-right".format(s)]
    strand = series["{}-strand".format(s)]

    if s_strand == "+":
        left = left + (ungapped_offset_nt - series["{}-offset".format(s)])
    else:
        right = right + (series["{}-offset".format(s)] - ungapped_offset_nt)

    return Label(
        Coordinates(
            left - 1, right - 1, strand
        ),
        seqname=series["{}-accession".format(s)]
    )




def search_for_start_for_msa_and_update_df(df, msa_t, sbsp_options):
    # type: (pd.DataFrame, MSAType, SBSPOptions) -> None
    predicted_at_step = ""

    pos_of_upstream_in_msa = compute_position_of_upstream_in_msa_for_query(df, msa_t)
    at_least_until = None
    if pos_of_upstream_in_msa is not None and pos_of_upstream_in_msa >= 0:
        at_least_until = pos_of_upstream_in_msa + sbsp_options["block-region-length-aa"]

    # get all candidates before conserved block
    candidate_positions = get_all_candidates_before_conserved_block(
        msa_t, sbsp_options,
        at_least_until=at_least_until
    )

    # step A: check if LORF
    start_position_in_msa = step_a_check_if_at_lorf(candidate_positions)

    if start_position_in_msa is None:
        # Step U: Upstream 1,4 rule
        start_position_in_msa = select_by_upstream_1_4_rule(msa_t, sbsp_options, pos_of_upstream_in_msa)

        if start_position_in_msa is None:
            # Step B: find candidate with strong 5' end score
            start_position_in_msa = step_b_find_first_candidate_with_strong_5prime_score(
                msa_t, candidate_positions, sbsp_options, pos_of_upstream_in_msa=pos_of_upstream_in_msa
            )

            if start_position_in_msa is None:
                # Step C:
                start_position_in_msa = step_c_find_rightmost_by_standard_aa_score(
                    msa_t, candidate_positions, sbsp_options, pos_of_upstream_in_msa=pos_of_upstream_in_msa
                )
                if start_position_in_msa is not None:
                    predicted_at_step = "C"
            else:
                predicted_at_step = "B"
        else:
            predicted_at_step = "U"
    else:
        predicted_at_step = "A"

    if start_position_in_msa is not None and start_position_in_msa < 0:
        logger.critical("Somehow, start position {} < 0".format(start_position_in_msa))
        start_position_in_msa = None

    # if all steps failed
    if start_position_in_msa is None:
        df.drop(df.index, inplace=True)
        return  # FIXME: implement recovery strategy



    # get label of new start in genome
    q_label_sbsp = get_label_from_start_position_in_msa(
        df.iloc[0],
        msa_t,
        start_position_in_msa,
        s="q"
    )       # type: Label

    msa_t.add_marker(MSASinglePointMarker(start_position_in_msa, msa_t.alignment_length(), name="selected"))
    msa_t.add_marker(MSASinglePointMarker(pos_of_upstream_in_msa, msa_t.alignment_length(), name="q-3prime", mark="*"))

    df["predicted-at-step"] = predicted_at_step
    df["start-position-in-msa"] = start_position_in_msa
    df["q-left-sbsp"] = q_label_sbsp.left() + 1
    df["q-right-sbsp"] = q_label_sbsp.right() + 1
    df["q-strand-sbsp"] = q_label_sbsp.strand()
    df["msa"] = msa_t


def perform_msa_on_df_with_single_query(env, df, sbsp_options, **kwargs):
    # type: (Environment, pd.DataFrame, SBSPOptions, Dict[str, Any]) -> pd.DataFrame
    """

    Assumption: only a single query exists in this data frame
    :param env:
    :param df:
    :param sbsp_options:
    :param kwargs:
        - inplace: if set to True, df is updated with values, otherwise, a copy is created and returned
    :return: A data frame (either a copy or the input df, based on value of inplace argument)
    """
    inplace = get_value(kwargs, "inplace", False)

    if not inplace:
        df = df.copy()

    if len(df) == 0:
        return df

    # construct msa and filter (if necessary)
    while True:
        msa_t_aa, msa_t_nt = construct_msa_from_df(env, df, sbsp_options, **kwargs)

        # pairwise kimura filter
        targets_before = len(df)

        filter_df_based_on_msa(df, msa_t_aa, msa_t_nt, sbsp_options, inplace=True)

        # if nothing has been filtered or if everything has been filtered, we're done
        if targets_before == len(df) or len(df) == 0:
            break

    if len(df) > 0:
        logger.debug("Searching for start on {} targets".format(len(df)))
        search_for_start_for_msa_and_update_df(df, msa_t_aa, sbsp_options)

    return df

def write_msa_to_directory(df, pd_msa, **kwargs):
    # type: (pd.DataFrame, str, Dict[str, Any]) -> None
    from shutil import copyfile

    msa_number = get_value(kwargs, "msa_number", 0)
    fn_tmp_prefix = get_value(kwargs, "fn_tmp_prefix", 0)

    for msa_t, df_group in df.groupby("msa", as_index=False):
        pf_msa = os.path.join(pd_msa, "msa_{}_{}.txt".format(fn_tmp_prefix, msa_number))

        r = df_group.iloc[0]
        msa_t[0].id = "{};{};{};{}".format(r["q-left"], r["q-right"], r["q-strand"], r["predicted-at-step"])

        # add distance
        for i in range(1, msa_t.number_of_sequences()):
            msa_t[i].id = "{};{}".format(msa_t[i].id, round(df_group.iloc[i-1]["distance"], 4))

        msa_t.to_file(pf_msa)
        df.loc[df_group.index, "pf-msa-output"] = pf_msa

        msa_number += 1


def find_start_for_query_blast_record(env, r, sbsp_options, **kwargs):
    # type: (Environment, Record, SBSPOptions, Dict[str, Any]) -> pd.DataFrame
    """
    Find the 5' end location of the query, and return a data frame with all information, such
    as used orthologs, location of MSA file, etc...
    :param env: Environment
    :param r: Blast record for query
    :param sbsp_options: Options for SBSP
    :param kwargs:
        - pd_msa_final: Path to directory where MSA files are copied (useful for transferring files from
        compute nodes to head node).
    :return: Data frame containing all information. Data frame is empty if not prediction is made.
    """

    msa_number = get_value(kwargs, "msa_number", 0)
    stats = get_value(kwargs, "stats", init=dict)
    fn_tmp_prefix = get_value(kwargs, "msa_output_start", None)
    num_processors = get_value(kwargs, "num_processors", None)

    pd_msa_final = get_value(kwargs, "pd_msa_final", None)

    # read targets and filter out what isn't needed - construct data frame ready for MSA
    df = create_data_frame_for_msa_search_from_blast_results(r, sbsp_options, **kwargs)

    # FIXME  REMOVE
    # if len(df) > 0 and df.iloc[0]["q-left"] == 178270:
    #     pass
    # else:
    #     df.drop(df.index, inplace=True)
    #     return df

    logger.debug("Number of targets after filtering: {}".format(len(df)))

    # run MSA(s) and find gene-start
    perform_msa_on_df_with_single_query(
        env, df, sbsp_options, inplace=True,
        msa_output_start=msa_number,
        msa_number=msa_number, stats=stats,
        fn_tmp_prefix=msa_number,
        num_processors=num_processors
    )

    # for each query in blast
    if pd_msa_final is not None:
        try:
            write_msa_to_directory(df, pd_msa_final, fn_tmp_prefix=fn_tmp_prefix, msa_number=msa_number)
        except Exception:
            pass

    return df


def append_data_frame_to_csv(df, pf_output):
    # type: (pd.DataFrame, str) -> None
    if df is not None and len(df) > 0:
        if not os.path.isfile(pf_output):
            df.to_csv(pf_output, index=False)
        else:
            df.to_csv(pf_output, mode="a", index=False, header=False)


def thread_safe_find_start_and_save_to_csv(env, r, sbsp_options, msa_number, pf_output, lock, **kwargs):
    # type: (Environment, Record, SBSPOptions, int, str, , Dict[str, Any]) -> None
    df_result = find_start_for_query_blast_record(env, r, sbsp_options, msa_number=msa_number, **kwargs)

    append_data_frame_to_csv(df_result, pf_output)


def process_find_start_for_multiple_query_blast_record(lock, process_number, env, records, sbsp_options, pf_output, **kwargs):
    # type: (Lock, int, Environment, List[Record], SBSPOptions, str, Dict[str, Any]) -> None

    msa_number = 0
    for r in records:
        df_result = find_start_for_query_blast_record(env, r, sbsp_options, msa_number="{}_{}".format(process_number, msa_number), **kwargs)

        lock.acquire()
        try:
            append_data_frame_to_csv(df_result, pf_output)
        finally:
            lock.release()

        msa_number += 1


def run_sbsp_steps(env, data, pf_t_db, pf_output, sbsp_options, **kwargs):
    # type: (Environment, Dict[str, Seq], str, str, SBSPOptions, Dict[str, Any]) -> str

    num_processors = get_value(kwargs, "num_processors", None)

    q_sequences = data

    remove_p(pf_output)                     # start clean

    # Run blast
    pf_blast_output = os.path.join(env["pd-work"], "blast_output.xml")
    remove_p(pf_blast_output)
    try:
        run_blast_on_sequences(env, q_sequences, pf_t_db, pf_blast_output, sbsp_options, **kwargs)
        pass
    except ValueError:
        raise ValueError("Couldn't run blast successfully")

    # open blast stream
    try:
        f_blast_results = open(pf_blast_output, "r")
    except OSError:
        raise ValueError("Could not open blast results file: {}".format(pf_blast_output))

    records = NCBIXML.parse(f_blast_results)
    #num_processors = None

    if num_processors is None or num_processors == 0:
        msa_number = 0
        # for each query, find start
        for r in records:
            df_result = find_start_for_query_blast_record(env, r, sbsp_options, msa_number=msa_number, **kwargs)
            append_data_frame_to_csv(df_result, pf_output)
            msa_number += 1
    else:

        logger.debug("Run in parallel mode with {} processors".format(num_processors))

        num_simultaneous_records_per_process = 4
        active_processes = dict()
        worker_id = 0
        kwargs_duplicated = kwargs.copy()
        kwargs_duplicated["num_processors"] = 1

        # run separate process on each split
        lock = Lock()

        while True:

            no_more_records = False

            while len(active_processes) < num_processors:
                list_records = list()           # type: List[Record]

                for r in records:
                    list_records.append(r)
                    if len(list_records) == num_simultaneous_records_per_process:
                        break

                if len(list_records) == 0:
                    no_more_records = True
                    break

                p = Process(target=process_find_start_for_multiple_query_blast_record,
                            args=(lock, worker_id, env, list_records, sbsp_options, pf_output),
                            kwargs={**kwargs_duplicated}
                            )

                logger.debug("Starting process {}".format(worker_id))
                p.start()
                active_processes[worker_id] = p
                worker_id += 1

            if no_more_records:
                break

            # wait until all processes are done
            if len(active_processes) > 0:
                while True:
                    completed_process_ids = set()
                    for i, p in active_processes.items():
                        if not p.is_alive():
                            completed_process_ids.add(i)
                            logger.debug("Done running process {}".format(i))

                    # clean up
                    if len(completed_process_ids) > 0:
                        for i in completed_process_ids:
                            del active_processes[i]

                        break
                    else:
                        time.sleep(1)

        # wait for remaining processes
        for p in active_processes.values():
            p.join()





        # active_processes = set()
        # msa_number = 0
        # parsed_all_records = False
        # manager = Manager()
        # return_dict = manager.dict()
        #
        # # for each query, find start
        # while True:
        #
        #     while not parsed_all_records and len(active_processes) < num_processors:
        #         r = next(records)
        #         if not r:
        #             parsed_all_records = True
        #             break
        #
        #         # create new process
        #         p = Process(target=process_find_start_for_query_blast_record, args=(msa_number, return_dict),
        #                 kwargs={"env": env, "r": r, "sbsp_options": sbsp_options,
        #                         "msa_number": msa_number, **kwargs})
        #         active_processes.add(p)
        #         msa_number += 1
        #
        #     while len(active_processes) > 0:
        #         completed_threads = set()
        #         for p in active_processes:
        #             if not p.is_alive():
        #                 df_result = return_dict[p.]
        #                 completed_threads.add(p)



    return pf_output


def run_sbsp_steps_DEPRECATED(env, data, pf_t_db, pf_output, sbsp_options, **kwargs):
    # type: (Environment, Dict[str, Seq], str, str, SBSPOptions, Dict[str, Any]) -> str

    hsp_criteria = get_value(kwargs, "hsp_criteria", None)
    msa_output_start = get_value(kwargs, "msa_output_start", 0)
    pd_msa_final = get_value(kwargs, "pd_msa_final", env["pd-work"])

    distance_min = sbsp_options.safe_get("distance-min")
    distance_max = sbsp_options.safe_get("distance-max")

    sequences = data

    elapsed_times = dict()
    stats = dict()

    # write sequences to a file
    pf_q_sequences = os.path.join(env["pd-work"], "query_sequences_{}.fasta".format(msa_output_start))
    write_fasta_hash_to_file(sequences, pf_q_sequences)

    # run blast
    curr_time = timeit.default_timer()
    pf_blast_output = os.path.join(env["pd-work"], "blast_output.xml")

    # start clean
    remove_p(pf_blast_output)

    try:
        run_blast_on_sequence_file(env, pf_q_sequences, pf_t_db, pf_blast_output, sbsp_options=sbsp_options)
    except ValueError:
        raise ValueError("Couldn't run blast")

    elapsed_times["1-blastp"] = timeit.default_timer() - curr_time

    # open blast stream
    try:
        f_blast_results = open(pf_blast_output, "r")
    except OSError:
        raise ValueError("Could not open blast results file: {}".format(pf_blast_output))

    # start clean
    remove_p(pf_output)

    matrix = matlist.blosum62
    import sbsp_alg.phylogeny
    sbsp_alg.phylogeny.add_stop_codon_to_blosum(matrix)

    records = NCBIXML.parse(f_blast_results)

    elapsed_times["2-read-filter-per-query"] = 0
    elapsed_times["3-msa-per-query"] = 0
    num_queries = 0
    msa_number = 0
    stats["num-queries-with-support-before-filtering"] = 0
    stats["num-queries-with-support-after-filtering"] = 0
    stats["num-queries-with-support-after-msa"] = 0
    stats["num-queries-with-support-before-pairwise-filtering"] = 0
    stats["num-queries-with-support-after-pairwise-filtering"] = 0
    stats["num-queries-with-support-before-gaps-filtering"] = 0
    stats["num-queries-with-support-after-gaps-filtering"] = 0

    stats_per_query = dict()
    lost_query_at_step = {
            x: 0 for x in ["blast", "kimura", "pairwise-kimura", "msa-gaps", "start-search"]
    }

    # for each blast query
    for r in records:

        query_info = unpack_fasta_header(r.query)

        list_entries = list()
        num_queries += 1
        stats_per_query[num_queries] = 0

        curr_time = timeit.default_timer()
        one_if_one_target = 0

        if len(r.alignments) == 0:
            lost_query_at_step["blast"] += 1

        # for each alignment to a target protein for the current query
        for alignment in r.alignments:

            stats_per_query[num_queries] += 1
            one_if_one_target = 1
            hsp = select_representative_hsp(alignment, hsp_criteria)

            target_info = unpack_fasta_header(alignment.title)

            # get nucleotide sequence that corresponds to protein
            original_q_nt = Seq(query_info["lorf_nt"][query_info["offset"]:])
            original_t_nt = Seq(target_info["lorf_nt"][target_info["offset"]:])

            distance = compute_distance_based_on_local_alignment(query_info, target_info, hsp,
                                                                 original_q_nt=original_q_nt,
                                                                 original_t_nt=original_t_nt,
                                                                 **kwargs)

            original_q_aa = original_q_nt.translate()
            original_t_aa = original_t_nt.translate()

            #global_distance, global_length, global_length_without_gaps = compute_distance_based_on_global_alignment_from_sequences(
            #    original_q_aa, original_t_aa, original_q_nt, original_t_nt, matrix
            #)
            global_distance = global_length = global_length_without_gaps = 0

            # FIXME: thresholds should be from input configuration files
            if distance > distance_min and distance < distance_max:
            #if True:

                output_info = create_info_for_query_target_pair(
                    query_info, target_info, hsp,
                    distance_blast=distance,
                    distance=distance,
                    global_distance=global_distance,
                    global_length=global_length,
                    global_length_without_gaps=global_length_without_gaps,
                    local_distance=distance,
                    local_length=hsp.align_length,
                    local_length_without_gaps=sum([
                        1 for i in range(len(hsp.query)) if hsp.query[i] != "-" and hsp.sbjct[i] != "-"
                    ])
                )

                # clean up
                output_info["q-prot-pos-5prime-in-frag-msa"] = query_info["offset"] / 3
                output_info["q-nucl-pos-5prime-in-frag-msa"] = query_info["offset"]
                output_info["q-prot-position-of-5prime-in-msa-fragment-no-gaps"] = query_info["offset"] / 3
                output_info["q-nucl-position-of-5prime-in-msa-fragment-no-gaps"] = query_info["offset"]
                output_info["t-prot-pos-5prime-in-frag-msa"] = target_info["offset"] / 3
                output_info["t-nucl-pos-5prime-in-frag-msa"] = target_info["offset"]
                output_info["t-prot-position-of-5prime-in-msa-fragment-no-gaps"] = target_info["offset"] / 3
                output_info["t-nucl-position-of-5prime-in-msa-fragment-no-gaps"] = target_info["offset"]

                output_info["q-prot-msa"] = Seq(query_info["lorf_nt"]).translate()._data
                output_info["t-prot-msa"] = Seq(target_info["lorf_nt"]).translate()._data

                output_info["q-nucl-msa"] = Seq(query_info["lorf_nt"])._data
                output_info["t-nucl-msa"] = Seq(target_info["lorf_nt"])._data

                list_entries.append(output_info)


        stats["num-queries-with-support-before-filtering"] += one_if_one_target
        elapsed_times["2-read-filter-per-query"] += timeit.default_timer() - curr_time
        # run MSA on remaining targets

        if len(list_entries) == 0:
            if one_if_one_target != 0:
                lost_query_at_step["kimura"] += 1

            continue

        stats["num-queries-with-support-after-filtering"] += 1

        print("{};{};{};{}".format(query_info["accession"], query_info["left"], query_info["right"], query_info["strand"]))

        curr_time = timeit.default_timer()

        df_entries = pd.DataFrame(list_entries)
        df_results = perform_msa_on_df(env, df_entries, msa_options=sbsp_options, msa_output_start=msa_output_start,
                                       msa_number=msa_number, stats=stats)
        elapsed_times["3-msa-per-query"] += timeit.default_timer() - curr_time
        msa_number += 1


        if stats["num-queries-with-support-before-pairwise-filtering"] > stats["num-queries-with-support-after-pairwise-filtering"]:
            lost_query_at_step["pairwise-kimura"] += 1
        if stats["num-queries-with-support-before-gaps-filtering"] > stats["num-queries-with-support-after-gaps-filtering"]:
            lost_query_at_step["msa-gaps"] += 1

        # for each query in blast
        if pd_msa_final is not None:
            try:
                move_files_using_scp(df_results, pd_msa_final)
            except Exception:
                pass

        # write/append result to file
        if df_results is not None and len(df_results) > 0:
            stats["num-queries-with-support-after-msa"] += 1
            if not os.path.isfile(pf_output):
                df_results.to_csv(pf_output, index=False)
            else:
                df_results.to_csv(pf_output, mode="a", index=False, header=False)

    for key in sorted(elapsed_times.keys()):
        logging.critical("Timer: {},{}".format(key, elapsed_times[key] / 60))

    for key in sorted(stats.keys()):
        logging.critical("Stats: {},{}".format(key, stats[key]))

    for key in sorted(stats_per_query.keys()):
        logging.critical("Stats per query: {},{}".format(key, stats_per_query[key]))


    for key in lost_query_at_step.keys():
        logging.critical("Lost-query: {},{}".format(key, lost_query_at_step[key]))

    return pf_output


def sbsp_steps(env, pipeline_options):
    # type: (Environment, PipelineSBSPOptions) -> Dict[str, Any]
    """
    Given a list of query and target database, run all SBSP steps
    """

    logger.debug("Running: sbsp steps")

    # read input sequences
    q_gil = GenomeInfoList.init_from_file(pipeline_options["pf-q-list"])

    mkdir_p(env["pd-work"])
    pf_aa = os.path.join(env["pd-work"], "query.faa")
    extract_labeled_sequences_for_genomes(env, q_gil, pf_aa,
                                          ignore_frameshifted=True, reverse_complement=True, ignore_partial=True)
    q_sequences = read_fasta_into_hash(pf_aa, stop_at_first_space=False)

    if pipeline_options.use_pbs():
        pbs_options = duplicate_pbs_options_with_updated_paths(env, pipeline_options["pbs-options"], keep_on_head=False)

        if pbs_options.safe_get("pd-data-compute"):
            env = env.duplicate({"pd-data": pbs_options["pd-data-compute"]})

        pbs = PBS(env, pbs_options,
                  splitter=split_dict,
                  merger=merge_identity
                  )

        if pipeline_options.perform_step("build-msa"):

            pd_msa = os.path.join(pbs_options["pd-head"], "msa")
            mkdir_p(pd_msa)

            output = pbs.run(
                data={"dict": q_sequences,
                      "pf_output_template": os.path.join(pbs_options["pd-head"],
                                                         pipeline_options["fn-msa"] + "_{}")},
                func=run_sbsp_steps,
                func_kwargs={
                    "env": env,
                    "pf_t_db": pipeline_options["pf-t-db"],
                    "sbsp_options": pipeline_options["sbsp-options"],
                    "clean": True,
                    "pd_msa_final": pd_msa,
                    "num_processors": pbs_options["num-processors"]
                }
            )

        else:
            # read data from file
            list_pf_output_packages = read_rows_to_list(os.path.join(env["pd-work"], "pbs-summary.txt"))
            output = pbs.merge_output_package_files(list_pf_output_packages)
    else:
        raise NotImplementedError("Only PBS version supported.s")

    return output


def sbsp_step_accuracy(env, pipeline_options, list_pf_previous):
    # type: (Environment, PipelineSBSPOptions, List[str]) -> List[str]
    """
    Given a list of query and target genomes, find the set of related genes
    for each query
    """

    logger.debug("Running: sbsp_step_accuracy")

    mkdir_p(env["pd-work"])

    df = pd.concat([pd.read_csv(f, header=0) for f in list_pf_previous], ignore_index=True)

    df = pipeline_step_compute_accuracy(env, df, pipeline_options)

    df.to_csv(pipeline_options["pf-output"])


    # copy labels
    add_true_starts_to_msa_output(env, df, fn_q_labels_true=pipeline_options["fn-q-labels-true"])
    # add_true_starts_to_msa_output(env, df, msa_nt=True, fn_q_labels_true=pipeline_options["fn-q-labels-true"])
    separate_msa_outputs_by_stats(env, df, pipeline_options["dn-msa-output"])

    return list()




