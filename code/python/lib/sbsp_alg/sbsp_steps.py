import copy
import random
import time
import timeit
from multiprocessing import Process, Lock
from random import shuffle

from Bio.Align import MultipleSeqAlignment
from Bio.Blast import NCBIXML, Record
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from numpy import isclose
from tqdm import tqdm

import sbsp_ml
import sbsp_ml.msa_features
from sbsp_alg.msa import should_count_in_neighbor, filter_by_pairwise_kimura_from_msa
from sbsp_alg.shelf import run_msa_on_sequences
from sbsp_general.labels import Label, Coordinates
from sbsp_container.msa import MSAType, MSASinglePointMarker
from sbsp_general.shelf import append_data_frame_to_csv
from sbsp_io.general import mkdir_p, remove_p
from sbsp_general.general import except_if_not_in_set, os_join
from sbsp_alg.ortholog_finder import extract_labeled_sequences_for_genomes, \
    unpack_fasta_header, select_representative_hsp, create_info_for_query_target_pair, \
    compute_distance_based_on_local_alignment, run_blast_on_sequence_file, is_valid_start
from sbsp_alg.sbsp_compute_accuracy import pipeline_step_compute_accuracy, separate_msa_outputs_by_stats, df_print_labels
from sbsp_general import Environment
from sbsp_io.general import read_rows_to_list
from sbsp_io.msa_2 import add_true_starts_to_msa_output
from sbsp_io.sequences import read_fasta_into_hash, write_fasta_hash_to_file
from sbsp_ml.msa_features import ScoringMatrix
from sbsp_options.parallelization import ParallelizationOptions
from sbsp_options.sbsp import SBSPOptions
from sbsp_options.pipeline_sbsp import PipelineSBSPOptions
from sbsp_parallelization.pbs import PBS
from sbsp_pbs_data.mergers import merge_identity
from sbsp_pbs_data.splitters import *

logger = logging.getLogger(__name__)


def duplicate_parallelization_options_with_updated_paths(env, prl_options, **kwargs):
    # type: (Environment, ParallelizationOptions, Dict[str, Any]) -> ParallelizationOptions
    keep_on_head = get_value(kwargs, "keep_on_head", False, default_if_none=True)

    prl_options = copy.deepcopy(prl_options)
    prl_options["pbs-pd-head"] = os.path.abspath(env["pd-work"])

    if keep_on_head:
        prl_options["pbs-pd-root-compute"] = os.path.abspath(env["pd-work"])
    elif prl_options["pbs-pd-root-compute"] is None:
        prl_options["pbs-pd-root-compute"] = os.path.abspath(env["pd-work"])

    return prl_options


def run_blast_on_sequences(env, q_sequences, pf_t_db, pf_blast_output, sbsp_options, **kwargs):
    # type: (Environment, Dict[str, Seq], str, str, SBSPOptions, Dict[str, Any]) -> None

    fn_tmp_prefix = get_value(kwargs, "fn_tmp_prefix", None)

    # write sequences to a file
    pf_q_sequences = os.path.join(env["pd-work"], "{}query_sequences.fasta".format(fn_tmp_prefix))
    write_fasta_hash_to_file(q_sequences, pf_q_sequences)

    # start clean
    remove_p(pf_blast_output)

    block_size = 0.5

    blast_successful = False
    max_attempts = 3
    attempt = 0

    while not blast_successful and attempt < max_attempts:
        attempt += 1
        try:
            logger.info("Running Diamond Blastp")
            run_blast_on_sequence_file(env, pf_q_sequences, pf_t_db, pf_blast_output, sbsp_options=sbsp_options,
                                       block_size=block_size)
            blast_successful = True
            break
        except ValueError:
            block_size = block_size / 2.0
            logger.info("Blast failed. Trying again with block size {}".format(block_size))

    remove_p(pf_q_sequences)

    if not blast_successful:
        raise ValueError("Couldn't run blast")


def quick_filter_alignments(list_alignments, query_info, **kwargs):
    # type: (List, Dict[str, Any]) -> List
    if len(list_alignments) < 2000:
        return list_alignments

    threshold = 0.6
    threshold_int = int(threshold * 10)
    # binary search your way
    original_q_nt = query_info["lorf_nt"]

    begin = 0
    end = len(list_alignments) - 1

    index_closest = end
    diff_closest = float("inf")

    while (end - begin) > 1:
        mid = int((end + begin) / 2)

        alignment = list_alignments[mid]

        target_info = unpack_fasta_header(alignment.title)
        hsp = select_representative_hsp(alignment, "")  # get reference hit for target

        # get nucleotide sequence that corresponds to proteins

        original_t_nt = target_info["lorf_nt"]

        distance = compute_distance_based_on_local_alignment(query_info, target_info, hsp,
                                                             original_q_nt=original_q_nt,
                                                             original_t_nt=original_t_nt,
                                                             original_q_nt_offset=query_info["offset"],
                                                             original_t_nt_offset=target_info["offset"],
                                                             **kwargs)

        distance_int = int(distance * 10)

        diff = distance_int - threshold_int

        if diff > 0 and diff < diff_closest:
            index_closest = mid
            diff_closest = diff

        if diff == 0:
            index_closest = mid
            break

        if distance_int > threshold_int:
            end = mid
        elif distance_int < threshold_int:
            begin = mid

    return list_alignments[:index_closest]


def create_data_frame_for_msa_search_from_blast_results(r, sbsp_options, **kwargs):
    # type: (Record, SBSPOptions, Dict[str, Any]) -> pd.DataFrame

    df = pd.DataFrame()
    query_info = unpack_fasta_header(r.query)

    # if query_info["left"] != 178270:
    #     logger.debug("{} != 178270".format(query_info["left"]))
    #     return pd.DataFrame()

    distance_min = sbsp_options.safe_get("distance-min")
    distance_max = sbsp_options.safe_get("distance-max")
    rng = get_value(kwargs, "rng", None)
    max_targets = sbsp_options.safe_get("filter-max-number-orthologs")

    filter_orthologs_with_equal_kimura_to_query = sbsp_options.safe_get("filter-orthologs-with-equal-kimura")
    set_target_kimuras = set()

    fsf = False  # filter source found?
    if len(r.alignments) == 0:
        logger.debug("Query Filtered: No blast hits")
        fsf = True

    # for each alignment to a target protein for the current query
    list_entries = list()
    key = "{};{};{};{}".format(query_info["accession"], query_info["left"], query_info["right"], query_info["strand"])
    logger.debug("{}: Reading {} targets from blast".format(key, len(r.alignments)))
    # shuffled_alignments = [a for a in r.alignments]
    before = len(r.alignments)
    shuffled_alignments = quick_filter_alignments(r.alignments, query_info, **kwargs)
    shuffle(shuffled_alignments, random=rng.random)
    logger.debug("Quick filter: {} -> {}".format(before, len(shuffled_alignments)))

    if not fsf and len(shuffled_alignments) == 0:
        logger.debug("Query Filtered: Quick filter by BS")
        fsf = True

    original_q_nt = query_info["lorf_nt"]

    num_analyzed = 0
    acc_lengths = 0

    for alignment in shuffled_alignments:
        if len(list_entries) > max_targets:
            logger.debug("Reached limit on number of targets: {} from {}".format(max_targets, len(shuffled_alignments)))
            break

        target_info = unpack_fasta_header(alignment.title)
        hsp = select_representative_hsp(alignment, "")  # get reference hit for target

        # get nucleotide sequence that corresponds to proteins

        original_t_nt = target_info["lorf_nt"]

        distance = compute_distance_based_on_local_alignment(query_info, target_info, hsp,
                                                             original_q_nt=original_q_nt,
                                                             original_t_nt=original_t_nt,
                                                             original_q_nt_offset=query_info["offset"],
                                                             original_t_nt_offset=target_info["offset"],
                                                             **kwargs)

        acc_lengths += len(original_q_nt)
        num_analyzed += 1
        #        logger.debug("{}, {}".format(round(distance, 2), len(original_q_nt)))

        if distance > distance_min and distance < distance_max:
            if filter_orthologs_with_equal_kimura_to_query is not None:
                rounded = round(distance, filter_orthologs_with_equal_kimura_to_query)
                # if another target with equal (to decimal place) exists, skip this one.
                if rounded in set_target_kimuras:
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
        # sort by distance
        df.sort_values("distance", inplace=True)
        df.reset_index(inplace=True)

    if not fsf and len(list_entries) == 0:
        logger.debug("Query Filtered: Kimura at blast")
        fsf = True

    if num_analyzed > 0:
        logger.debug("Analyzed: {}, {}".format(num_analyzed, round(acc_lengths / float(num_analyzed), 2)))

    return df


def extract_sequences_from_df_for_msa(df):
    # type: (pd.DataFrame) -> List[Seq]

    list_sequences = list()

    if len(df) > 0:

        max_position = min(len(df.iloc[0]["q-lorf_nt"]), df.iloc[0]["q-offset"] + 9000)
        list_sequences.append(Seq(df.iloc[0]["q-lorf_nt"][:max_position]).translate())

        for _, row in df.iterrows():
            max_position = min(len(row["t-lorf_nt"]), row["t-offset"] + 9000)
            list_sequences.append(Seq(row["t-lorf_nt"][:max_position]).translate())

    return list_sequences


def convert_gapped_aa_to_gapped_nt(seq_aa, seq_nt_no_gaps):
    # type: (Seq, str) -> Seq

    seq_nt_with_gaps = ""
    pos_no_gaps = 0
    for i in range(len(seq_aa)):
        if seq_aa[i] == "-":
            seq_nt_with_gaps += "---"
        else:
            seq_nt_with_gaps += seq_nt_no_gaps[pos_no_gaps:pos_no_gaps + 3]
            pos_no_gaps += 3

    return Seq(seq_nt_with_gaps)


def convert_msa_aa_to_nt(msa_t_aa, df):
    # type: (MSAType, pd.DataFrame) -> MSAType

    seq_record_list = list()
    # query sequence
    seq_record_list.append(SeqRecord(convert_gapped_aa_to_gapped_nt(msa_t_aa[0].seq, df.iloc[0]["q-lorf_nt"])))

    # targets
    for i in range(1, msa_t_aa.number_of_sequences()):
        row = df.iloc[i - 1]  # -1 to account for query as first sequence
        seq_record_list.append(SeqRecord(convert_gapped_aa_to_gapped_nt(msa_t_aa[i].seq, row["t-lorf_nt"])))

    return MSAType(MultipleSeqAlignment(seq_record_list))


def lower_case_non_5prime_in_msa(msa_t_aa, msa_t_nt):
    # type: (MSAType, MSAType) -> MSAType

    seq_record_list = list()

    for i in range(msa_t_aa.number_of_sequences()):

        new_seq_aa = ""
        for j_aa in range(msa_t_aa.alignment_length()):
            j_nt = j_aa * 3

            if is_valid_start(msa_t_nt[i].seq._data[j_nt:j_nt + 3], "+"):
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

    fsf = False

    # pairwise Kimura
    if sbsp_options.safe_get("filter-by-pairwise-kimura-from-msa"):
        _, indices_to_keep = filter_by_pairwise_kimura_from_msa(
            [msa_t_nt[i].seq._data for i in range(msa_t_nt.number_of_sequences())], sbsp_options
        )

        indices_in_msa_to_remove = set(range(msa_t_nt.number_of_sequences())).difference(indices_to_keep)
        for i in indices_in_msa_to_remove:
            row_numbers_to_drop.add(i - 1)

        if not fsf and len(row_numbers_to_drop) == len(df):
            logger.debug("Query Filtered: Pairwise Kimura")
            fsf = True

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
                    sequences_that_contribute_to_block.append(j - 1)

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

    if not fsf and len(df) == 0:
        logger.debug("Query Filtered: MSA Gap")
        fsf = True

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


def check_each_column_below_gap_allowance(msa_t, start, end, max_frac_allowed_gaps, **kwargs):
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
    skip_gaps_in_query = get_value(kwargs, "skip_gaps_in_query", False)

    if start < 0 or start >= msa_t.alignment_length():
        raise ValueError("Start of region out of bounds: {} not in [{},{}]".format(start, 0, msa_t.alignment_length()))
    if end < start or end > msa_t.alignment_length():
        raise ValueError(
            "Start of region out of bounds: {} not in [{},{})".format(end, start, msa_t.alignment_length()))

    num_positions_to_analyze = end - start

    curr_pos = start if direction == "downstream" else end - 1
    total_number_of_computations = 0
    pos_score = 0
    num_sequences = msa_t.number_of_sequences()

    for n in range(num_positions_to_analyze):

        if curr_pos is None:
            raise ValueError("Not enough region to compute score")

        num_gaps_in_pos = sum([1 for k in range(msa_t.number_of_sequences()) if msa_t[k][curr_pos] == "-"])

        frac_gaps_in_pos = float(num_gaps_in_pos) / num_sequences

        if frac_gaps_in_pos > max_frac_allowed_gaps:
            return False

        curr_pos = get_next_position_in_msa(curr_pos, msa_t, direction, skip_gaps_in_query)

    return True


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
    score_on_all_pairs = True  # get_value(kwargs, "score_on_all_pairs", False) #FIXME
    skip_gaps_in_query = get_value(kwargs, "skip_gaps_in_query", False)

    if start < 0 or start >= msa_t.alignment_length():
        raise ValueError("Start of region out of bounds: {} not in [{},{}]".format(start, 0, msa_t.alignment_length()))
    if end < start or end > msa_t.alignment_length():
        raise ValueError(
            "Start of region out of bounds: {} not in [{},{})".format(end, start, msa_t.alignment_length()))

    num_positions_to_analyze = end - start

    curr_pos = start if direction == "downstream" else end - 1
    total_number_of_computations = 0
    pos_score = 0
    num_sequences = msa_t.number_of_sequences()

    for n in range(num_positions_to_analyze):

        if curr_pos is None:
            raise ValueError("Not enough region to compute score")

        if score_on_all_pairs:
            for i in range(num_sequences):
                element_i = msa_t[i][curr_pos]
                for j in range(i + 1, num_sequences):
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

    region_length = sbsp_options["block-region-length-aa"]  # get_value(kwargs, "region_length", 10)
    threshold = 0.5  # FIXME
    score_on_all_pairs = False  # get_value(kwargs, "score_on_all_pairs", False)

    if at_least_until is None:
        at_least_until = 0

    scorer = ScoringMatrix("identity")  # FIXME: get from sbsp options

    # find the positions of the first two candidate starts in the query
    start = 0
    end = msa_t.alignment_length() - region_length

    candidates = list()  # type: List[int]
    passed_column_of_no_gaps = False
    for i in range(start, end):

        if not passed_column_of_no_gaps:
            if sum([1 for k in range(msa_t.number_of_sequences()) if
                    msa_t[k][i] != "-"]) == msa_t.number_of_sequences():
                passed_column_of_no_gaps = True

        if msa_t[0][i].isupper():
            candidates.append(i)
        elif i >= at_least_until and len(candidates) > 0 and passed_column_of_no_gaps:
            # compute conservation of block upstream of candidate
            try:
                conservation = compute_conservation_in_region(
                    msa_t, i, i + region_length,
                    scorer=scorer,
                    direction="downstream",
                    skip_gaps_in_query=True,
                    score_on_all_pairs=score_on_all_pairs,
                )

                columns_saturated = check_each_column_below_gap_allowance(msa_t, i, i + region_length, 0.2,
                                                                          direction="downstream",
                                                                          skip_gaps_in_query=True
                                                                          )
                first_column_conservation = compute_conservation_in_region(
                    msa_t, i, i + 1,
                    scorer=scorer,
                    direction="downstream",
                    skip_gaps_in_query=False,
                    score_on_all_pairs=score_on_all_pairs,
                )

                # if block not conserved, add candidate
                if conservation > threshold and columns_saturated and first_column_conservation > threshold:
                    break
            except ValueError:
                pass

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


def find_first_5prime_that_satisfies_5prime_threshold(msa_t, sbsp_options, begin, radius_aa, direction,
                                                      skip_gaps_in_query):
    # type: (MSAType, SBSPOptions, int, int, str, bool) -> Union[int, None]

    start_position_in_msa = None
    threshold = 0.5  # FIXME
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
    radius_aa = 4

    if pos_of_upstream_in_msa is not None and pos_of_upstream_in_msa >= 0:
        # check upstream of position
        start_position_in_msa = find_first_5prime_that_satisfies_5prime_threshold(
            msa_t, sbsp_options, pos_of_upstream_in_msa, radius_aa + 1, "upstream", True
        )

        # if not found, try downstream of position
        if start_position_in_msa is None and pos_of_upstream_in_msa < msa_t.alignment_length() - 1:
            start_position_in_msa = find_first_5prime_that_satisfies_5prime_threshold(
                msa_t, sbsp_options, pos_of_upstream_in_msa + 1, radius_aa, "downstream", True
            )

    return start_position_in_msa


def region_between_two_positions_is_conserved(msa_t, sbsp_options, pos_a, pos_b, **kwargs):
    # type: (MSAType, SBSPOptions, int, int, Dict[str, Any]) -> bool
    # compute conservation of block upstream of candidate
    scorer = ScoringMatrix()
    score_on_all_pairs = sbsp_options.safe_get("score-on-all-pairs")
    threshold = get_value(kwargs, "threshold", 0.5, default_if_none=True)  # FIXME

    conservation = compute_conservation_in_region(
        msa_t, pos_a, pos_b + 1,
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

    support_a_raw = num_5prime_a / float(msa_t.number_of_sequences())
    support_a = by_at_least + num_5prime_a / float(msa_t.number_of_sequences())
    support_b = num_5prime_b / float(msa_t.number_of_sequences())

    if isclose(support_a_raw, support_b):

        # run more stringent count
        copy_sbsp_options = copy.deepcopy(sbsp_options)
        copy_sbsp_options["search-neighbor"] = 0
        copy_sbsp_options["search-limit-gap-skips"] = 0

        num_5prime_a = count_number_of_5prime_candidates_at_position(msa_t, pos_a, copy_sbsp_options)
        num_5prime_b = count_number_of_5prime_candidates_at_position(msa_t, pos_b, copy_sbsp_options)

        support_a_raw = num_5prime_a / float(msa_t.number_of_sequences())
        support_a = by_at_least + num_5prime_a / float(msa_t.number_of_sequences())
        support_b = num_5prime_b / float(msa_t.number_of_sequences())

        if isclose(support_a_raw, support_b):
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

    if not region_between_two_positions_is_conserved(msa_t, sbsp_options, pos_a, pos_b, threshold=0.3):
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
            if candidate_b_has_better_support(msa_t, sbsp_options, pos_a, pos_b, by_at_least=0.3):
                selected = pos_b
            else:
                selected = pos_a

    return selected


def step_b_find_first_candidate_with_strong_5prime_score(msa_t, candidate_positions, sbsp_options,
                                                         pos_of_upstream_in_msa):
    # type: (MSAType, List[int], SBSPOptions, int) -> Union[int, None]

    threshold = 0.5  # FIXME

    # skip over
    begin = 0
    if pos_of_upstream_in_msa is not None and pos_of_upstream_in_msa >= 0:
        begin = pos_of_upstream_in_msa

    idx_first_valid_candidate = 0
    num_sequences = msa_t.number_of_sequences()
    start_position_in_msa = None

    # skip over all candidates up until 'begin' position
    while idx_first_valid_candidate < len(candidate_positions) and candidate_positions[
        idx_first_valid_candidate] < begin:
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
        region_begin = start_position_in_msa + 1
        region_end = min(region_begin + sbsp_options["search-better-downstream-aa"], msa_t.alignment_length())
        region_len = region_end - region_begin
        remaining_region = region_len
        cursor = region_begin

        pos_a = curr_pos
        pos_b = find_first_5prime_that_satisfies_5prime_threshold(
            msa_t, sbsp_options, cursor, remaining_region,
            "downstream", True
        )

        if pos_b is not None:
            cursor = pos_b + 1

        while pos_b is not None:
            new_pos = select_from_two_neighboring_candidates(
                msa_t, sbsp_options, start_position_in_msa, pos_b,
            )
            start_position_in_msa = new_pos

            pos_a = new_pos
            remaining_region = remaining_region - (pos_b - pos_a + 1)
            if remaining_region <= 0:
                break

            pos_b = find_first_5prime_that_satisfies_5prime_threshold(
                msa_t, sbsp_options, cursor, remaining_region,
                "downstream", True
            )
            if pos_b is not None:
                cursor = pos_b + 1
        # downstream_start_position_in_msa = find_first_5prime_that_satisfies_5prime_threshold(
        #     msa_t, sbsp_options, start_position_in_msa+1, sbsp_options["search-better-downstream-aa"],
        #     "downstream", True
        # )

        # pos_a = curr_pos
        # pos_b = downstream_start_position_in_msa

        # if downstream_start_position_in_msa is not None:
        #     start_position_in_msa = select_from_two_neighboring_candidates(
        #         msa_t, sbsp_options, start_position_in_msa, downstream_start_position_in_msa,
        #     )

    return start_position_in_msa


def step_c_find_rightmost_by_standard_aa_score(msa_t, candidate_positions, sbsp_options, pos_of_upstream_in_msa):
    # type: (MSAType, List[int], SBSPOptions, int) -> Union[int, None]
    threshold = sbsp_options["search-skip-by-standard-aa-score"]  # type: float

    start_position_in_msa = None

    logger.debug("Func: find-rightmost-by-standard-aa-score")
    for i in reversed(candidate_positions):

        penalized_start_score = sbsp_ml.msa_features.compute_simple_saas(msa_t, i)

        logger.debug("Candidate {}, SAAS = {}".format(i, penalized_start_score))

        if penalized_start_score < threshold:
            logger.debug("SAAS < threshold ({}). Select it".format(penalized_start_score))
            start_position_in_msa = i
            break

    logger.debug("No candidate found with low SAAS")

    return start_position_in_msa


def compute_position_of_upstream_in_lorf_nt(series, s):
    # type: (pd.Series, str) -> Union[int, None]

    except_if_not_in_set(s, {"q", "t"})

    # no upstream label
    if series["{}-upstream_left".format(s)] == -1 or series["{}-upstream_right".format(s)] == -1:
        return None

    # only consider on the same strand
    if series["{}-upstream_strand".format(s)] != series["{}-strand".format(s)]:
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

    pos_of_upstream_in_lorf_aa = int(pos_of_upstream_in_lorf_nt / 3)  # approximate to nearest AA

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
        # Step B: Upstream 1,4 rule
        start_position_in_msa = select_by_upstream_1_4_rule(msa_t, sbsp_options, pos_of_upstream_in_msa)

        if start_position_in_msa is None:
            # Step C: find candidate with strong 5' end score
            start_position_in_msa = step_b_find_first_candidate_with_strong_5prime_score(
                msa_t, candidate_positions, sbsp_options, pos_of_upstream_in_msa=pos_of_upstream_in_msa
            )

            if start_position_in_msa is None:
                # Step C:
                # start_position_in_msa = step_c_find_rightmost_by_standard_aa_score(
                #    msa_t, candidate_positions, sbsp_options, pos_of_upstream_in_msa=pos_of_upstream_in_msa
                # #)
                # copy_sbsp_options = copy.deepcopy(sbsp_options)
                # copy_sbsp_options["search-neighbor"] = 4
                # copy_sbsp_options["search-limit-gap-skips"] = 2
                #
                # start_position_in_msa = step_b_find_first_candidate_with_strong_5prime_score(
                #     msa_t, candidate_positions, copy_sbsp_options, pos_of_upstream_in_msa=pos_of_upstream_in_msa
                # )
                # if start_position_in_msa is not None:
                #     predicted_at_step = "D"
                pass
            else:
                predicted_at_step = "C"
        else:
            predicted_at_step = "B"
    else:
        predicted_at_step = "A"

    if start_position_in_msa is not None and start_position_in_msa < 0:
        logger.critical("Somehow, start position {} < 0".format(start_position_in_msa))
        start_position_in_msa = None

    # if all steps failed
    if start_position_in_msa is None:
        logger.debug("Search could not find start. Support: {}".format(len(df)))
        df.drop(df.index, inplace=True)
        return  # FIXME: implement recovery strategy

    logger.debug("Step {}: S5 = {}".format(predicted_at_step,
                                          count_number_of_5prime_candidates_at_position(msa_t, start_position_in_msa,
                                                                                        sbsp_options) / float(
                                              msa_t.number_of_sequences())
                                          ))

    # get label of new start in genome
    q_label_sbsp = get_label_from_start_position_in_msa(
        df.iloc[0],
        msa_t,
        start_position_in_msa,
        s="q"
    )  # type: Label

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

    qkey = df.iloc[0]["q-key"] if "q-key" in df else None

    fsf = False

    # construct msa and filter (if necessary)
    while True:
        curr_time = timeit.default_timer()
        msa_t_aa, msa_t_nt = construct_msa_from_df(env, df, sbsp_options, **kwargs)
        logger.debug("MSA: Time (min): {}, Support: {}, Key: {}".format(
            round((timeit.default_timer() - curr_time) / 60.0, 2), len(df), qkey
        ))

        # pairwise kimura filter
        targets_before = len(df)

        curr_time = timeit.default_timer()
        filter_df_based_on_msa(df, msa_t_aa, msa_t_nt, sbsp_options, inplace=True)
        logger.debug("Filter: Time (min): {}, Support: {}, Key: {}".format(
            round((timeit.default_timer() - curr_time) / 60.0, 2), len(df), qkey
        ))

        # if nothing has been filtered or if everything has been filtered, we're done
        if targets_before == len(df) or len(df) == 0:
            break

    if len(df) > 0:
        logger.debug("Searching for start on {} targets".format(len(df)))
        curr_time = timeit.default_timer()
        search_for_start_for_msa_and_update_df(df, msa_t_aa, sbsp_options)
        logger.debug("Search: Time (min): {}, Support: {}, Key: {}".format(
            round((timeit.default_timer() - curr_time) / 60.0, 2), len(df), qkey
        ))

        if not fsf and len(df) == 0:
            logger.debug("Query Filtered: Start search")
            fsf = True

    return df


def write_msa_to_directory(df, pd_msa, **kwargs):
    # type: (pd.DataFrame, str, Dict[str, Any]) -> None

    msa_number = get_value(kwargs, "msa_number", 0)
    fn_tmp_prefix = get_value(kwargs, "fn_tmp_prefix", 0)

    for msa_t, df_group in df.groupby("msa", as_index=False):
        pf_msa = os.path.join(pd_msa, "msa_{}_{}.txt".format(fn_tmp_prefix, msa_number))

        r = df_group.iloc[0]
        msa_t[0].id = "{};{};{};{}".format(r["q-left"], r["q-right"], r["q-strand"], r["predicted-at-step"])

        # add distance
        for i in range(1, msa_t.number_of_sequences()):
            msa_t[i].id = "{};{}".format(msa_t[i].id, round(df_group.iloc[i - 1]["distance"], 4))

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
    rng = get_value(kwargs, "rng", None)

    pd_msa_final = get_value(kwargs, "pd_msa_final", None)

    # read targets and filter out what isn't needed - construct data frame ready for MSA
    curr_time = timeit.default_timer()
    df = create_data_frame_for_msa_search_from_blast_results(r, sbsp_options, **kwargs)
    elapsed_time = round((timeit.default_timer() - curr_time) / 60.0, 2)
    # logger.info("CDFFMSFBR ({}): {} (min) for {} orthologs".format(msa_number, elapsed_time, len(r.alignments)))

    # FIXME  REMOVE
    # if len(df) > 0 and df.iloc[0]["q-left"] == 178270:
    #     pass
    # else:
    #     df.drop(df.index, inplace=True)
    #     return df

    # logger.debug("Number of targets after filtering: {}".format(len(df)))

    # run MSA(s) and find gene-start
    curr_time = timeit.default_timer()
    perform_msa_on_df_with_single_query(
        env, df, sbsp_options, inplace=True,
        msa_output_start=msa_number,
        msa_number=msa_number, stats=stats,
        fn_tmp_prefix=msa_number,
        num_processors=num_processors
    )
    elapsed_time = round((timeit.default_timer() - curr_time) / 60.0, 2)
    # logger.info("PMODWSQ ({}): {} (min) for {} orthologs".format(msa_number, elapsed_time, len(df)))

    # for each query in blast
    if pd_msa_final is not None:
        try:
            write_msa_to_directory(df, pd_msa_final, fn_tmp_prefix=fn_tmp_prefix, msa_number=msa_number)
        except Exception:
            pass

    if "q-lorf_nt" in df.columns and "t-lorf_nt" in df.columns:
        df.drop(["q-lorf_nt", "t-lorf_nt"], axis=1, inplace=True)
    return df


def process_find_start_for_multiple_query_blast_record(lock, process_number, env, records, sbsp_options, pf_output,
                                                       **kwargs):
    # type: (Lock, int, Environment, List[Record], SBSPOptions, str, Dict[str, Any]) -> None

    # local_rng = np.random.RandomState(sbsp_options.safe_get("random-seed"))
    local_rng = random.Random(sbsp_options.safe_get("random-seed"))

    msa_number = 0
    for r in records:
        df_result = find_start_for_query_blast_record(
            env, r, sbsp_options, msa_number="{}_{}".format(process_number, msa_number), rng=local_rng, **kwargs
        )

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
    # REMOVE
    # q_sequences = debug_filter_queries(q_sequences)
    # num_processors = None

    remove_p(pf_output)  # start clean

    # Run blast
    pf_blast_output = os.path.join(env["pd-work"], "blast_output.xml")
    remove_p(pf_blast_output)
    try:
        curr_time = timeit.default_timer()
        run_blast_on_sequences(env, q_sequences, pf_t_db, pf_blast_output, sbsp_options, **kwargs)
        logger.info("Blast runtime (min): {:.2f}".format((timeit.default_timer() - curr_time) / float(60)))
    except ValueError:
        remove_p(pf_blast_output)
        raise ValueError("Couldn't run blast successfully")

    # open blast stream
    try:
        f_blast_results = open(pf_blast_output, "r")
    except OSError:
        raise ValueError("Could not open blast results file: {}".format(pf_blast_output))

    records = NCBIXML.parse(f_blast_results)

    # REMOVE

    if num_processors is None or num_processors == 0:
        msa_number = 0
        # for each query, find start
        for r in tqdm(records, total=len(records)):
            # REMOVE
            # query_info = unpack_fasta_header(r.query)
            # if  int(query_info["right"]) not in {449870}:
            #    logger.debug("Skipping: {}".format(query_info["right"]))
            #    continue

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
                list_records = list()  # type: List[Record]

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

    remove_p(pf_blast_output)


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
    fn_labels = pipeline_options["fn-q-labels"]
    extract_labeled_sequences_for_genomes(env, q_gil, pf_aa,
                                          ignore_frameshifted=True, reverse_complement=True, ignore_partial=True,
                                          fn_labels=fn_labels)
    q_sequences = read_fasta_into_hash(pf_aa, stop_at_first_space=False)

    if pipeline_options.use_pbs():
        prl_options = duplicate_parallelization_options_with_updated_paths(env, pipeline_options["prl-options"],
                                                                           keep_on_head=False)

        if prl_options.safe_get("pd-data-compute"):
            env = env.duplicate({"pd-data": prl_options["pd-data-compute"]})

        pbs = PBS(env, prl_options,
                  splitter=split_dict,
                  merger=merge_identity
                  )

        if pipeline_options.perform_step("prediction"):

            pd_msa = os.path.join(prl_options["pbs-pd-head"], "msa")
            mkdir_p(pd_msa)

            output = pbs.run(
                data={"dict": q_sequences,
                      "pf_output_template": os.path.join(prl_options["pbs-pd-head"],
                                                         pipeline_options["fn-msa"] + "_{}")},
                func=run_sbsp_steps,
                func_kwargs={
                    "env": env,
                    "pf_t_db": pipeline_options["pf-t-db"],
                    "sbsp_options": pipeline_options["sbsp-options"],
                    "clean": True,
                    "pd_msa_final": pd_msa,
                    "num_processors": prl_options["pbs-ppn"]
                }
            )

        else:
            # read data from file
            list_pf_output_packages = read_rows_to_list(os.path.join(env["pd-work"], "pbs-summary.txt"))
            output = pbs.merge_output_package_files(list_pf_output_packages)
    else:
        pd_msa = os_join(env["pd-work"], "msa")
        mkdir_p(pd_msa)

        if pipeline_options.perform_step("prediction"):
            run_sbsp_steps(
                env, q_sequences,
                pipeline_options["pf-t-db"], pipeline_options["pf-output"], pipeline_options["sbsp-options"],
                clean=True,
                pd_msa_final=pd_msa,
                num_processors=pipeline_options["prl-options"]["num-processors"],
            )

        output = [pipeline_options["pf-output"]]

    return output


def sbsp_step_compare(env, pipeline_options, list_pf_previous):
    # type: (Environment, PipelineSBSPOptions, List[str]) -> List[str]
    """
    Given a list of query and target genomes, find the set of related genes
    for each query
    """

    mkdir_p(env["pd-work"])


    if len(list_pf_previous) == 0:
        raise ValueError("Cannot produce results: {}".format(pipeline_options["pf-q-list"]))

    # read data
    df = pd.concat([pd.read_csv(f, header=0) for f in list_pf_previous], ignore_index=True)

    # get labels
    df_print_labels(env, df, "q", suffix_coordinates="sbsp",
                                          suffix_fname="")

    if pipeline_options.perform_step("comparison"):
        df = pipeline_step_compute_accuracy(env, df, pipeline_options)

        df.to_csv(pipeline_options["pf-output"])

        # copy labels
        add_true_starts_to_msa_output(env, df, fn_q_labels_true=pipeline_options["fn-q-labels-compare"])
        # add_true_starts_to_msa_output(env, df, msa_nt=True, fn_q_labels_true=pipeline_options["fn-q-labels-true"])
        separate_msa_outputs_by_stats(env, df, pipeline_options["dn-msa-output"])
    return list()
