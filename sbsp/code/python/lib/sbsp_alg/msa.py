from __future__ import print_function

import os
import Bio.Align
import sbsp_general.dataframe
import numpy as np
from sbsp_general.general import get_value
import sbsp_general.labels
from sbsp_options.msa import MSAOptions
from typing import *
import pandas as pd
import copy
import sbsp_ml.msa_features_2
from sbsp_general.msa_2 import MSAType

import logging
logger = logging.getLogger(__name__)



def msa(sequences, **kwargs):
    # type: (list[str], Dict[str, Any]) -> [list[str], str, str]

    tags = get_value(kwargs, "tags", None)
    pd_work = get_value(kwargs, "pd_work", ".")
    suffix_fname = get_value(kwargs, "suffix_fname", "")

    if suffix_fname != "":
        suffix_fname = "-" + str(suffix_fname)


    if len(sequences) == 0:
        return None

    if tags is not None:
        if len(tags) != len(sequences):
            raise ValueError("Length of tags list ({}) != number of sequences ({})".format(len(tags), len(sequences)))

    # convert list of sequences to string (with fasta)
    data = ""
    for i in range(len(sequences)):
        if i > 0:
            data += "\n"

        curr_tag = ""
        if tags:
            curr_tag = "-" + str(tags[i])
        data += ">{}{}\n{}".format(i, curr_tag, sequences[i].replace("J", ""))

    # enter into muscle
    # from Bio.Align.Applications import MuscleCommandline
    # muscle_cline = MuscleCommandline(clwstrict=True)           # stable keeps input order in output
    # stdout, stderr = muscle_cline(stdin=data)

    from Bio import AlignIO
    # from StringIO import StringIO
    # align = AlignIO.read(StringIO(stdout), "clustal")

    # enter into clustal
    from Bio.Align.Applications import ClustalwCommandline
    import sbsp_io.general
    pf_tmp = os.path.join(pd_work, "tmp{}.msa".format(suffix_fname))
    pf_aln = os.path.join(pd_work, "tmp{}.aln".format(suffix_fname))

    sbsp_io.general.write_string_to_file(data,pf_tmp )
    clustalw_cline = ClustalwCommandline("clustalw2", infile=pf_tmp, gapopen=20)
    stdout, stderr = clustalw_cline()


    align = AlignIO.read(pf_aln, "clustal")

    list_alignments = [a for a in align]
    list_alignments.sort(key=lambda r: float(r.id.strip().split("-")[0]))

    align._records = list_alignments
    stdout = align.format("clustal")

    output_list = [a.seq._data for a in list_alignments]
    return [output_list, stdout, stderr, align]



# TODO
# convert AA alignment to nucleotide alignment
#    account for stops
#

def add_stops_to_aa_alignment(aa_alignment_no_stops, aa_sequence_with_stops):
    # type: (str, str) -> str

    pos_in_aa_sequence = 0

    aa_alignment_with_stops = ""

    for pos_in_aa_alignment in range(len(aa_alignment_no_stops)):

        c_align = aa_alignment_no_stops[pos_in_aa_alignment]

        # gaps
        if c_align == "-":
            aa_alignment_with_stops += "-"          # add gap
            continue

        c_seq = aa_sequence_with_stops[pos_in_aa_sequence]

        while c_seq == "*":
            aa_alignment_with_stops += "*"
            pos_in_aa_sequence += 1
            c_seq = aa_sequence_with_stops[pos_in_aa_sequence]

        if c_seq.upper() != c_align.upper():
            raise ValueError("Sequences don't match...")

        aa_alignment_with_stops += c_seq

        pos_in_aa_sequence += 1

    return aa_alignment_with_stops



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


def convert_using_marked(aa_alignment, marks):
    # type: (str, str) -> str
    import string

    result = ""

    for i in range(len(marks)):

        # if reached start
        if marks[i] != " ":
            if aa_alignment[i] == "-":
                raise ValueError("Mismatch between marks and alignment")
            else:
                result += aa_alignment[i].upper()
        # if not start, lower it
        else:
            if aa_alignment[i] == "-":
                result += "-"
            else:
                result += aa_alignment[i].lower()

    return result

def mark_all_candidate_starts(seq_aa_with_gaps, seq_nt_no_gaps):

    seq_candidates = ""

    pos_in_no_gap_aa = 0

    for i in range(len(seq_aa_with_gaps)):

        if seq_aa_with_gaps[i] == "-":
            seq_candidates += " "
            continue

        pos_in_no_gap_nt = pos_in_no_gap_aa * 3

        codon = seq_nt_no_gaps[pos_in_no_gap_nt:pos_in_no_gap_nt+3]

        if codon == "ATG":
            seq_candidates += "A"
        elif codon == "GTG":
            seq_candidates += "G"
        elif codon == "TTG":
            seq_candidates += "T"
        else:
            seq_candidates += " "

        pos_in_no_gap_aa += 1

    return seq_candidates




def lower_case_everything_except_starts(alignments_aa, list_sequences_aa, list_sequences_nt):
    # type: (Bio.Align.MultipleSeqAlignment, list) -> Bio.Align.MultipleSeqAlignment

    i = 0

    for alignment in alignments_aa:

        aa_align = alignment.seq._data
        nt_seq = list_sequences_nt[i]
        aa_seq = list_sequences_aa[i]       # type: str

        aa_seq = aa_seq.replace("J", "")

        # add stops to alignment
        try:
            aa_align_with_stops = add_stops_to_aa_alignment(aa_align, aa_seq)
        except:
            pass

        # mark candidate starts
        marked_alignment = mark_all_candidate_starts(aa_align_with_stops, nt_seq)

        # convert aa alignment to lowercase except marks
        aa_align_with_stops_lowered = convert_using_marked(aa_align_with_stops, marked_alignment)

        # remove all stops again
        aa_align_no_stops_lowered = aa_align_with_stops_lowered.replace("*", "")

        alignment.seq._data = aa_align_no_stops_lowered

        i += 1

    return alignments_aa


# deprecated
def lower_case_everything_except_starts_deprecated(list_alignments_aa, list_sequences_nt):
    # type: (list, list) -> list

    list_result = list()

    for i in range(len(list_alignments_aa)):

        aa_align = list_alignments_aa[i]
        nt_seq = list_sequences_nt[i]

        # add stops to alignment
        aa_align_with_stops = add_stops_to_aa_alignment(aa_align, nt_seq)

        # mark candidate starts
        # from sbsp_viz.phylogeny import mark_all_candidate_starts
        marked_alignment = mark_all_candidate_starts(aa_align_with_stops, nt_seq)

        # convert aa alignment to lowercase except marks
        aa_align_with_stops_lowered = convert_using_marked(aa_align_with_stops, marked_alignment)

        # remove all stops again
        aa_align_no_stops_lowered = aa_align_with_stops_lowered.replace("*", "")

        list_result.append(aa_align_no_stops_lowered)

    return list_result


def convert_alignments_to_str(alignments):
    return alignments.format("clustal")


def should_count_in_neighbor(pos, alignment, msa_options, q_curr_type):
    # type: (int, str, MSAOptions, str) -> bool

    limit_gap_skips = msa_options.safe_get("search-limit-gap-skips")

    if msa_options["search-neighbor"] > 0:
        # upstream

        upstream_should_count = False
        downstream_should_count = False

        for r in range(1, msa_options["search-neighbor"]+1):
            r_upstr = r
            gap_skips = 0
            while pos - r_upstr >= 0 and alignment[pos - r_upstr] == "-":       # skip over all gaps
                if limit_gap_skips is not None:
                    if gap_skips >= limit_gap_skips:
                        break
                r_upstr += 1
                gap_skips += 1

            if pos - r_upstr >= 0:
                t_curr_type = alignment[pos - r_upstr]

                if t_curr_type.isupper():

                    upstream_should_count = True

                    if msa_options["search-ignore-m-to-l-mutation"]:
                        if q_curr_type == "M" and t_curr_type == "L":
                            upstream_should_count = False

            r_dnstr = r
            # skip gaps
            gap_skips = 0
            while pos + r_dnstr < len(alignment) and alignment[pos + r_dnstr] == "-":
                if limit_gap_skips is not None:
                    if gap_skips >= limit_gap_skips:
                        break
                r_dnstr += 1
                gap_skips += 1

            if pos + r_dnstr < len(alignment):
                t_curr_type = alignment[pos + r_dnstr]
                if t_curr_type.isupper():
                    downstream_should_count = True

                    if msa_options["search-ignore-m-to-l-mutation"]:
                        if q_curr_type == "M" and t_curr_type == "L":
                            downstream_should_count = False

            if upstream_should_count or downstream_should_count:
                return True

    return False

def count_num_upper(list_alignments, i, msa_options):
    num_upper = 0
    q_curr_type = list_alignments[0][i]

    for j in range(len(list_alignments)):

        should_count = False

        letter_at_i_j = list_alignments[j][i]

        if list_alignments[j][i].isupper():

            should_count = True

            t_curr_type = list_alignments[j][i]

            if msa_options["search-ignore-m-to-l-mutation"]:
                if q_curr_type == "M" and t_curr_type == "L":
                    should_count = False

            if msa_options["search-ignore-l-to-m-mutation"]:
                if q_curr_type == "L" and t_curr_type == "M":
                    should_count = False

        # if current position isn't upper, check neighbors
        else:
            if msa_options["search-neighbor"]:
                should_count = should_count_in_neighbor(i, list_alignments[j], msa_options, q_curr_type)

        if should_count:
            num_upper += 1

        # penalize
        if msa_options.safe_get("search-penalize-standard-aa") is not None:
            if letter_at_i_j in {"v", "l", "i"}:
                num_upper -= msa_options.safe_get("search-penalize-standard-aa")

        if msa_options.safe_get("search-penalize-no-sequence") is not None:
            if letter_at_i_j == "-":
                if list_alignments[j][0:i].count("-") == i:
                    num_upper -= msa_options.safe_get("search-penalize-no-sequence")

    return num_upper

def get_aligned_letters_at_position(list_alignments, i, upstream_region=0, downstream_region=0):

    list_items = list()

    for j in range(len(list_alignments)):

        should_count = False

        t_curr_type = list_alignments[j][i]

        if not t_curr_type.isupper():
            # search neighborhood
            for radius in range(1, max(upstream_region, downstream_region)+1):

                # check upstream
                if radius <= upstream_region:
                    new_pos = i - radius
                    if new_pos >= 0 and list_alignments[j][new_pos].isupper():
                        t_curr_type = list_alignments[j][new_pos]
                        break

                if radius <= downstream_region:
                    new_pos = i + radius
                    if new_pos < len(list_alignments[0]) and list_alignments[j][new_pos].isupper():
                        t_curr_type = list_alignments[j][new_pos]
                        break

        list_items.append(t_curr_type)

    return list_items



def compute_conservation_in_region(list_alignments, begin, end, direction="upstream", skip_gaps=False, only_full_length=False,
                                   max_fraction_of_gaps_in_pos=None, scorer=None, **kwargs):
    # type: (List[str], int, int, str, bool, bool, bool, ScoringFunction) -> float

    """

    :param alignments:
    :param begin: (inclusive)
    :param end: (exclusive)
    :return:
    """

    score_on_all_pairs = get_value(kwargs, "score_on_all_pairs", False)

    def gap_fraction_in_pos(list_alignments, pos):
        total = len(list_alignments)
        if total == 0:
            return 0.0

        return sum(1 for i in range(len(list_alignments)) if list_alignments[i][pos] == "-") / float(total)

    number_of_positions = end - begin

    total = 0
    matched = 0

    if max_fraction_of_gaps_in_pos is None:
        max_fraction_of_gaps_in_pos = 1

    skipped_gaps = 0

    for i in range(int(number_of_positions)):

        curr_pos = begin + skipped_gaps + i
        if direction == "upstream":     # go backwards
            curr_pos = end - 1 - skipped_gaps - i

        if skip_gaps:
            # keep going until no gap
            if direction == "upstream":
                while curr_pos >= 0 and (list_alignments[0][curr_pos] == "-" or gap_fraction_in_pos(list_alignments, curr_pos) > max_fraction_of_gaps_in_pos):
                    curr_pos -= 1
                    skipped_gaps += 1

                if curr_pos < 0:
                    if only_full_length:
                        raise ValueError("Need enough region")
                    break
            else:
                while curr_pos < end and (list_alignments[0][curr_pos] == "-" or gap_fraction_in_pos(list_alignments, curr_pos) > max_fraction_of_gaps_in_pos):
                    curr_pos += 1
                    skipped_gaps += 1

                if curr_pos == len(list_alignments[0]):
                    if only_full_length:
                        raise ValueError("Need enough region")
                    break

        if score_on_all_pairs:
            for row_a in range(len(list_alignments) - 1):
                curr_letter_a = list_alignments[row_a][curr_pos]
                for row_b in range(row_a + 1, len(list_alignments)):
                    curr_letter_b = list_alignments[row_b][curr_pos]
                    if scorer is None:
                        if curr_letter_a == curr_letter_b and curr_letter_a != "-":
                            matched += 1
                        total += 1
                    else:
                        matched += scorer.score(curr_letter_a, curr_letter_b)
                        total += 1
        else:
            # for all sequences in that position
            curr_letter = list_alignments[0][curr_pos]

            for j in range(1, len(list_alignments)):
                if scorer is None:
                    if curr_letter == list_alignments[j][curr_pos]:
                        matched += 1
                    total += 1

                else:
                    matched += scorer.score(curr_letter, list_alignments[j][curr_pos])
                    total += 1


    if total == 0:
        return 0

    return matched / float(total)


def select_start_position_from_msa_for_lorf(alignments, **kwargs):
    """

    :param alignments:
    :type alignments: Bio.Align.MultipleSequeuence
    :param kwargs:
    :type kwargs: Dict[str, Any]
    :return:
    """

    logger.debug("Func: select-start-positions-from-msa-for-lorf")

    region_length = get_value(kwargs, "region_length", 10)
    threshold = get_value(kwargs, "threshold", 0.5)
    scorer = get_value(kwargs, "scorer", None)

    msa_options = get_value(kwargs, "msa_options", None)
    pos_of_upstream_in_msa = get_value(kwargs, "pos_of_upstream_in_msa", None)



    score_on_all_pairs = get_value(kwargs, "score_on_all_pairs", False)

    list_alignments = [a.seq._data for a in alignments]

    # find the positions of the first two candidate starts in the query
    start = 0
    end = len(list_alignments[0])

    candidates = list()
    for i in range(start, end):

        if list_alignments[0][i].isupper():
            candidates.append(i)

    logger.debug("Total number of candidate 5prime ends: {}".format(len(candidates)))

    if len(candidates) == 1:
        logger.debug("Only one candidate. Selecting position {}".format(candidates[0]))
        return candidates[0]

    # if not enough, return None
    distance_between_candidates = len(list_alignments[0][candidates[0]:candidates[1]].replace("-", ""))
    if distance_between_candidates < region_length:
        logger.debug("Distance between first and second candidate is {} < {}".format(distance_between_candidates, region_length))
        return None

    # otherwise, compute conservation
    conservation = None
    if msa_options is not None and msa_options.safe_get("search-use-upstream-for-lorf") and pos_of_upstream_in_msa is not None:
        # only consider candidate if it is not overlapping and region above isn't overlapping either

        if candidates[1] > pos_of_upstream_in_msa + region_length:
            conservation = compute_conservation_in_region(list_alignments, candidates[1] - region_length, candidates[1],
                                                          direction="upstream", skip_gaps=True, only_full_length=True,
                                                          score_on_all_pairs=score_on_all_pairs,
                                                          scorer=scorer
                                                          )
    else:
        conservation = compute_conservation_in_region(list_alignments, candidates[1]-region_length, candidates[1],
                                                      direction="upstream", skip_gaps=True, only_full_length=True,
                                                      score_on_all_pairs=score_on_all_pairs,
                                                      scorer=scorer
                                                      )

    logger.debug("Conservation upstream of second candidate is {} (threshold={})".format(conservation, threshold))

    # if conservation exists, return
    if conservation is not None and conservation > threshold:
        logger.debug("Conservation > threshold. Selecting LORF at position {}".format(candidates[0]))
        return candidates[0]

    return None


def get_candidates_before_any_conservation_block(alignments, **kwargs):
    """

    :param alignments:
    :type alignments: Bio.Align.MultipleSequence
    :param kwargs:
    :return: List of candidate positions before any that have a block of upstream conservation
    :rtype: List[int]
    """

    region_length = get_value(kwargs, "region_length", 10)
    threshold = get_value(kwargs, "threshold", 0.9)
    max_fraction_of_gaps_in_pos = get_value(kwargs, "max_fraction_of_gaps_in_pos", None)

    list_alignments = [a.seq._data for a in alignments]

    # find the positions of the first two candidate starts in the query
    start = 0
    end = len(list_alignments[0])

    candidates = list()

    # reach first candidate
    pos_of_first_candidate = 0
    for i in range(start, end):

        if list_alignments[0][i].isupper():
            pos_of_first_candidate = i
            break

    def gap_fraction_in_pos(list_alignments, pos):
        total = len(list_alignments)
        if total == 0:
            return 0.0

        return sum(1 for i in range(len(list_alignments)) if list_alignments[i][pos] == "-") / float(total)

    def gap_fraction_in_positions(list_alignments, pos_start, pos_end):
        # positions end is exclusive
        total = len(list_alignments) * (pos_end - pos_start)
        if total == 0:
            return 0.0

        numerator = 0
        for pos in range(pos_start, pos_end):
            numerator += sum(1 for i in range(len(list_alignments)) if list_alignments[i][pos] == "-")

        return numerator / float(total)


    # find: first block without many gaps
    curr_pos = pos_of_first_candidate
    for i in range(curr_pos, end):

        if gap_fraction_in_positions(list_alignments, i, i + region_length) < max_fraction_of_gaps_in_pos:
            if gap_fraction_in_pos(list_alignments, i) > max_fraction_of_gaps_in_pos:
                curr_pos = i
                break

        if list_alignments[0][i].isupper():
            candidates.append(i)



    # compute conservation at this position
    for i in range(curr_pos, end):

        if list_alignments[0][i].isupper():
            candidates.append(i)

        if list_alignments[0][i] == "-" or gap_fraction_in_pos(list_alignments, i) > max_fraction_of_gaps_in_pos:
            continue

        try:
            conservation = compute_conservation_in_region(list_alignments, i, i+region_length,
                                                          direction="downstream", skip_gaps=True, only_full_length=True,)

            if conservation > threshold:
                break

        except ValueError:
            # not enough sequence
            pass

    return candidates




def get_candidates_without_upstream_conservation(alignments, **kwargs):
    """

    :param alignments:
    :type alignments: Bio.Align.MultipleSequence
    :param kwargs:
    :return: List of candidate positions before any that have a block of upstream conservation
    :rtype: List[int]
    """

    logger.debug("Func: get-candidates-without-upstream-conservation")

    region_length = get_value(kwargs, "region_length", 10)
    threshold = get_value(kwargs, "threshold", 0.5)
    score_on_all_pairs = get_value(kwargs, "score_on_all_pairs", False)
    scorer = get_value(kwargs, "scorer", None)

    msa_options = get_value(kwargs, "msa_options", None)
    pos_of_upstream_in_msa = get_value(kwargs, "pos_of_upstream_in_msa", None)

    if pos_of_upstream_in_msa is not None and pos_of_upstream_in_msa < 0:
        pos_of_upstream_in_msa = None


    list_alignments = [a.seq._data for a in alignments]

    # find the positions of the first two candidate starts in the query
    start = 0
    end = len(list_alignments[0])

    candidates = list()
    for i in range(start, end):

        if list_alignments[0][i].isupper():
            try:
                # if candidate in upstream region, add it and move on to next
                if msa_options is not None and msa_options.safe_get(
                        "search-use-upstream-for-candidate-selection") and pos_of_upstream_in_msa is not None and pos_of_upstream_in_msa > -1 and i <= pos_of_upstream_in_msa:
                    logger.debug("Candidate in upstream gene. Add and move on")
                    candidates.append(i)
                else:
                    conservation = compute_conservation_in_region(list_alignments, i - region_length, i,
                                                              direction="upstream", skip_gaps=True, only_full_length=True,
                                                                  score_on_all_pairs=score_on_all_pairs,
                                                                  scorer=scorer)

                    logger.debug("Candidate pos {}, conservation = {}".format(i, conservation))
                    if conservation < threshold:
                        logger.debug("Conservation < threshold ({}). Add to list".format(threshold))
                        candidates.append(i)
                    else:
                        logger.debug("Conservation >= threshold ({}). Stop".format(threshold))
                        break

            except ValueError:
                # not enough sequence
                logger.debug("Not enough upstream for candidate {}. Add to list".format(i))
                candidates.append(i)

    return candidates













def select_start_position_from_msa(alignments, msa_options, **kwargs):
    # type: (Bio.Align.MultipleSeqAlignment, MSAOptions, Dict[str, Any]) -> int

    logger.debug("Func: Select start position from MSA")

    list_alignments = [a.seq._data for a in alignments]     # type: list[str]

    score_on_all_pairs = get_value(kwargs, "score_on_all_pairs", False)

    threshold = get_value(kwargs, "threshold", 0.5)
    start_from_pos = get_value(kwargs, "start_from_pos", 0)
    end_at_pos = get_value(kwargs, "end_at_pos", len(list_alignments[0]))

    pos_of_upstream_in_msa = get_value(kwargs, "pos_of_upstream_in_msa", None)

    upstream_block_scorer = get_value(kwargs, "upstream_block_scorer", None)

    if msa_options is None:
        raise ValueError("MSA Options cannot be None")

    pos_to_num = dict()

    start_pos = None
    curr_support = 0

    # check upstream gene
    if msa_options.safe_get("search-use-upstream-for-start-selection") and pos_of_upstream_in_msa is not None:
        if pos_of_upstream_in_msa > 0:
            curr_msa_options = copy.deepcopy(msa_options)
            curr_msa_options["search-penalize-no-sequence"] = 1
            curr_msa_options["search-start-selection-threshold"] = 0.3
            # check if there is a candidate very close to upstream gene
            # check upstream
            aa_checked = 0
            curr_pos = pos_of_upstream_in_msa
            while aa_checked < 2:

                if list_alignments[0][curr_pos].isupper():
                    logger.debug("Found at {} from upstream".format(pos_of_upstream_in_msa - curr_pos))
                    if curr_msa_options.safe_get("search-5prime-near-upstream-is-conserved"):
                        num_upper = count_num_upper(list_alignments, curr_pos, curr_msa_options)
                        if  num_upper / float(len(list_alignments)) > curr_msa_options["search-start-selection-threshold"]:
                            return curr_pos
                    else:
                        return curr_pos
                aa_checked += 1

                while curr_pos >= 0 and list_alignments[0][curr_pos] == "-":
                    curr_pos -= 1
                curr_pos -= 1

            aa_checked = 0
            curr_pos = pos_of_upstream_in_msa
            while aa_checked < 2:

                if list_alignments[0][curr_pos].isupper():
                    logger.debug("Found at {} from upstream".format(pos_of_upstream_in_msa - curr_pos))
                    if curr_msa_options.safe_get("search-5prime-near-upstream-is-conserved"):
                        num_upper = count_num_upper(list_alignments, curr_pos, curr_msa_options)
                        if  num_upper / float(len(list_alignments)) > curr_msa_options["search-start-selection-threshold"]:
                            return curr_pos
                    else:
                        return curr_pos
                aa_checked += 1
                curr_pos += 1
                while curr_pos < len(list_alignments[0]) and list_alignments[0][curr_pos] == "-":
                    curr_pos += 1


    for i in range(int(start_from_pos),int( min(end_at_pos, len(list_alignments[0])))):

        num_upper = 0

        if list_alignments[0][i].isupper():

            num_upper = count_num_upper(list_alignments, i, msa_options)

            logger.debug("Position {}: total seqs = {}, num upper = {}".format(i, len(list_alignments), num_upper))

            pos_to_num[i] = num_upper



            if num_upper / float(len(list_alignments)) > threshold:
                start_pos = i
                curr_num_upper = num_upper
                curr_support = num_upper / float(len(list_alignments))
                logger.debug("Position {}: {} > {}".format(i, num_upper / float(len(list_alignments)), threshold))
                break

    consider_upstream_conservation = True
    # try msa options
    if start_pos is not None:

        # search downstream
        if msa_options["search-better-downstream-aa"] is not None:
            if msa_options["search-better-downstream-aa"] > 0:
                logger.debug("Search for better downstream: width = {}".format(msa_options["search-better-downstream-aa"]))

                msa_options_tmp = copy.deepcopy(msa_options)
                # msa_options_tmp["search-better-downstream-aa"] = 0
                msa_options_tmp["search-upstream-of-conserved-region"] = None
                new_sp = select_start_position_from_msa(
                    alignments, msa_options_tmp,
                    threshold=threshold,
                    start_from_pos=start_pos + 1,
                    end_at_pos=start_pos + 1 + msa_options["search-better-downstream-aa"]
                )

                logger.debug("Potential new position: {}".format(new_sp))

                if new_sp is not None:
                    new_num_upper = count_num_upper(list_alignments, new_sp, msa_options)
                    logger.debug(
                        "Position {}: total seqs = {}, num upper = {}".format(
                            new_sp, len(list_alignments), new_num_upper)
                    )

                    new_support = new_num_upper / float(len(list_alignments))

                    current_sp = start_pos

                    if msa_options["search-favor-m"]:

                        logger.debug("Position {}: {}".format(current_sp, list_alignments[0][current_sp]))
                        logger.debug("Position {}: {}".format(new_sp, list_alignments[0][new_sp]))

                        # if new position is M, select it
                        if list_alignments[0][start_pos] != "M" and list_alignments[0][new_sp] == "M":
                            start_pos = new_sp
                        # if current is M and new isn't M, keep original
                        elif list_alignments[0][start_pos] == "M" and list_alignments[0][new_sp] != "M":
                            pass
                        # otherwise, select based on support
                        else:
                            if new_support > curr_support:
                                logger.debug("Better support, choosing {}".format(new_sp))
                                start_pos = new_sp

                            if msa_options.safe_get("search-stringent-check-when-equal-support") and new_num_upper == curr_num_upper:
                                logger.debug("Positions {} and {} have same support. Computing support with neighbor=0".format(
                                    current_sp, new_sp
                                ))
                                tmp_msa_options = copy.deepcopy(msa_options)
                                tmp_msa_options["search-neighbor"] = 0
                                start_pos = new_sp if count_num_upper(list_alignments, new_sp, tmp_msa_options) > count_num_upper(list_alignments, current_sp, tmp_msa_options) else current_sp

                    else:
                        if new_support > curr_support:
                            logger.debug("Better support, choosing {}".format(new_sp))
                            start_pos = new_sp

                        if msa_options.safe_get("search-stringent-check-when-equal-support") and new_support == curr_support:
                            logger.debug(
                                "Positions {} and {} have same support. Computing support with neighbor=0".format(
                                    current_sp, new_sp
                                ))
                            tmp_msa_options = copy.deepcopy(msa_options)
                            tmp_msa_options["search-neighbor"] = 0
                            start_pos = new_sp if count_num_upper(list_alignments, new_sp,
                                                                  tmp_msa_options) > count_num_upper(list_alignments,
                                                                                                     current_sp,
                                                                                                     tmp_msa_options) else current_sp

                    # if start position was changed, check if region between is conserved
                    if msa_options.safe_get("search-better-downstream-aa-if-not-conserved"):
                        if current_sp != start_pos and start_pos - current_sp > 1:
                            consider_upstream_conservation = False
                            logger.debug("Better downstream selected. Check conservation")
                            conservation = compute_conservation_in_region(
                                list_alignments,
                                current_sp+1, start_pos,
                                direction="downstream",
                                score_on_all_pairs=score_on_all_pairs,
                                scorer=upstream_block_scorer
                            )

                            logger.debug("Conservation = {}".format(conservation))

                            if msa_options.safe_get("search-better-downstream-by-margin"):

                                margin = msa_options.safe_get("search-better-downstream-by-margin")
                                logger.debug("Search better downstream by margin: {}".format(margin))

                                def search_better_downstream_by_margin(l_new_support, l_curr_support, start_pos):
                                    if True or new_num_upper == curr_num_upper:
                                        tmp_msa_options = copy.deepcopy(msa_options)
                                        tmp_msa_options["search-neighbor"] = 0
                                        l_new_support = count_num_upper(list_alignments, new_sp, tmp_msa_options) / float(len(list_alignments))
                                        l_curr_support = count_num_upper(list_alignments, current_sp,
                                                                        tmp_msa_options) / float(len(list_alignments))
                                    if l_new_support > l_curr_support + margin:
                                        logger.debug("{} > {} + {}".format(l_new_support, l_curr_support, margin))
                                        start_pos = new_sp
                                        logger.debug("Position = {}".format(new_sp))
                                        pass
                                    else:
                                        logger.debug("{} < {} + {}".format(l_new_support, l_curr_support, margin))
                                        if conservation > 0.5:
                                            logger.debug("{} > 0.5".format(conservation))
                                            start_pos = current_sp
                                            logger.debug(
                                                "Conservation too high - go back to origin: {}".format(current_sp))

                                    return start_pos


                                start_pos = search_better_downstream_by_margin(new_support, curr_support, start_pos)

                            elif conservation > 0.5:
                                start_pos = current_sp
                                logger.debug("Conservation too high - go back to origin: {}".format(current_sp))

        # now search upstream: if current start has a strong conserved upstream region,
        # then loop back and see if there's a better candidate upstream by relaxing
        # decision
        if consider_upstream_conservation and msa_options["search-upstream-of-conserved-region"] is not None \
                and msa_options["search-upstream-of-conserved-region-threshold"] is not None:

            logger.debug("Search upstream of conserved region")
            region_length = msa_options["search-upstream-of-conserved-region"]
            conservation_threshold = msa_options["search-upstream-of-conserved-region-threshold"]

            # do this until no upstream region is conserved or no more candidates
            done = False

            curr_start_pos = start_pos
            selection_has_conserved_upstream = False
            while not done:

                begin = curr_start_pos - region_length
                end = curr_start_pos
                try:
                    conservation = compute_conservation_in_region(list_alignments, begin, end, skip_gaps=True, only_full_length=True,
                                                                  score_on_all_pairs=score_on_all_pairs,
                                                                  scorer=upstream_block_scorer) # FIXME
                except ValueError:
                    conservation = 0.0

                if conservation >= conservation_threshold:
                    selection_has_conserved_upstream = True
                    # find closest upstream start
                    i = curr_start_pos - 1

                    found_new_start_at = None
                    while i > 0:
                        if list_alignments[0][i].isupper():
                            found_new_start_at = i
                            break
                        i -= 1

                    # if no upstream start is found, we're done
                    if found_new_start_at is None:
                        done = True
                    else:

                        msa_options_tmp = copy.deepcopy(msa_options)
                        msa_options_tmp["search-better-downstream-aa"] = 0
                        msa_options_tmp["search-neighbor"] = 4
                        if msa_options.safe_get("search-neighbor-relaxed") is not None:
                            msa_options_tmp["search-neighbor"] = msa_options.safe_get("search-neighbor-relaxed")
                        num_upper = count_num_upper(list_alignments, i, msa_options_tmp)

                        pos_to_num[i] = num_upper

                        if num_upper / float(len(list_alignments)) > threshold:
                            start_pos = i
                            curr_start_pos = start_pos
                        # if score isn't good
                        else:
                            # FIXME: here what's the effect?
                            # return None
                            curr_start_pos = i
                            pass

                else:
                    selection_has_conserved_upstream = False
                    done = True

            if msa_options.safe_get("search-never-stay-at-conserved-upstream"):
                if selection_has_conserved_upstream:
                    return None
    return start_pos



def add_start_position_to_msa_alignments(alignments, position, position_of_upstream=None):
    # type: (Bio.Align.MultipleSeqAlignment, int, Union[int, None]) -> Bio.Align.MultipleSeqAlignment

    if position is None:
        start_sequence = "-" * alignments.get_alignment_length()
    else:
        start_sequence = "-" * (position) + "M" + "-" * (alignments.get_alignment_length() - position - 1)

    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq

    seq_records = [SeqRecord(Seq(start_sequence), id="#selected")]

    # if position_of_upstream is not None:
    #     if position_of_upstream < 0:
    #         upstream_sequence = "-" * alignments.get_alignment_length()
    #     else:
    #         upstream_sequence = "-" * (position_of_upstream) + "X" + "-" * (alignments.get_alignment_length() - position_of_upstream - 1)
    #
    #     seq_records.append(SeqRecord(Seq(upstream_sequence), id="#upstream"))

    for a in alignments:
        seq_records.append(a)

    return Bio.Align.MultipleSeqAlignment(seq_records)


def select_and_add_start_position_alignments(alignments, msa_options):
    # type: (Bio.Align.MultipleSeqAlignment, MSAOptions) -> Bio.Align.MultipleSeqAlignment

    position = select_start_position_from_msa(alignments, msa_options=msa_options)

    alignments = add_start_position_to_msa_alignments(alignments, position)

    return alignments


def filter_df_by_upstream_to_msa(df, msa_options, **kwargs):
    # type: (pd.DataFrame, MSAOptions, Dict[str, Any]) -> pd.DataFrame

    suffix_upstream = get_value(kwargs, "suffix_upstream", "upstream-distance-msa")
    suffix_gene_length = get_value(kwargs, "suffix_gene_length", "gene-length-msa")
    filter_stats = get_value(kwargs, "filter_stats", None)

    column_q_upstream = "q-{}".format(suffix_upstream)
    column_t_upstream = "t-{}".format(suffix_upstream)

    if msa_options["filter-min-upstream-distance"] is not None:

        before = len(df)

        df = df[df[column_q_upstream] >= msa_options["filter-min-upstream-distance"]]
        df = df[df[column_t_upstream] >= msa_options["filter-min-upstream-distance"]]

        after = len(df)

        if filter_stats:
            filter_stats["filter-min-upstream-distance"] = before - after

    if msa_options["filter-gene-length-max-percent-difference"] is not None:

        column_q_gene_length = "q-{}".format(suffix_gene_length)
        column_t_gene_length = "t-{}".format(suffix_gene_length)

        before = len(df)

        val = msa_options["filter-gene-length-max-percent-difference"]

        df = df[((100.0 * abs(df[column_q_gene_length] - df[column_t_gene_length ])) / df[column_q_gene_length]) < val]

        after = len(df)

        if filter_stats:
            filter_stats["filter-gene-length-max-percent-difference"] = before - after

    return df



def filter_df(df, msa_options, **kwargs):
    # type: (pd.DataFrame, MSAOptions, Dict[str, Any]) -> pd.DataFrame

    column_distance = get_value(kwargs, "column_distance", "k2p-distance")
    filter_stats = get_value(kwargs, "filter_stats", None)
    filter_non_group_only = get_value(kwargs, "filter_non_group_only", False)

    sbsp_general.general.df_add_5prime_3prime_key(df, "q-")

    # remove frameshifted
    before = len(df)
    not_frame_shifted = ((df["q-right"]-df["q-left"]+1) % 3 == 0) & ((df["t-right"]-df["t-left"]+1) % 3 == 0)
    df = df[not_frame_shifted]
    after = len(df)

    if filter_stats is not None:
        filter_stats["filter-frame-shifted"] = before - after
    else:
        logger.critical("Filter ({}): {} - {} = {}".format("frame-shifted", before, after, before - after))


    # max distance
    if len(df) > 0 and msa_options["filter-max-distance"] is not None:
        if column_distance not in df.columns.values:
            logger.warning("Can't filter by distance. Unknown column ({})".format(column_distance))
        else:
            before = len(df)
            filter_indeces = df[column_distance] < msa_options["filter-max-distance"]
            df = df[filter_indeces]
            after = len(df)

            if filter_stats is not None:
                filter_stats["filter-max-distance"] = before - after
            else:
                logger.critical("Filter ({}): {} - {} = {}".format("filter-max-distance", before, after, before - after))

    # min distance
    if len(df) > 0 and msa_options["filter-min-distance"] is not None:
        if column_distance not in df.columns.values:
            logger.warning("Can't filter by distance. Unknown column ({})".format(column_distance))
        else:
            before = len(df)
            filter_indeces = df[column_distance] > msa_options["filter-min-distance"]
            df = df[df[column_distance] > msa_options["filter-min-distance"]]
            after = len(df)

            if filter_stats is not None:
                filter_stats["filter-min-distance"] = before - after
            else:
                logger.critical("Filter ({}): {} - {} = {}".format("filter-min-distance", before, after, before - after))

    # equal distances
    if len(df) > 0 and msa_options["filter-orthologs-with-equal-kimura"] is not None and column_distance in df.columns.values:

        column_distance_rounded = "column_distance_rounded"
        df[column_distance_rounded] = df[column_distance].round(msa_options["filter-orthologs-with-equal-kimura"])
        # df.round({column_distance_rounded: msa_options["filter-orthologs-with-equal-kimura"]})

        before = len(df)

        sbsp_general.general.df_add_3prime_key(df, "q-")
        df = df.groupby(["q-3prime", column_distance_rounded], as_index=False).agg("first").reset_index(drop=True)

        df.drop("q-3prime", inplace=True, axis=1)
        df.drop(column_distance_rounded, inplace=True, axis=1)

        after = len(df)

        if filter_stats is not None:
            filter_stats["filter-orthologs-with-equal-kimura"] = before - after
        else:
            logger.critical("Filter ({}): {} - {} = {}".format("filter-orthologs-with-equal-kimura", before, after, before - after))

    if len(df) > 0 and msa_options["filter-max-number-orthologs"] is not None and not filter_non_group_only:

        number = msa_options["filter-max-number-orthologs"]

        before = len(df)

        sbsp_general.general.df_add_3prime_key(df, "q-")
        # df = df.groupby(["q-3prime"], as_index=False).head(number)
        print ("Length = {}".format(len(df)))

        if len(df) > number:
            sample_method = msa_options.safe_get("filter-sample-from-set")
            if sample_method is None:
                sample_method = "uniform"

            if sample_method == "uniform":
                df = df.groupby(["q-3prime"], as_index=False).apply(lambda x: x.sample(n=number, random_state=1) if len(x) > number else x)
            elif sample_method == "closest-distance":
                logger.debug("Closest samples by distance")
                df = df.groupby(["q-3prime"], as_index=False).apply(lambda x: x.sort_values("kimura").head(number)).reset_index(drop=True)


        df.drop("q-3prime", inplace=True, axis=1)

        after = len(df)
        if filter_stats is not None:
            filter_stats["filter-max-number-orthologs"] = before - after
        else:
            logger.critical(
                "Filter ({}): {} - {} = {}".format("filter-max-number-orthologs", before, after, before - after))

    if len(df) > 0:
        sbsp_general.general.df_add_3prime_key(df, "q-")

    #  # filter based on upstream
    # filter-min-upstream-distance: null
    # filter-close-max-upstream-distance: null
    # filter-far-min-upstream-distance: null
    #
    # # filter by gene length
    # filter-gene-length-max-percent-difference: null

    # if msa_options["filter-min-upstream-distance"]


    return df


def setup_directory_for_msa_outputs(env, dn_msa_output):
    # type: (Dict[str, Any], str) -> str
    pd_msa_output = None
    if dn_msa_output is not None:
        pd_msa_output = os.path.join(env["pd-work"], dn_msa_output)

        try:
            import sbsp_io.general
            sbsp_io.general.mkdir_p(pd_msa_output)
        except OSError as e:
            e.message = e.message + "\nCouldn't create MSA output directory"
            raise e
    else:
        pd_msa_output = env["pd-work"]

    return pd_msa_output


def get_label_from_start_position_in_msa(df, q_sequence_aa, q_alignment_aa,
                                         q_pos_of_previous_start_in_msa,
                                         q_pos_of_current_start_in_msa,
                                         source="q",
                                         row_num=0):
    # type: (pd.DataFrame, str, str, int, int) -> sbsp_general.labels.Label

    num_gaps_before_pos = q_alignment_aa[:q_pos_of_current_start_in_msa].count("-")
    pos_in_non_gapped = q_pos_of_current_start_in_msa - num_gaps_before_pos + \
                        q_sequence_aa[:q_pos_of_current_start_in_msa - num_gaps_before_pos].count("*")

    # position relative to 5 prime (negative is upstream, positive is downstream)
    pos_relative_to_5prime = pos_in_non_gapped - q_pos_of_previous_start_in_msa

    # update data frame rows
    q_label = sbsp_general.dataframe.df_get_label_from_row(
        df, df.index[row_num], source, None)  # all indeces should have the same query label # FIXME: None for suffix

    # shift 5 prime
    sbsp_general.labels.shift_5prime(q_label, pos_relative_to_5prime * 3)

    return q_label


def print_alignment_to_file(alignments, start_position_in_msa, df, pf_msa_out, pos_of_upstream_in_msa=None):
    # type: (Bio.Align.MultipleSeqAlignment, int, pd.DataFrame, str, str) -> None

    """
    Print alignments to file and update dataframe.
    The number of sequences in the alignments should be len(df) + 1, where the first alignment sequence
    is the query gene. Furthermore, the order of entries in df should correspond to the order
    of sequences in alignments[1:]
    :param alignments:
    :param start_position_in_msa:
    :param df:
    :return:
    """

    alignments = add_start_position_to_msa_alignments(alignments, start_position_in_msa, position_of_upstream=pos_of_upstream_in_msa)
    stdout = convert_alignments_to_str(alignments)
    indeces_of_orthologs = list()
    if df is not None:
        indeces_of_orthologs = df.index

    with open(pf_msa_out, "w") as f_msa_out:
        f_msa_out.write(stdout)

        # add file
        pos_num = 1
        for row_num in indeces_of_orthologs:
            df.at[row_num, "pf-msa-output"] = pf_msa_out
            df.at[row_num, "t-pos-in-msa-output"] = pos_num

            pos_num += 1


def compute_match_in_area(list_alignments, pos_of_5prime, region_aa, direction,
                          frame=None):

    from sbsp_general.general import except_if_not_in_set

    except_if_not_in_set(direction, {"upstream", "downstream"})

    def next_position(pos, direction, step=1):
        if direction == "upstream":
            pos -= step
        else:
            pos += step
        return pos

    curr_pos = pos_of_5prime

    # shift to frame if set
    step = 1
    if frame is not None:
        step = 3
        curr_pos += frame
        region_aa /= 3

    alignment_length = len(list_alignments[0])
    num_none_gapped = 0
    match = 0
    total = 0
    while num_none_gapped < region_aa:

        done = False            # set to true when no more positions are left to analyze

        # skip all gaps in query and make sure we are in bounds
        while True:
            curr_pos = next_position(curr_pos, direction, step=step)

            if curr_pos < 0 or curr_pos >= alignment_length:
                done = True
                break

            if list_alignments[0][curr_pos] != "-":
                break

        # if no positions left, break
        if done:
            break

        # otherwise, count matches to query
        q_element = list_alignments[0][curr_pos].lower()

        curr_match = 0
        for i in range(1, len(list_alignments)):
            if q_element == list_alignments[i][curr_pos].lower():
                match += 1
                curr_match += 1

            total += 1

        num_none_gapped += 1

    if total == 0:
        return 0.0
    else:
        return match / float(total)






# def add_features()
def compute_post_alignment_features(
        df,
        list_alignments_aa,
        list_alignments_nt,
        pos_aa,
        pos_nt,
        msa_options,
        region_aa
):
    # type:

    support = count_num_upper(list_alignments_aa, pos_aa, msa_options)

    downstream_match = compute_match_in_area(list_alignments_aa, pos_aa, region_aa, "downstream")
    upstream_match = compute_match_in_area(list_alignments_aa, pos_aa, region_aa, "upstream")

    df["f-support"] = support
    df["f-conservation-downstream-match"] = downstream_match
    df["f-conservation-upstream-match"] = upstream_match

    for frame in range(3):
        df["f-conservation-downstream-match-{}".format(frame)] = compute_match_in_area(list_alignments_nt, pos_nt, region_aa*3, "downstream", frame=frame)
        df["f-conservation-upstream-match-{}".format(frame)] = compute_match_in_area(list_alignments_nt, pos_nt, region_aa*3, "upstream", frame=frame)

    # ratios
    for frame_1 in range(3):
        for frame_2 in range(3):
            # downstream
            df["f-conservation-downstream-match-{}-div-{}".format(frame_1, frame_2)] = \
                df["f-conservation-downstream-match-{}".format(frame_1)] / df[
                    "f-conservation-downstream-match-{}".format(frame_2)]

            # upstream
            df["f-conservation-upstream-match-{}-div-{}".format(frame_1, frame_2)] = \
                df["f-conservation-upstream-match-{}".format(frame_1)] / df[
                    "f-conservation-upstream-match-{}".format(frame_2)]


        # downstream / upstream
        df["f-conservation-ratio-match-{}".format(frame_1)] = \
            df["f-conservation-downstream-match-{}".format(frame_1)] / df[
                "f-conservation-upstream-match-{}".format(frame_1)]






def convert_msa_from_aa_to_nt(alignments_aa, list_sequences_nt):
    # type: (Bio.Align.MultipleSeqAlignment, List[str]) -> Bio.Align.MultipleSeqAlignment

    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq

    seq_records = list()
    i = 0
    for a in alignments_aa:
        align_nt = sbsp_general.general.add_gaps_to_nt_based_on_aa(
            list_sequences_nt[i],      # remove stop
            a.seq._data,
            preserve_case=True
        )

        seq_records.append(SeqRecord(Seq(align_nt), id=a.id))

        i += 1

    return Bio.Align.MultipleSeqAlignment(seq_records)


def compute_kimura_matrix(list_sequences_aligned_nt):
    # type: (List[str]) -> np.ndarray

    from sbsp_alg.phylogeny import k2p_distance
    num_sequences = len(list_sequences_aligned_nt)

    output = np.zeros((num_sequences, num_sequences), dtype=float)

    # compute pairwise kimura
    # diagonal is zero, and matrix is symmetric (I think :))

    for i in range(1, num_sequences):
        for j in range(i):

            val = k2p_distance(list_sequences_aligned_nt[i], list_sequences_aligned_nt[j])

            output[i, j] = val
            output[j, i] = val

    return output

def get_indices_after_filtering_random(edge_mat, min_range, max_range):
    # type: (np.ndarray, Union[float, None], Union[float, None]) -> List[int]

    rows, cols = edge_mat.shape

    removed_nodes = set()

    for i in range(2, rows):            # skip query

        if i in removed_nodes:
            continue

        for j in range(1, i):           # skip query

            if i == j:
                continue

            if j in removed_nodes:
                continue

            remove = False
            if min_range is not None:
                if edge_mat[i][j] < min_range:
                    remove = True

            if max_range is not None:
                if edge_mat[i][j] > max_range:
                    remove = True

            if remove:
                # delete entire row column for the ith
                # edge_mat[i, :] = 0
                # edge_mat[:, i] = 0
                removed_nodes.add(i)

    return list(set(range(rows)) - removed_nodes)

def get_indices_after_filtering(edge_mat, min_range=None, max_range=None, strategy="random"):
    # type: (np.ndarray, Union[float, None], Union[float, None], str) -> List[int]

    from sbsp_general.general import except_if_not_in_set

    except_if_not_in_set(strategy, ["random", "min-number-of-deletions", "max-number-of-deletions"])

    if strategy in ["min-number-of-deletions", "max-number-of-deletions"]:
        raise NotImplementedError("Strategy ({}) not yet implemented".format(strategy))

    if strategy == "random":
        return get_indices_after_filtering_random(edge_mat, min_range, max_range)




def filter_by_pairwise_kimura_from_msa(list_sequences_aligned_nt, msa_options):
    # type: (List[str], MSAOptions) -> Tuple(List[str], List[int])

    output = list()

    kimura_mat = compute_kimura_matrix(list_sequences_aligned_nt)

    min_val = 0.001
    max_val = 0.4
    if msa_options.safe_get("filter-pairwise-kimura-min-max"):
        values = msa_options.safe_get("filter-pairwise-kimura-min-max")

        min_val = float(values[0]) if values[0] is not None else min_val
        max_val = float(values[1]) if values[1] is not None else max_val

    indices_to_keep = get_indices_after_filtering(kimura_mat, min_val, max_val, "random")

    for i in sorted(indices_to_keep):
        output.append(list_sequences_aligned_nt[i])

    return output, sorted(indices_to_keep)


def filter_sequences_that_introduce_gaps_in_msa(list_sequences_aligned, msa_options):
    """

    :param list_sequences_aligned:
    :type list_sequences_aligned: List[str]
    :param msa_options:
    :type msa_options: MSAOptions
    :return:
    """


    # sketch: get to the first start codon in methianine
    # from then (up to 50% of alignment length), if a large gap observed, remove sequences that caused it

    params = msa_options.safe_get("filter-remove-sequences-that-introduce-gaps")
    gap_width = params[0]
    seq_frac = params[1]

    alignment_length = len(list_sequences_aligned[0])

    num_sequences_aligned = len(list_sequences_aligned)

    # get to first start codon
    first_start_codon_position = 0
    for i in range(alignment_length):
        if list_sequences_aligned[0][i].isupper():
            first_start_codon_position = i
            break

    end_search = int(alignment_length/2.0) - gap_width
    for i in range(first_start_codon_position, end_search):

        # if block
        block_detected = list_sequences_aligned[0][i:i+gap_width].count("-") == gap_width

        sequences_that_contribute_to_block = list()

        if block_detected:

            for j in range(1, num_sequences_aligned):

                if list_sequences_aligned[j][i:i+gap_width].count("-") == 0:
                    sequences_that_contribute_to_block.append(j)

            num_sequences_that_contribute_to_block = len(sequences_that_contribute_to_block)
            fraction = float(num_sequences_that_contribute_to_block) / num_sequences_aligned

            if fraction != 0 and fraction < seq_frac:
                indices_to_keep = set(range(num_sequences_aligned)).difference(sequences_that_contribute_to_block)
                output = [list_sequences_aligned[i] for i in indices_to_keep]

                return output, sorted(indices_to_keep)

    return list_sequences_aligned, range(num_sequences_aligned)


def get_elements_at_indices(original_list, index_list):
    # type: (List[Any], List[int]) -> List[Any]
    return [original_list[i] for i in index_list]


def find_rightmost_by_standard_aa_score(alignments, candidates, threshold=0.5):
    # type: (Bio.Align.MultipleSeqAlignment, List[int], float) -> Union[int, None]

    logger.debug("Func: find-rightmost-by-standard-aa-score")
    import sbsp_ml.msa_features
    for i in reversed(candidates):

        penalized_start_score = sbsp_ml.msa_features.compute_simple_saas(alignments, i, 0, 0)

        logger.debug("Candidate {}, SAAS = {}".format(i, penalized_start_score))

        if penalized_start_score < threshold:
            logger.debug("SAAS < threshold ({}). Select it".format(penalized_start_score))
            return i

    logger.debug("No candidate found with low SAAS")
    return None


def compute_position_of_upstream_of_query_in_msa(msa_t, q_pos_5prime_in_msa_no_gap, distance_to_upstream_nt, index_in_msa=0):
    # type: (MSAType, int, float, int) -> int

    """
    Computes the position of the upstream gene in the MSA
    :param msa_t:
    :param q_pos_5prime_in_msa:
    :param distance_to_upstream_nt:
    :return: the position of the upstream gene in the MSA; -1 if it isn't in the MSA
    """

    q_pos_5prime_in_msa = 0

    # skip all empty gaps
    while q_pos_5prime_in_msa < msa_t.alignment_length() and msa_t[index_in_msa][q_pos_5prime_in_msa] == "-":
        q_pos_5prime_in_msa += 1

    remaining = q_pos_5prime_in_msa_no_gap
    while remaining > 0:
        while q_pos_5prime_in_msa < msa_t.alignment_length() and msa_t[index_in_msa][q_pos_5prime_in_msa] == "-":
            q_pos_5prime_in_msa += 1

        q_pos_5prime_in_msa += 1
        remaining -= 1

    curr_pos = q_pos_5prime_in_msa
    distance_to_upstream_aa = int(distance_to_upstream_nt/3)

    if distance_to_upstream_aa > 0:         # no overlap

        remaining_aa = distance_to_upstream_aa
        while remaining_aa > 0:

            # move position upstream by 1 aa (i.e. ignoring all gaps)
            curr_pos -= 1
            while curr_pos >= 0 and msa_t[index_in_msa][curr_pos] == "-":
                curr_pos -= 1

            if curr_pos < 0:
                break

            remaining_aa -= 1

    elif distance_to_upstream_aa < 0:

        remaining_aa = -distance_to_upstream_aa

        while remaining_aa > 0:
            # move position downstream by 1 aa (ignoring all gaps)
            # move position upstream by 1 aa (i.e. ignoring all gaps)
            curr_pos += 1
            while curr_pos < msa_t.alignment_length() and msa_t[index_in_msa][curr_pos] == "-":
                curr_pos += 1

            if curr_pos >= msa_t.alignment_length():
                curr_pos = -1
                break

            remaining_aa -= 1

    return int(curr_pos)








    #
    # # type: (Dict[str, Any], pd.DataFrame, Dict[str, Any]) -> None
    #
    # suffix_coordinates = get_value(kwargs, "suffix_corrdinates", None)
    # suffix_upstream_distance = get_value(kwargs, "suffix_upstream_distance", "upstream-distance")
    #
    # # for each row
    # for index, row in df.iterrows():
    #     # get position of query in msa
    #     q_pos_of_5prime_in_msa = row["q-prot-pos-5prime-in-frag-msa"]
    #
    #     # get number of AA between query 5prime and upstream gene (zero or negative means overlapping)
    #     distance_aa = row["q-upstream-distance-msa"]
    #
    #     # find location of upstream gene in msa (-1 means it's not there)
    #
    #     pos_upstream_in_msa = move_by_and_ignore_gaps()


def perform_msa_on_ortholog_group(env, df_group, msa_options, **kwargs):
    # type: (Dict[str, Any], pd.DataFrame, MSAOptions, Dict[str, Any]) -> (pd.DataFrame, sbsp_general.labels.Label)

    order_by = get_value(kwargs, "order_by", None)
    column_k2p_distance = get_value(kwargs, "k2p_distance", "k2p-distance")
    ortholog_group_id = get_value(kwargs, "ortholog_group_id", 0)
    pf_msa_out = get_value(kwargs, "pf_msa_out", None)
    suffix_fname = get_value(kwargs, "suffix_fname", None)
    labels_info = get_value(kwargs, "labels_info", None)
    filter_stats = get_value(kwargs, "filter_stats", dict())

    column_q_msa_aa = "q-prot-msa"
    column_t_msa_aa = "t-prot-msa"

    column_q_msa_nt = "q-nucl-msa"
    column_t_msa_nt = "t-nucl-msa"

    column_q_pos_5_prime_msa = "q-prot-pos-5prime-in-frag-msa"
    column_t_pos_5_prime_msa = "t-prot-pos-5prime-in-frag-msa"

    if msa_options.safe_get("column-distance"):
        column_k2p_distance = msa_options.safe_get("column-distance")

    # if labels_info is None:
    #     labels_info = sbsp_general.dataframe.df_get_index_of_label_per_genome(env, df_group, "both")

    if order_by is not None:
        if order_by in df_group.columns.values:
            df_group.sort_values(order_by, inplace=True)
        else:
            logger.warning("Could not order MSA input by column ({})".format(order_by))

    list_sequences_aa = [df_group[column_q_msa_aa].iloc[0]] + list(df_group[column_t_msa_aa])
    list_sequences_nt = [df_group[column_q_msa_nt].iloc[0]] + list(df_group[column_t_msa_nt])

    list_distances = [0] + list(df_group[column_k2p_distance])
    if "kimura3" in df_group:
        list_distances = ["{}-{}".format(x, y) for x, y in zip(list_distances, [0] + list(df_group["kimura3"]))]


    if msa_options.safe_get("search-msa-max-sequence-length-aa") is not None:
        val = msa_options["search-msa-max-sequence-length-aa"]
        list_sequences_aa = [l[0:val] for l in list_sequences_aa]
        list_sequences_nt = [l[0:val*3] for l in list_sequences_nt]


    suffix_msa_fname = "{}".format(ortholog_group_id)
    if suffix_fname is not None:
        suffix_msa_fname = "{}-{}".format(suffix_msa_fname, suffix_fname)

    # run alignment
    [list_sequences_aa_aligned, stdout, stderr, alignments] = msa(
        list_sequences_aa,
        tags=list_distances,
        pd_work=env["pd-work"],
        suffix_fname=suffix_msa_fname
    )



    if msa_options["filter-by-pairwise-kimura-from-msa"]:
        filtered_list_sequeces_aa_aligned, indices_of_remaining = filter_by_pairwise_kimura_from_msa(
            list_sequences_aa_aligned, msa_options
        )


        filter_stats["filter-by-pairwise-kimura-from-msa"] = len(list_sequences_aa_aligned) - len(filtered_list_sequeces_aa_aligned)

        # if some sequences were filtered, rerun alignment
        if len(filtered_list_sequeces_aa_aligned) < len(list_sequences_aa_aligned):



            list_sequences_aa = get_elements_at_indices(list_sequences_aa, indices_of_remaining)
            list_distances = get_elements_at_indices(list_distances, indices_of_remaining)
            list_sequences_nt = get_elements_at_indices(list_sequences_nt, indices_of_remaining)

            # remove from dataframe
            df_group = df_group.iloc[[x-1 for x in indices_of_remaining if x != 0]]

            [list_sequences_aa_aligned, stdout, stderr, alignments] = msa(
                list_sequences_aa,
                tags=list_distances,
                pd_work=env["pd-work"],
                suffix_fname=suffix_msa_fname
            )

    if msa_options.safe_get("filter-remove-sequences-that-introduce-gaps"):

        before_filtering = len(list_sequences_aa_aligned)

        while True:
            alignments = lower_case_everything_except_starts(alignments, list_sequences_aa,
                                                             list_sequences_nt)

            list_sequences_aa_aligned = [a.seq._data for a in alignments]

            filtered_list_sequeces_aa_aligned, indices_of_remaining = filter_sequences_that_introduce_gaps_in_msa(
                list_sequences_aa_aligned, msa_options
            )

            # if nothing was filtered, stop
            if len(filtered_list_sequeces_aa_aligned) == len(list_sequences_aa_aligned):
                break

            list_sequences_aa = get_elements_at_indices(list_sequences_aa, indices_of_remaining)
            list_distances = get_elements_at_indices(list_distances, indices_of_remaining)
            list_sequences_nt = get_elements_at_indices(list_sequences_nt, indices_of_remaining)

            # remove from dataframe
            df_group = df_group.iloc[[x - 1 for x in indices_of_remaining if x != 0]]

            [list_sequences_aa_aligned, stdout, stderr, alignments] = msa(
                list_sequences_aa,
                tags=list_distances,
                pd_work=env["pd-work"],
                suffix_fname=suffix_msa_fname
            )


        filter_stats["filter-remove-sequences-that-introduce-gaps"] = before_filtering - len(
            filtered_list_sequeces_aa_aligned)

    # choose start
    # lower case everything, except
    alignments = lower_case_everything_except_starts(alignments, list_sequences_aa,
                                                     list_sequences_nt)

    try:
        alignments_nt = convert_msa_from_aa_to_nt(alignments, list_sequences_nt)
    except:
        pass
    list_sequences_nt_aligned = [a.seq._data for a in alignments_nt]


    score_on_all_pairs = msa_options.safe_get("search-score-on-all-pairs")

    upstream_block_scorer_name = get_value(kwargs, "search-upstream-block-scorer", None)
    upstream_block_scorer = None
    if upstream_block_scorer_name is not None:
        upstream_block_scorer = sbsp_ml.msa_features_2.ScoringMatrix(name=upstream_block_scorer_name)

    # # add distance to upstream gene
    # sbsp_general.dataframe.df_add_distance_to_upstream_gene_DEPRECATED(env, df_group, "q",
    #                                                                    suffix_coordinates=None,
    #                                                                    suffix_upstream_distance="upstream-distance-original",
    #                                                                    labels_info=labels_info)
    # sbsp_general.dataframe.df_add_distance_to_upstream_gene_DEPRECATED(env, df_group, "t",
    #                                                                    suffix_coordinates=None,
    #                                                                    suffix_upstream_distance="upstream-distance-original",
    #                                                                    labels_info=labels_info)

    # add distance to upstream gene
    sbsp_general.dataframe.df_add_distance_to_upstream_gene(env, df_group, "q",
                                                                       suffix_coordinates=None,
                                                                       suffix_upstream_distance="upstream-distance-original",
                                                                       labels_info=labels_info)
    sbsp_general.dataframe.df_add_distance_to_upstream_gene(env, df_group, "t",
                                                                       suffix_coordinates=None,
                                                                       suffix_upstream_distance="upstream-distance-original",
                                                                       labels_info=labels_info)

    sbsp_general.dataframe.df_add_position_of_upstream_gene_in_msa_no_gaps(df_group)

    # compute position of upstream gene in msa
    pos_of_upstream_in_msa = compute_position_of_upstream_of_query_in_msa(
        MSAType(alignments),
        int(df_group.iloc[0]["q-prot-pos-5prime-in-frag-msa"]),
        distance_to_upstream_nt=int(df_group.iloc[0]["q-upstream-distance-original"])
    )
    df_group["q-upstream-pos-in-frag-msa"] = pos_of_upstream_in_msa


    thresh = msa_options["search-start-selection-threshold"]
    start_position_in_msa = None
    if "search-select-lorf" in msa_options and msa_options["search-select-lorf"]:
        start_position_in_msa = select_start_position_from_msa_for_lorf(alignments,
                                                                        score_on_all_pairs=score_on_all_pairs,
                                                                        scorer=upstream_block_scorer,
                                                                        msa_options=msa_options,
                                                                        pos_of_upstream_in_msa=pos_of_upstream_in_msa)

    if start_position_in_msa is None:
        params = msa_options.safe_get("search-candidates-without-upstream-conservation")
        region_length = threshold = None
        if params is not None:
            region_length, threshold = params

        candidate_positions = get_candidates_without_upstream_conservation(alignments, region_length=region_length,
                                                                           threshold=threshold,
                                                                           score_on_all_pairs=score_on_all_pairs,
                                                                           scorer=upstream_block_scorer,
                                                                           msa_options=msa_options,
                                                                           pos_of_upstream_in_msa=pos_of_upstream_in_msa)
        if msa_options.safe_get("search-candidates-before-conserved-block"):
            max_fraction_of_gaps_in_pos = msa_options.safe_get("search-max-fraction-of-gaps-in-pos")
            region_length = threshold = None
            if msa_options.safe_get("search-candidates-before-conserved-block"):
                region_length, threshold = msa_options.safe_get("search-candidates-before-conserved-block")
            candidate_positions = get_candidates_before_any_conservation_block(alignments,
                                                                               max_fraction_of_gaps_in_pos=max_fraction_of_gaps_in_pos,
                                                                               region_length=region_length,
                                                                               threshold=threshold)

        if len(candidate_positions) == 0:
            logger.debug("No filtered candidates. Run search on entire gene")
            start_position_in_msa = select_start_position_from_msa(alignments, threshold=thresh, msa_options=msa_options,
                                                                   score_on_all_pairs=score_on_all_pairs,
                                                                   upstream_block_scorer=upstream_block_scorer,
                                                                   pos_of_upstream_in_msa=pos_of_upstream_in_msa)
        else:
            logger.debug("Run search for candidate positions: {}".format(candidate_positions))
            start_position_in_msa = select_start_position_from_msa(alignments, threshold=thresh,
                                                                   msa_options=msa_options,
                                                                   end_at_pos=candidate_positions[-1]+1,
                                                                   score_on_all_pairs=score_on_all_pairs,
                                                                   upstream_block_scorer=upstream_block_scorer,
                                                                   pos_of_upstream_in_msa=pos_of_upstream_in_msa)

            # FIXME: remove from else block
            if start_position_in_msa is None and msa_options.safe_get("search-skip-by-standard-aa-score"):
                logger.debug("No candidate found. Select from right, but skip with high standard aa score")
                start_position_in_msa = find_rightmost_by_standard_aa_score(
                    alignments,
                    candidate_positions,
                    threshold=msa_options.safe_get("search-skip-by-standard-aa-score"),
                )

            if start_position_in_msa is None and msa_options.safe_get("search-select-closest-to-coding-on-fail"):
                logger.debug("No candidate found. Select rightmost")
                start_position_in_msa = candidate_positions[-1]



    q_pos_of_previous_start_in_msa = df_group[column_q_pos_5_prime_msa].iloc[0]

    # count skipped MSA where no start position is above threshold
    if start_position_in_msa is None:
        if filter_stats is not None:
            key = "filter-no-start-position-passes-threshold"
            if key not in filter_stats:
                filter_stats[key] = 0
            filter_stats[key] += 1

    # if a start is found
    q_label_msa = None
    if start_position_in_msa is not None:

        # get label of new start
        q_label_msa = get_label_from_start_position_in_msa(
            df_group,
            list_sequences_aa[0],
            list_sequences_aa_aligned[0],
            q_pos_of_previous_start_in_msa,
            start_position_in_msa)

        # add to dataframe
        suffix_msa_updated = "msa"
        for index, row in df_group.iterrows():
            sbsp_general.dataframe.df_coordinates_to_row(df_group, index, q_label_msa, "q", suffix_msa_updated)


        # add distance to upstream gene
        sbsp_general.dataframe.df_add_distance_to_upstream_gene_DEPRECATED(env, df_group, "q",
                                                                           suffix_coordinates=suffix_msa_updated,
                                                                           suffix_upstream_distance="upstream-distance-msa",
                                                                           labels_info=labels_info)
        sbsp_general.dataframe.df_add_distance_to_upstream_gene_DEPRECATED(env, df_group, "t",
                                                                           suffix_coordinates=suffix_msa_updated,
                                                                           suffix_upstream_distance="upstream-distance-msa",
                                                                           labels_info=labels_info)

        # add position of upstream gene in msa


        # add gene lengths
        target_number = 1
        for index, row in df_group.iterrows():
            # print row
            t_label_msa = get_label_from_start_position_in_msa(
                df_group,
                list_sequences_aa[target_number],
                list_sequences_aa_aligned[target_number],
                df_group[column_t_pos_5_prime_msa].iloc[target_number-1],       # -1 to account for added query sequence
                start_position_in_msa,
                source="t",
                row_num=target_number-1
            )

            sbsp_general.dataframe.df_coordinates_to_row(df_group, index, t_label_msa, "t", suffix_msa_updated)

            df_group.at[index, "q-gene-length-msa"] = q_label_msa.length()
            df_group.at[index, "t-gene-length-msa"] = t_label_msa.length()

            target_number += 1

        # compute 5prime conservation scores
        compute_post_alignment_features(
                df_group,
                list_sequences_aa_aligned,
                list_sequences_nt_aligned,
                start_position_in_msa,
                start_position_in_msa*3,
                msa_options,
                10
        )


        # FILTERING
        while True:
            num_orthologs_before_filtering = len(df_group)
            df_group = filter_df_by_upstream_to_msa(df_group, msa_options, suffix_upstream="upstream-distance-msa",
                                                    filter_stats=filter_stats)
            num_orthologs_after_filtering = len(df_group)

            if num_orthologs_after_filtering == 0:
                return df_group, None

            if msa_options["filter-min-number-orthologs"] is not None and \
                    num_orthologs_after_filtering < msa_options["filter-min-number-orthologs"]:
                if "filter-min-number-orthologs" not in filter_stats:
                    filter_stats["filter-min-number-orthologs"] = 0
                filter_stats["filter-min-number-orthologs"] += 1
                return df_group, None

            # if no filtering is done, we're done :)
            if num_orthologs_after_filtering == num_orthologs_before_filtering:
                break

            # otherwise, rerun msa
            df_group, q_label_msa = perform_msa_on_ortholog_group(env, df_group, msa_options, **kwargs)
            pf_msa_out = None


        # format alignment and print to file
        if pf_msa_out is not None:
            print_alignment_to_file(alignments, start_position_in_msa, df_group, pf_msa_out, pos_of_upstream_in_msa=pos_of_upstream_in_msa)


            print_alignment_to_file(alignments_nt, start_position_in_msa*3, None,
                                    "{}_nt".format(pf_msa_out), pos_of_upstream_in_msa=pos_of_upstream_in_msa*3)
    # otherwise no start is found - print the alignment to file
    else:
        if pf_msa_out is not None:
            pf_msa_out = "{}_nopred".format(pf_msa_out)
            print_alignment_to_file(alignments, None, df_group, pf_msa_out)

            print_alignment_to_file(alignments_nt, None, None,
                                    "{}_nt".format(pf_msa_out))

    i = 1
    for index, _ in df_group.iterrows():
        df_group.at[index, "t-pos-in-msa-output"] = i
        i += 1

    def set_upstream_position_for_each_target_sequence(curr_df_group, alignments):
        # type: (pd.DataFrame, MSAType) -> None

        for anindex, arow in curr_df_group.iterrows():
            seq_index_in_msa = int(arow["t-pos-in-msa-output"])
            curr_df_group.at[anindex, "t-upstream-pos-in-frag-msa"] = compute_position_of_upstream_of_query_in_msa(
                alignments,
                int(arow["t-prot-pos-5prime-in-frag-msa"]),
                distance_to_upstream_nt=int(arow["t-upstream-distance-original"]),
                index_in_msa=seq_index_in_msa
            )

    set_upstream_position_for_each_target_sequence(df_group, MSAType(alignments))

    return df_group, q_label_msa


def perform_msa_on_df_and_compute_conservation_features(env, df, **kwargs):
    # type: (Dict[str, Any], pd.DataFrame, Dict[str, Any]) -> pd.DataFrame
    pass

def print_filter_stats(meine_filter_stats):
    for k in sorted(meine_filter_stats.keys()):
        logger.critical(
            "Filter ({}): {}".format(k, meine_filter_stats[k]))

def print_filter_stats_to_file(meine_filter_stats, pf_filter_stats):

    output = ""
    for k in sorted(meine_filter_stats.keys()):
        output += "Filter ({}): {}".format(k, meine_filter_stats[k]) + "\n"

    from sbsp_io.general import write_string_to_file
    write_string_to_file(output, pf_filter_stats)


def perform_msa_on_df(env, df, **kwargs):
    # type: (Dict[str, Any], pd.DataFrame, Dict[str, Any]) -> pd.DataFrame

    msa_options = get_value(kwargs, "msa_options", MSAOptions(env))
    msa_output_start = get_value(kwargs, "msa_output_start", 0)
    dn_msa_output = get_value(kwargs, "dn_msa_output", None)
    column_k2p_distance = get_value(kwargs, "k2p_distance", "k2p-distance")
    suffix_coordinates = get_value(kwargs, "suffix_coordinates", None)
    tag_msa = get_value(kwargs, "tag_msa", "msa")

    upstream_length_nt = get_value(kwargs, "upstream_length_nt", None)
    downstream_length_nt = get_value(kwargs, "downstream_length_nt", None)

    logger.info("Performing msa on df with parameters:\n{}".format(msa_options.to_string()))

    if msa_options.safe_get("column-distance"):
        column_k2p_distance = msa_options.safe_get("column-distance")

    # setup directory
    pd_msa_output = setup_directory_for_msa_outputs(env, dn_msa_output)

    filter_stats = {x: 0 for x in msa_options.option_names() if x.startswith("filter")}

    # sbsp_general.general.df_add_5prime_3prime_key(df, "q-")
    #
    # Step 1: Filter data points
    logger.debug("Pipeline Step 1: Filter")
    df = filter_df(df, msa_options, column_distance=column_k2p_distance, filter_stats=filter_stats, filter_non_group_only=False )

    limit_upstream_to_first_candidate = msa_options.safe_get("search-limit-upstream-to-first-candidate")


    # Step 2: Add sequences for each data point (to perform alignment on)
    logger.debug("Pipeline Step 2: Add sequences for each genome")
    df = sbsp_general.dataframe.df_add_labeled_sequences(env, df,
                                                         source="both",
                                                         suffix_coordinates=suffix_coordinates,
                                                         suffix_gene_sequence=tag_msa,
                                                         upstream_length_nt=upstream_length_nt,
                                                         downstream_length_nt=downstream_length_nt,
                                                         limit_upstream_to_first_candidate=limit_upstream_to_first_candidate)

    # Step 3: Perform MSA on each ortholog group
    logger.debug("Pipeline Step 3: Perform MSA on each ortholog group")

    # get information for labels to use in computation of upstream distances
    labels_info = None # sbsp_general.dataframe.df_get_index_of_label_per_genome(env, df, "both")

    df_result = pd.DataFrame()

    if column_k2p_distance not in df.columns.values:
        df[column_k2p_distance] = df["evalue"]

    # get all file labels in order to compute

    # for each ortholog group
    msa_number = 0
    for name, df_group in df.groupby("q-5prime-3prime"):

        # if name == "Escherichia_coli_K_12_substr__MG1655_uid57779;NC_000913;2080789;2081262;-":
        #     pass
        # else:
        #     continue
        if len(df_group)+1 < msa_options["filter-min-number-orthologs"]:
            filter_stats["filter-min-number-orthologs"] += 1
            continue

        df_group.sort_values(by=[column_k2p_distance], inplace=True)

        # print (pd_msa_output, msa_output_start, msa_number)
        pf_msa_out = os.path.join(pd_msa_output, "msa_{}_{}.out".format(msa_output_start, msa_number))
        [df_group_msa, q_label] = perform_msa_on_ortholog_group(env, df_group, msa_options, pf_msa_out=pf_msa_out,
                                                                labels_info=labels_info,
                                                                ortholog_group_id=msa_number,
                                                                suffix_fname=msa_output_start,
                                                                filter_stats=filter_stats)

        if q_label is not None:
            df_result = df_result.append(df_group_msa)
            msa_number += 1

    # update msa outputs
    # sbsp_general.dataframe.df_add_distance_to_upstream_gene(env, df_result, "q", suffix_coordinates=tag_msa, suffix_upstream_distance="upstream-distance-msa")
    # sbsp_general.dataframe.df_add_distance_to_upstream_gene(env, df_result, "t", suffix_coordinates=tag_msa, suffix_upstream_distance="upstream-distance-msa")



    print_filter_stats(filter_stats)

    if len(df_result) == 0:
        return df_result

    import sbsp_io.msa_2
    sbsp_io.msa_2.update_msa_outputs(df_result,
                                              column_q_upstream_distance="q-upstream-distance-msa",
                                              column_t_upstream_distance="t-upstream-distance-msa")

    df_result["q-left-msa"] = df_result["q-left-msa"].astype(int)
    df_result["q-right-msa"] = df_result["q-right-msa"].astype(int)

    return df_result



def select_start_for_msa_from_file(env, pf_msa, **kwargs):
    # type: (Dict[str, Any], str, Dict[str, Any]) -> None

    suffix = get_value(kwargs, "suffix", None)
    msa_options = get_value(kwargs, "msa_options", MSAOptions(env))

    # read alignment





