from typing import Any

import pandas as pd
import numpy as np
import re
from typing import *
import copy

import sbsp_alg.msa
from sbsp_general.msa_2 import MSAType, MultipleSeqAlignment
from sbsp_options.sbsp import SBSPOptions
from sbsp_general.general import get_value
from sbsp_io.general import read_rows_to_list
import Bio.SubsMat.MatrixInfo


class ScoringMatrix:

    def __init__(self, name="identity", ignore_case=True):
        self._name = name
        self._ignore_case = ignore_case

        if name.startswith("blosum"):
            self._blosum_matrix = getattr(Bio.SubsMat.MatrixInfo, name)
            ScoringMatrix.extend_blosum(self._blosum_matrix)
            self._scorer = self._score_blosum

        elif name == "identity":
            self._scorer = self._score_identity

        else:
            raise ValueError("Unknown scoring function")

    @staticmethod
    def extend_blosum(matrix, stop_to_aa=-4, stop_to_stop=1, gap_and_aa=-4, gap_and_gap=-1):
        # type: (Dict[(str, str), float], int, int, int, int) -> None
        """
        Extends blosum by making it symmetric (why the hell isn't it!!), and adding mapping to
        stop codons (i.e. *)
        :param matrix:
        :param stop_to_aa:
        :param stop_to_stop:
        :return:
        """
        # get unique aa
        unique_aa = set(item for sublist in matrix for item in sublist)     # get unique set of AA

        copy_matrix = copy.deepcopy(matrix)

        # make symmetric
        for a in copy_matrix:
            matrix[(a[1], a[0])] = copy_matrix[a]

        # for each, add penalty score
        for aa in unique_aa:
            matrix[(aa, "*")] = stop_to_aa
            matrix[("*", aa)] = stop_to_aa

            matrix[(aa, "-")] = gap_and_aa
            matrix[("-", aa)] = gap_and_aa

        matrix[("*", "*")] = stop_to_stop
        matrix[("-", "-")] = gap_and_gap


    def _score_identity(self, a, b):
        # type: (str, str) -> int
        return 1 if a == b and a != "-" else 0

    def _score_blosum(self, a, b):
        # type: (str, str) -> float
        return self._blosum_matrix[(a, b)]

    def score(self, a, b):
        # type: (str, str) -> float
        a = a.upper() if self._ignore_case else a
        b = b.upper() if self._ignore_case else b
        return self._scorer(a, b)


def compute_upstream_score(msa_t, position, msa_options, **kwargs):
    # type: (MSAType, int, SBSPOptions, Dict[str, Any]) -> float

    require_full_length = get_value(kwargs, "require_full_length", False)
    ignore_gaps_in_query = get_value(kwargs, "ignore_gaps_in_query", False)
    score_on_all_pairs = get_value(kwargs, "score_on_all_pairs", False)

    scoring_function = get_value(kwargs, "scoring_function", ScoringMatrix("identity"), default_if_none=True)

    region_length = msa_options["search-upstream-of-conserved-region"]

    begin = position - region_length        # inclusive
    end = position                          # exclusive (don't count start)

    if begin < 0:
        if require_full_length:
            raise ValueError("Not enough upstream region")
        begin = 0

    score = sbsp_alg.msa.compute_conservation_in_region(
        [x.seq._data for x in msa_t.list_alignment_sequences],          # TODO: make compatible
        begin,
        end,
        skip_gaps=ignore_gaps_in_query,
        only_full_length=require_full_length,
        direction="upstream",
        scorer=scoring_function,
        score_on_all_pairs=score_on_all_pairs
    )

    return score


def compute_downstream_score(msa_t, position, msa_options, **kwargs):
    # type: (MSAType, int, SBSPOptions, Dict[str, Any]) -> float

    require_full_length = get_value(kwargs, "require_full_length", False)
    ignore_gaps_in_query = get_value(kwargs, "ignore_gaps_in_query", False)

    scoring_function = get_value(kwargs, "scoring_function", ScoringMatrix("identity"), default_if_none=True)

    region_length = msa_options["search-upstream-of-conserved-region"]

    begin = position + 1                # inclusive (don't count start)
    end = position + 1 + region_length  # exclusive

    if end >= msa_t.alignment_length():
        if require_full_length:
            raise ValueError("Not enough downstream region")
        end = msa_t.alignment_length()

    score = sbsp_alg.msa.compute_conservation_in_region(
        [x.seq._data for x in msa_t.list_alignment_sequences],          # TODO: make compatible
        begin,
        end,
        skip_gaps=ignore_gaps_in_query,
        only_full_length=require_full_length,
        direction="downstream",
        scorer=scoring_function
    )

    return score


def compute_simple_saas(msa_t, i):
    # type: (MSAType, int) -> float

    num_vli = sum(1 for j in range(msa_t.number_of_sequences()) if msa_t[j][i] in {"v", "l", "i", "-"})

    return float(num_vli) / msa_t.number_of_sequences()

def compute_5prime_score(msa_t, position, msa_options, **kwargs):
    # type: (MSAType, int, SBSPOptions, Dict[str, Any]) -> float

    num_upper = sbsp_alg.msa.count_num_upper(msa_t.list_alignment_sequences, position, msa_options)
    start_identity = num_upper / float(msa_t.number_of_sequences())

    return start_identity


def f_stop_after(msa, pos):
    # type: (MSAType, int) -> float

    total = msa.number_of_sequences()

    def count_non_gap_up_to(msa, seq_num, pos_inclusive):
        # type: (MSAType, int, int) -> int
        return (pos_inclusive+1) - msa[seq_num][:pos_inclusive+1].seq._data.count("-")

    count_stop_after = 0
    for i in range(total):

        num_non_gap = count_non_gap_up_to(msa, i, pos)

        if num_non_gap == 0:
            count_stop_after += 1

    return count_stop_after / float(total)


def is_lorf(msa_t, pos):
    # type: (MSAType, int) -> bool

    lorf_status = True

    curr_pos = pos - 1
    while curr_pos >= 0:
        if msa_t[0][curr_pos].isupper():
            lorf_status = False
            break
        curr_pos -= 1

    return lorf_status










def compute_diversity_score(alignments, pos, radius=0, elements=None, elements_count=None, elements_ignore=None,
                            limit_radius_by_next_methianine=False):
    # type: (List[SeqRecord], int, int, Set[str], Set[str], Set[str], bool) -> float

    """

    :param alignments:
    :param pos:
    :param radius:
    :param elements:
    :param elements_count:
    :param elements_ignore:
    :param limit_radius_by_next_methianine: If set, the radius is limited by the next "methianine" in the query
    :return:
    """
    import Bio.Alphabet.IUPAC

    if elements is None:
        elements = set(Bio.Alphabet.IUPAC.IUPACProtein().letters)

    if elements_count is None:
        elements_count = elements

    if elements_ignore is None:
        elements_ignore = elements - elements_count

    num_sequences = len(alignments)
    alignment_length = len(alignments[0])

    denom = num_sequences * (num_sequences - 1)

    counts = dict()

    # for each sequence, check if element should be counted exists
    for n in range(num_sequences):

        letter = alignments[n][pos]

        count_it = False
        letter_to_count = None

        if letter in elements_count:
            count_it = True
            letter_to_count = letter
        else:
            stop_looking_upstream = False
            stop_looking_downstream = False
            for j in range(1, radius):
                query_letter_at_j_downstream = None
                query_letter_at_j_upstream = None

                if pos + j < alignment_length:
                    query_letter_at_j_downstream = alignments[0][pos + j]
                if pos - j >= 0:
                    query_letter_at_j_upstream = alignments[0][pos - j]


                if limit_radius_by_next_methianine:

                    if query_letter_at_j_downstream is not None and query_letter_at_j_downstream.isupper():
                        stop_looking_downstream = True

                    if query_letter_at_j_upstream is not None and query_letter_at_j_upstream.isupper():
                        stop_looking_upstream = True


                # check upstream
                if not stop_looking_upstream and query_letter_at_j_upstream is not None:
                    if query_letter_at_j_upstream in elements_count:
                        count_it = True
                        letter_to_count = query_letter_at_j_upstream
                        break

                # check downstream
                if not stop_looking_downstream and query_letter_at_j_downstream is not None:
                    if query_letter_at_j_downstream in elements_count:
                        count_it = True
                        letter_to_count = query_letter_at_j_downstream
                        break

        if count_it:
            if letter_to_count not in counts.keys():
                counts[letter_to_count] = 0

            counts[letter_to_count] += 1


    # compute score

    if len(counts) == 0:
        return 0.0

    score = 1.0 / float(denom)

    num = 0
    for l in counts:
        num += counts[l] * (counts[l] - 1)

    score *= num

    return score


def compute_saas(alignments, pos, upstream_region=0, downstream_region=0):
    # type: (List[SeqRecord], int) -> float

    total = len(alignments)

    list_aa = ["M", "V", "L", "v", "l", "other"]

    letters = sbsp_alg.msa.get_aligned_letters_at_position(alignments, pos, upstream_region, downstream_region)

    # get counts at position
    counts = {x: 0 for x in list_aa}
    for letter in letters:
        if letter in counts.keys():
            counts[letter] += 1
        else:
            counts["other"] += 1

    numerator = 0
    numerator += counts["V"] * counts["v"]
    numerator += counts["L"] * counts["l"]
    # numerator -= counts["M"] * counts["V"]
    # numerator -= counts["M"] * counts["L"]
    # numerator -= counts["V"] * counts["L"]
    # numerator -= counts["M"] * counts["V"] * counts["L"]

    # additional numberator
    # numerator += 3 * pow(total / 3.0, 2) + pow(total / 3.0, 3)

    denominator = 3 * pow(total / 3.0, 2) #+ pow(total / 3.0, 3) + pow(total / 2.0, 2)

    return numerator / denominator



def compute_mds(alignments, pos, upstream_region=0, downstream_region=0):
    # type: (List[SeqRecord], int) -> float

    list_aa = ["M", "V", "L", "other"]

    letters = sbsp_alg.msa.get_aligned_letters_at_position(alignments, pos, upstream_region, downstream_region)
    # get counts at position
    # get counts at position
    counts = {x: 0 for x in list_aa}
    for letter in letters:
        if letter in counts.keys():
            counts[letter] += 1
        else:
            counts["other"] += 1

    total = len(alignments)

    numerator = 0
    numerator += counts["M"] * counts["V"]
    numerator += counts["M"] * counts["L"]
    numerator += counts["V"] * counts["L"]
    numerator += counts["M"] * counts["V"] * counts["L"]

    denominator = 3 * pow(total / 3.0, 2) + pow(total / 3.0, 3)

    return numerator / denominator

