import logging
from typing import *

import Bio.Alphabet

import sbsp_alg.msa

log = logging.getLogger(__name__)


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