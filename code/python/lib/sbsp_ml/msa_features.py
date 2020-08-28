from typing import *
import copy

import sbsp_alg.msa
from sbsp_container.msa import MSAType
from sbsp_options.sbsp import SBSPOptions
from sbsp_general.general import get_value
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










