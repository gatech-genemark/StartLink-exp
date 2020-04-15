import random
import logging

from typing import *

from sbsp_general.general import get_value
from sbsp_general.labels import Label

logger = logging.getLogger(__name__)


def get_possible_positions_for_candidate_starts_in_non_coding(sequence, label, upstream_length_nt):
    # type: (str, Label, int) -> List[Tuple[int, str]]

    list_positions = list()

    if label.strand() == "+":

        left = label.left()

        curr_pos = left - 3
        min_pos = max(left - upstream_length_nt, 0)

        while curr_pos >= min_pos:

            codon = sequence[left:left+3]

            # if one substitution away from ATG, GTG, or TTG
            if codon[:2] in {"AT", "TT", "GT"}:
                list_positions.append((curr_pos, codon[:2] + "G"))
            elif codon[1:3] == "TG":
                list_positions.append((curr_pos, "ATG"))
            curr_pos -= 3

    else:
        raise NotImplementedError("Only positive strand supported for now")

    return list_positions


def get_possible_positions_for_candidate_starts_in_coding(sequence, label, downstream_length_nt):
    # type: (str, Label, int) -> List[Tuple[int, str]]

    list_positions = list()

    if label.strand() == "+":

        left = label.left()

        curr_pos = left + 3
        max_pos = min(left + downstream_length_nt, label.right())

        while curr_pos <= max_pos:

            codon = sequence[left:left+3]

            # synonymous Valine
            if codon in {"GTC", "GTT", "GTA"}:
                list_positions.append((curr_pos, "GTG"))
            # synonymous Leucine
            elif codon in {"TTA", "CTA", "CTT", "CTC", "CTG"}:
                list_positions.append((curr_pos, "TTG"))
            # something close to methionine (Isoleucine)
            elif codon in {"ATT", "ATC", "ATA"}:
                list_positions.append((curr_pos, "ATG"))

            curr_pos += 3
    else:
        raise NotImplementedError("Only positive strand supported for now")

    return list_positions


def add_candidate_starts(sequence, label, **kwargs):
    # type: (str, Label, Dict[str, Any]) -> str

    list_seq = list(sequence)

    num_starts_upstream = get_value(kwargs, "num_starts_upstream", 0, default_if_none=True)
    num_starts_downstream = get_value(kwargs, "num_starts_downstream", 0, default_if_none=True)
    upstream_length_nt = get_value(kwargs, "upstream_length_nt", default=180, value_type=int, default_if_none=True)
    downstream_length_nt = get_value(kwargs, "downstream_length_nt", default=180, value_type=int, default_if_none=True)

    def update_sequence_with_codons(list_sequence, list_position_codon_pairs):
        # type: (List[str], List[Tuple[int, str]]) -> None
        for pc_pair in list_position_codon_pairs:
            pos, codon = pc_pair
            list_sequence[pos:pos+len(codon)] = list(codon)

    if num_starts_upstream > 0:
        list_positions = get_possible_positions_for_candidate_starts_in_non_coding(sequence, label, upstream_length_nt)

        if num_starts_upstream < len(list_positions):
            list_positions = random.sample(list_positions, num_starts_upstream)

        update_sequence_with_codons(list_seq, list_positions)

    if num_starts_downstream > 0:
        list_positions = get_possible_positions_for_candidate_starts_in_non_coding(sequence, label,
                                                                                   downstream_length_nt)

        if num_starts_downstream < len(list_positions):
            list_positions = random.sample(list_positions, num_starts_upstream)

        update_sequence_with_codons(list_seq, list_positions)

    return "".join(list_seq)

