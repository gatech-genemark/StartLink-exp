import logging
from typing import *

from sbsp_general.general import get_value
from sbsp_container.msa import MSAType

logger = logging.getLogger(__name__)


def get_positions_of_query_candidate_starts_in_msa(msa_t):
    # type: (MSAType) -> List[int]

    list_positions = list()

    for pos in range(msa_t.alignment_length()):

        if msa_t[0][pos].isupper():
            list_positions.append(pos)

    return list_positions


def get_positions_of_query_candidate_starts_in_msa_with_class(msa_t, **kwargs):
    # type: (MSAType, Dict[str, Any]) -> [List[int], List[str]]

    reference_marker_name = get_value(kwargs, "reference_marker_name", "ref", default_if_none=True)

    max_number_downstream = get_value(kwargs, "max_number_downstream", None)    # if set, select closest to reference
    max_number_upstream = get_value(kwargs, "max_number_upstream", None)        # if set, select closest to reference

    marker = msa_t.get_marker(reference_marker_name)
    if marker is None:
        raise ValueError("Invalid marker with name: {}".format(reference_marker_name))

    ref_position = marker.mark_position

    list_positions = list()
    list_class = list()

    num_upstream = 0
    num_downstream = 0

    def get_class(reference_position, current_position):
        # type: (int, int) -> str
        if current_position < reference_position:
            return "Non-coding"
        if current_position == reference_position:
            return "True Start"
        if current_position > reference_position:
            return "Coding"

    for pos in range(msa_t.alignment_length()):

        if msa_t[0][pos].isupper():

            curr_class = get_class(ref_position, pos)

            # early stopping for downstream
            if max_number_downstream is not None:
                if  curr_class == "Coding" and num_downstream >= max_number_downstream:
                    break

            list_positions.append(pos)
            list_class.append(curr_class)

            if curr_class == "Coding":
                num_downstream += 1
            if curr_class == "Non-coding":
                num_upstream += 1

    # filter upstream if set
    if max_number_upstream is not None:
        if num_upstream > max_number_upstream:
            list_positions = list_positions[num_upstream-max_number_upstream:]
            list_class = list_class[num_upstream-max_number_upstream:]

    return list_positions, list_class
