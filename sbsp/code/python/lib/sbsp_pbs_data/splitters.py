import math
import os
import logging
from typing import *

from sbsp_containers.genome_list import GenomeInfoList
from sbsp_general import Environment
from sbsp_general.general import get_value

log = logging.getLogger(__name__)

T = TypeVar('T')

def get_number_elements_per_split(total_elements, num_splits):
    # type: (int, int) -> int
    return math.ceil(total_elements / float(num_splits))


def generate_splits(elements, num_splits):
    # type: (Sized[T], int) -> Generator[T]

    num_per_split = get_number_elements_per_split(len(elements), num_splits)

    num_elements = len(elements)
    start = 0

    while start < num_elements:

        end = min(num_elements, start + num_per_split)
        yield elements[start:end]
        start = end


def split_query_genomes_target_genomes_one_vs_group(data, num_splits, pd_work, **kwargs):
    # type: (Dict[str, str], int, str, Dict[str, Any]) -> List[Dict[str, str]]

    fn_q_split_formatted = get_value(kwargs, "fn_q_split_formatted", "q_split_{}.list")
    fn_t_split_formatted = get_value(kwargs, "fn_t_split_formatted", "t_split_{}.list")

    pf_q_list = data["pf-q-list"]
    pf_t_list = data["pf-t-list"]

    list_pf_splits = list()

    q_list = GenomeInfoList.init_from_file(pf_q_list)
    t_list = GenomeInfoList.init_from_file(pf_t_list)

    split_number = 0
    for q_genome in q_list:

        # split
        for split in generate_splits(t_list, num_splits):
            pf_q_split = os.path.join(pd_work, fn_q_split_formatted.format(split_number))
            pf_t_split = os.path.join(pd_work, fn_t_split_formatted.format(split_number))

            q_split = GenomeInfoList([q_genome])
            t_split = GenomeInfoList(split)

            q_split.to_file(pf_q_split)
            t_split.to_file(pf_t_split)

            list_pf_splits.append({
                "pf-q-list": pf_q_split,
                "pf-t-list": pf_t_split
            })

    return list_pf_splits
