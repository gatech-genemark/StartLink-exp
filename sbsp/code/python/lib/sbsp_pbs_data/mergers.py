import math
import os
import logging
from typing import *

from sbsp_containers.genome_list import GenomeInfoList
from sbsp_general import Environment
from sbsp_general.general import get_value
from sbsp_options.pipeline_msa import PipelineMSAOptions

log = logging.getLogger(__name__)

T = TypeVar('T')


def merge_identity(list_output_data):
    # type: (List[T]) -> List[T]
    return list_output_data
