import argparse
from typing import *


def add_processor_parallelization_options(parser):
    # type: (argparse.ArgumentParser) -> None

    parser.add_argument('--num-processors', required=False, type=int, default=None, help="Number of processors")


# pbs options
def add_pbs_options(parser):
    # type: (argparse.ArgumentParser) -> None

    parser.add_argument('--pf-pbs-options', required=False, default=None, help="Configuration file for PBS options")

    parser.add_argument('--use-pbs', required=False, default=None, action='store_true', help="Use PBS job scheduler")
    parser.add_argument('--num-pbs-jobs', required=False, type=Union[int], default=None, help="Number of PBS jobs")
    add_processor_parallelization_options(parser)
    parser.add_argument('--pd-work-pbs-node', required=False, default=None,
                        help="Path to working directory on PBS node")

