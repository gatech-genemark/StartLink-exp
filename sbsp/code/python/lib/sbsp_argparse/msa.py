import argparse
from typing import *


def add_msa_options(parser):
    # type: (argparse.ArgumentParser) -> None

    parser.add_argument('--pf-msa-options', required=False, default=None, type=Union[str])

    parser.add_argument("--filter-max-distance", type=float, default=None)
    parser.add_argument("--filter-min-distance", type=float, default=None)
    parser.add_argument("--filter-min-number-orthologs", default=None, type=Union[int])
    parser.add_argument("--filter-max-number-orthologs", default=None, type=Union[int])
    parser.add_argument("--filter-orthologs-with-equal-kimura", default=None, type=Union[int],
                        help="If set, orthologs with equal distance (up to the specified decimal place)"
                             "will only have one of them kept")
    parser.add_argument("--filter-by-pairwise-kimura-from-msa", type=float, default=None)
    parser.add_argument("--filter-min-upstream-distance", type=Union[int], default=None)
    parser.add_argument("--filter-close-max-upstream-distance", type=Union[int], default=None)
    parser.add_argument("--filter-close-far-upstream-distance", type=Union[int], default=None)
    parser.add_argument("--filter-gene-length-max-percent-difference", type=float, default=None)

    # search
    parser.add_argument("--search-start-selection-threshold", default=None, type=float)
    parser.add_argument("--search-better-downstream-aa", default=None, type=Union[int])
    parser.add_argument("--search-favor-m", default=None, type=Union[int])
    parser.add_argument("--search-ignore-m-to-l-mutation", default=None, type=Union[bool])
    parser.add_argument("--search-ignore-l-to-m-mutation", default=None, type=Union[bool])
    parser.add_argument("--search-neighbor", default=None, type=Union[int])
    parser.add_argument("--search-msa-max-sequence-length-aa", default=None, type=Union[int])
    parser.add_argument("--search-upstream-of-conserved-region", default=None, type=Union[int])
    parser.add_argument("--search-upstream-of-conserved-region-threshold", default=None, type=float)
