# Karl Gemayel
# Georgia Institute of Technology
#
# Created: 3/10/20

import logging
import argparse
from typing import *

# noinspection All
from Bio.Seq import Seq

import pathmagic

# noinspection PyUnresolvedReferences
import sbsp_log  # runs init in sbsp_log and configures logger

# Custom imports
from sbsp_general import Environment

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from sbsp_general.labels import Label, Coordinates
from sbsp_io.sequences import read_fasta_into_hash

parser = argparse.ArgumentParser("Description of driver.")

parser.add_argument('--pf-sequence', required=True)
parser.add_argument('--accession', type=str, required=True)
parser.add_argument('--left', type=int, required=True)
parser.add_argument('--right', type=int, required=True)
parser.add_argument('--strand', type=str, required=True)


parser.add_argument('--pd-work', required=False, default=None, help="Path to working directory")
parser.add_argument('--pd-data', required=False, default=None, help="Path to data directory")
parser.add_argument('--pd-results', required=False, default=None, help="Path to results directory")
parser.add_argument("-l", "--log", dest="loglevel", choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                    help="Set the logging level", default='WARNING')

parsed_args = parser.parse_args()

# ------------------------------ #
#           Main Code            #
# ------------------------------ #

# Load environment variables
my_env = Environment(pd_data=parsed_args.pd_data,
                     pd_work=parsed_args.pd_work,
                     pd_results=parsed_args.pd_results)

# Setup logger
logging.basicConfig(level=parsed_args.loglevel)
logger = logging.getLogger("logger")  # type: logging.Logger


def get_fragment(sequences, label):
    # type: (Dict[str, Seq], Label) -> Seq

    if label.seqname() not in sequences:
        raise ValueError("Could not locate sequence with accession: {}".format(label.seqname()))

    seq = sequences[label.seqname()]

    if label.left() < 0 or label.right() >= len(seq):
        raise ValueError("Coordinates out of bound: [{}, {}] not in [0, {}]".format(
            label.left(), label.right(), len(seq)-1
        ))

    frag = seq[label.left():label.right()+1]

    if label.strand() == "-":
        frag = frag.reverse_complement()

    return frag


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None

    sequences = read_fasta_into_hash(args.pf_sequence)
    label = Label(
        Coordinates(args.left, args.right, args.strand),
        args.accession
    )

    try:
        frag = get_fragment(sequences, label)
        print(frag)
    except ValueError as e:
        print(e)

if __name__ == "__main__":
    main(my_env, parsed_args)
