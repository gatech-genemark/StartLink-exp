# Karl Gemayel
# Georgia Institute of Technology
#
# Created: 08/30/2020
import copy
import logging
import argparse
import os

import pandas as pd
from typing import *

# noinspection All
import seaborn

import pathmagic

# noinspection PyUnresolvedReferences
import sbsp_log  # runs init in sbsp_log and configures logger

# Custom imports
from sbsp_container.genome_list import GenomeInfoList, GenomeInfo
from sbsp_general import Environment
from sbsp_general.general import get_value, os_join, run_shell_cmd
from sbsp_general.labels import Labels
from sbsp_general.labels_comparison_detailed import LabelsComparisonDetailed
from sbsp_general.shelf import compute_gc_from_file, next_name
from sbsp_io.general import mkdir_p
from sbsp_io.labels import read_labels_from_file
from sbsp_io.objects import save_obj, load_obj
from sbsp_options.parallelization import ParallelizationOptions
from sbsp_parallelization.pbs import PBS
from sbsp_pbs_data.mergers import merge_identity
from sbsp_pbs_data.splitters import split_genome_info_list
import sbsp_argparse.parallelization
from Bio import SeqIO
# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from sbsp_viz.colormap import ColorMap
from sbsp_viz.general import FigureOptions
from sbsp_viz.sns import lmplot

parser = argparse.ArgumentParser("Run external prediction tools on genome list.")

parser.add_argument('--pf-sequence', required=True)
parser.add_argument('--pf-labels', required=True)

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


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None

    sequences = SeqIO.to_dict(SeqIO.parse(args.pf_sequence, "fasta"))
    labels = read_labels_from_file(args.pf_labels, shift=-1)

    for lab in labels:
        seq = sequences[lab.seqname()]

        left = int(lab.left())
        right = int(lab.right())

        if lab.strand() == "+":
            codon = seq[lab.left():lab.left()+3]
        else:
            codon = seq[lab.right()-2:lab.right() + 1] # type: SeqIO.SeqRecord
            codon = codon.reverse_complement()


        print(codon.seq._data)


if __name__ == "__main__":
    main(my_env, parsed_args)
