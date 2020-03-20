# Karl Gemayel
# Georgia Institute of Technology
#
# Created: 2/5/20

import logging
import argparse
import os

import pandas as pd
from typing import *

# noinspection All
from Bio.Seq import Seq

import pathmagic

# noinspection PyUnresolvedReferences
import sbsp_log  # runs init in sbsp_log and configures logger

# Custom imports
from sbsp_container.genome_list import GenomeInfoList, GenomeInfo
from sbsp_general import Environment

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from sbsp_general.general import get_value
from sbsp_general.labels import Labels, Label
from sbsp_io.labels import read_labels_from_file
from sbsp_io.sequences import read_fasta_into_hash

parser = argparse.ArgumentParser("Description of driver.")

parser.add_argument('--pf-genome-list', required=True, help="File containing genome information")
parser.add_argument('--pf-output', required=True, help="Output file")
parser.add_argument('--max-upstream-length-nt', required=False, type=int, help="Max length upstream region "
                                                                               "(until LORF) in which to count")

parser.add_argument('--max-downstream-length-nt', required=False, type=int, help="Max length downstream region "
                                                                                 "in which to count")

parser.add_argument('--fn-labels', required=False, default="verified.gff")


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

valid_starts_pos = ["ATG", "GTG", "TTG"]
valid_starts_neg = ["CAT", "CAC", "CAA"]

valid_stops_pos = ["TAA", "TGA", "TAG"]
valid_stops_neg = ["TTA", "TCA", "CTA"]


def is_valid_start(codon, strand):
    if strand == "+":
        return codon in valid_starts_pos
    else:
        return codon in valid_starts_neg


def is_valid_stop(codon, strand):
    if strand == "+":
        return codon in valid_stops_pos
    else:
        return codon in valid_stops_neg


def initialize_stats():
    # type: () -> Dict[str, int]
    stats = dict()

    for direction in ["downstream", "upstream"]:
        stats[direction] = 0
        for v in ["ATG", "GTG", "TTG"]:
            stats["{}-{}".format(direction, v.lower())] = 0
            stats["{}".format(v.lower())] = 0

    stats["num-upstream-until-lorf"] = False

    return stats


def count_candidates_on_positive_strand(sequence, label, **kwargs):
    # type: (Seq, Label, Dict[str, Any]) -> Union[None, Dict[str, int]\
    max_upstream_length_nt = get_value(kwargs, "max_upstream_length_nt", None)
    max_downstream_length_nt = get_value(kwargs, "max_downstream_length_nt", None)

    stats = initialize_stats()

    pos_5prime = label.left()

    codon = str(sequence[pos_5prime:pos_5prime+3])
    if not is_valid_start(codon, "+"):
        return None

    # at position
    stats["{}".format(
        str(sequence[pos_5prime:pos_5prime+3]).lower()
    )] = 1

    num_upstream_until_lorf = 0
    stop_counting = False
    # upstream
    curr_pos = pos_5prime - 3
    while curr_pos >= 0:
        if max_upstream_length_nt is not None:
            if pos_5prime - curr_pos + 1 > max_upstream_length_nt:
                stop_counting = True

        codon = str(sequence[curr_pos:curr_pos+3])
        if is_valid_stop(codon, "+"):
            break

        if not stop_counting and is_valid_start(codon, "+"):
            stats["upstream"] += 1
            stats["upstream-{}".format(codon.lower())] += 1

            num_upstream_until_lorf += 1

        curr_pos -= 3

    stats["num-upstream-until-lorf"] = num_upstream_until_lorf

    # downstream
    curr_pos = pos_5prime + 3
    while curr_pos <= label.right() - 3:
        if max_downstream_length_nt is not None:
            if curr_pos - pos_5prime + 1 > max_downstream_length_nt:
                break

        codon = str(sequence[curr_pos:curr_pos+3])
        if is_valid_start(codon, "+"):
            stats["downstream"] += 1
            stats["downstream-{}".format(codon.lower())] += 1

        curr_pos += 3

    return stats


def count_candidates_on_negative_strand(sequence, label, **kwargs):
    # type: (Seq, Label, Dict[str, Any]) -> Union[None, Dict[str, int]]
    max_upstream_length_nt = get_value(kwargs, "max_upstream_length_nt", None)
    max_downstream_length_nt = get_value(kwargs, "max_downstream_length_nt", None)

    stats = initialize_stats()

    pos_5prime = label.right()

    def convert(element):
        conversion = {
            "CAT": "ATG",
            "CAA": "TTG",
            "CAC": "GTG"
        }

        return conversion[str(element)]

    codon = str(sequence[pos_5prime-2:pos_5prime+1])
    if not is_valid_start(codon, "-"):
        return None

    # at position
    stats["{}".format(
        convert(str(sequence[pos_5prime-2:pos_5prime+1])).lower()
    )] = 1

    num_upstream_until_lorf = 0
    stop_counting = False

    # upstream
    curr_pos = pos_5prime + 3
    while curr_pos < len(sequence):
        if max_upstream_length_nt is not None:
            if curr_pos - pos_5prime + 1 > max_upstream_length_nt:
                stop_counting = True

        codon = str(sequence[curr_pos-2:curr_pos + 1])
        if is_valid_stop(codon, "-"):
            break
        if not stop_counting and is_valid_start(codon, "0"):
            stats["upstream"] += 1
            stats["upstream-{}".format(convert(codon).lower())] += 1

            num_upstream_until_lorf += 1

        curr_pos += 3

    stats["num-upstream-until-lorf"] = num_upstream_until_lorf

    # downstream
    curr_pos = pos_5prime - 3
    while curr_pos > label.left()+2:
        if max_downstream_length_nt is not None:
            if pos_5prime - curr_pos + 1 > max_downstream_length_nt:
                break

        codon = str(sequence[curr_pos-2:curr_pos+1])
        if is_valid_start(codon, "-"):
            stats["downstream"] += 1
            stats["downstream-{}".format(convert(codon).lower())] += 1

        curr_pos -= 3

    return stats


def count_candidates_per_gene(sequences, labels, **kwargs):
    # type: (Dict[str, Seq], Labels, Dict[str, Any]) -> pd.DataFrame

    list_stats = list()

    for label in labels:
        try:
            sequence = sequences[label.seqname()]
        except KeyError:
            logger.warning("Could not get sequence for label witn sequence name: {}".format(label.seqname()))
            continue

        if label.strand() == "+":
            stats = count_candidates_on_positive_strand(sequence, label, **kwargs)
        else:
            stats = count_candidates_on_negative_strand(sequence, label, **kwargs)

        if stats is not None:
            list_stats.append(stats)

    df = pd.DataFrame(list_stats)
    df = df.reindex(sorted(df.columns), axis=1)
    return df


def read_sequences_for_genome(env, genome_info):
    # type: (Environment, GenomeInfo) -> Dict[str, Seq]
    gcfid = genome_info.name
    pd_gcfid = os.path.join(env["pd-data"], gcfid)

    pf_sequences = os.path.join(pd_gcfid, "sequence.fasta")

    return read_fasta_into_hash(pf_sequences)


def read_labels_for_genome(env, genome_info, **kwargs):
    # type: (Environment, GenomeInfo, Dict[str, Any]) -> Labels
    fn_labels = get_value(kwargs, "fn_labels", "ncbi.gff")
    gcfid = genome_info.name
    pd_gcfid = os.path.join(env["pd-data"], gcfid)
    pf_ncbi = os.path.join(pd_gcfid, fn_labels)
    return read_labels_from_file(pf_ncbi)


def count_candidates_per_gene_for_genomes(env, gil, **kwargs):
    # type: (Environment, GenomeInfoList, Dict[str, Any]) -> pd.DataFrame

    list_df = list()
    for gi in gil:
        sequences = read_sequences_for_genome(env, gi)
        labels = read_labels_for_genome(env, gi, **kwargs)

        df_gi = count_candidates_per_gene(sequences, labels, **kwargs)
        df_gi["gcfid"] = gi.name
        df_gi["ancestor"] = gi.attributes["ancestor"]

        list_df.append(df_gi)

    return pd.concat(list_df, ignore_index=True)


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None

    gil = GenomeInfoList.init_from_file(args.pf_genome_list)

    df = count_candidates_per_gene_for_genomes(
        env, gil,
        max_upstream_length_nt=args.max_upstream_length_nt,
        max_downstream_length_nt=args.max_downstream_length_nt,
        fn_labels=args.fn_labels

    )

    df.to_csv(args.pf_output, index=False)

if __name__ == "__main__":
    main(my_env, parsed_args)
