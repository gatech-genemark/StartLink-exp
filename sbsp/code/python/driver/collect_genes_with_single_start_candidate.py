# Karl Gemayel
# Georgia Institute of Technology
#
# Created:

import os
import logging
import argparse
import pandas as pd
from typing import *
import random

# noinspection PyUnresolvedReferences
from Bio.Seq import Seq

import pathmagic                        # add path to custom library

# Custom library imports
import sbsp_general
from sbsp_alg.sequence_manipulationn import add_candidate_starts
from sbsp_container.genome_list import GenomeInfoList, GenomeInfo
from sbsp_general import Environment
from sbsp_general.general import get_value
from sbsp_general.labels import Label, Coordinates, Labels
from sbsp_io.general import mkdir_p, remove_p
from sbsp_io.labels import read_labels_from_file, write_labels_to_file
from sbsp_io.sequences import read_fasta_into_hash, write_fasta_hash_to_file

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #

parser = argparse.ArgumentParser("Description of driver.")

parser.add_argument('--pf-genome-list', required=True, help="List of genome names")
parser.add_argument('--upstream-length-nt', required=True, type=int, help="Length of upstream region. "
                                                                "Stop codons will be removed")

parser.add_argument('--add-candidates-in-upstream-length', type=int, default=99)
parser.add_argument('--add-candidates-in-downstream-length', type=int, default=99)
parser.add_argument('--num-starts-upstream', type=int, default=0)
parser.add_argument('--num-starts-downstream', type=int, default=0)

parser.add_argument('--pf-output', required=True, help="CSV output file")
parser.add_argument('--pd-output', required=True, help="Directory in which sequence.fasta and ncbi.gff will be created")

parser.add_argument('--pd-work', required=False, default=None, help="Path to working directory")
parser.add_argument('--pd-data', required=False, default=None, help="Path to data directory")
parser.add_argument('--pd-results', required=False, default=None, help="Path to results directory")
parser.add_argument("-l", "--log", dest="loglevel", choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                    help="Set the logging level", default='WARNING')
parser.add_argument("-L", "pf-log", default=None, required=False, help="The log file")


parsed_args = parser.parse_args()

# ------------------------------ #
#           Main Code            #
# ------------------------------ #

# Load environment variables
my_env = Environment(pd_data=parsed_args.pd_data,
                     pd_work=parsed_args.pd_work,
                     pd_results=parsed_args.pd_results)

# Setup logger
logging.basicConfig(level=parsed_args.loglevel, filename=parsed_args.pf_log)
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

def clean_up_stops(frag):
    # type:(str) -> str

    frag_list = list(frag)

    for i in range(0, len(frag)-10, 3):

        codon = frag[i:i+3]
        if is_valid_stop(codon, "+"):
            frag_list[i] = "C"

    return "".join(frag_list)

def add_starts(frag, pos_true):
    # type: (str, int) -> str

    frag_list = list(frag)
    pos_true_aa = int(pos_true / 3)

    for n in range(parsed_args.num_of_starts_upstream):

        try:
            new_pos = random.randint(0, pos_true_aa) * 3
        except ValueError:
            continue

        frag_list[new_pos] = "A"
        frag_list[new_pos+1] = "T"
        frag_list[new_pos+2] = "G"

    for n in range(parsed_args.num_of_starts_downstream):
        try:
            new_pos = random.randint(pos_true_aa+1, int(len(frag)/3 - 2)) * 3
        except ValueError:
            continue

        frag_list[new_pos] = "A"
        frag_list[new_pos + 1] = "T"
        frag_list[new_pos + 2] = "G"

    return "".join(frag_list)


def remove_starts_in_region(frag, begin, end):
    # type: (str, int, int) -> str
    frag_list = list(frag)

    for n in range(begin, end, 3):
        if is_valid_start(frag[n:n+3], "+"):
            frag_list[n] = "C"

    return "".join(frag_list)


def get_entry_for_label_if_single_candidate(sequences, label, tag, **kwargs):
    # type: (Dict[str, Seq], Label, int, Dict[str, Any]) -> Union[Tuple[str, Label], None]
    """Constructs a sequence/label pair for a gene if it has a single candidate start. Otherwise, returns None

    :param sequences:
    :param label:
    :param tag:
    :param kwargs:
    :return:
    """

    upstream_length_nt = get_value(kwargs, "upstream_length_nt", 180, default_if_none=True)
    acs_kws = get_value(kwargs, "acs_kws", init=dict)
    min_gene_length = 48

    if upstream_length_nt % 3 != 0:
        raise ValueError("Upstream length must be a multiple of 3")

    num_candidates = 1
    seq = sequences[label.seqname()]

    left = label.coordinates().left
    right = label.coordinates().right
    strand = label.strand()


    if strand == "+":

        # check that label actually represents a valid gene
        if not (is_valid_start(seq[left:left+3], strand) and is_valid_stop(seq[right-2:right+1], strand)):
            return None

        # search upstream for candidate
        curr_pos = left - 3
        while curr_pos >= 0:
            codon = seq[curr_pos:curr_pos+3]
            if is_valid_start(codon, strand):
                num_candidates += 1
            if is_valid_stop(codon, strand):
                break
            curr_pos -= 3

        # search downstream for candidate c
        curr_pos = left + 3
        while curr_pos <= right - 3 - min_gene_length:
            codon = seq[curr_pos:curr_pos+3]
            if is_valid_start(codon, strand):
                num_candidates += 1
            if is_valid_stop(codon, strand):
                break
            curr_pos += 3

    if strand == "-":
        if not (is_valid_stop(seq[left:left+3], strand) and is_valid_start(seq[right-2:right+1], strand)):
            return None

        # upstream
        curr_pos = right + 3
        while curr_pos < len(seq):
            codon = seq[curr_pos-2:curr_pos + 1]
            if is_valid_start(codon, strand):
                num_candidates += 1
            if is_valid_stop(codon, strand):
                break

            curr_pos += 3

        # downstream
        curr_pos = right - 3
        while curr_pos > left + 3 + min_gene_length:
            codon = seq[curr_pos - 2:curr_pos + 1]
            if is_valid_start(codon, strand):
                num_candidates += 1
            if is_valid_stop(codon, strand):
                break
            curr_pos -= 3

    if num_candidates == 1:

        if strand == "+":
            frag_left = left - (upstream_length_nt)
            if frag_left >= 0:
                frag = seq[frag_left:right+1]
            else:
                return None
        else:
            frag_right = right + (upstream_length_nt)
            if frag_right >= len(seq):
                return None
            else:
                frag = seq[left:frag_right+1]
                frag = frag.reverse_complement()

        if len(frag) == 0:
            return None

        frag = remove_starts_in_region(frag, 0, upstream_length_nt)

        new_label = Label(
                Coordinates(upstream_length_nt, len(frag) - 1, "+"),
                seqname=label.seqname() + "_tag{}".format(tag),
                attributes=label._attributes
        )

        # remove stop and add candidate starts
        frag = clean_up_stops(frag)
        frag = add_candidate_starts(frag, new_label, **acs_kws)

        result = (
            frag.strip(),
            new_label
        )

        return result

    return None


def get_single_candidate_genes(env, genome_info, **kwargs):
    # type: (Environment, GenomeInfo, Dict[str, Any]) -> List[Tuple[str, Label]]

    gcfid = genome_info.name
    pd_gcfid = os.path.join(env["pd-data"], gcfid)

    pf_sequences = os.path.join(pd_gcfid, "sequence.fasta")
    pf_ncbi = os.path.join(pd_gcfid, "ncbi.gff")

    sequences = read_fasta_into_hash(pf_sequences)
    labels = read_labels_from_file(pf_ncbi)

    result = list()         # type: List[Tuple[str, Label]]

    counter = 0
    for lab in labels:

        # skip partial, frame-shifted, and hypothetical genes
        if lab.is_partial() or lab.is_frameshifted() or lab.is_hypothetical():
            continue

        entry = get_entry_for_label_if_single_candidate(sequences, lab, tag=counter, **kwargs)

        if entry is not None:
            result.append(entry)
            counter += 1

    return result


def append_data_frame_to_csv(df, pf_output):
    # type: (pd.DataFrame, str) -> None
    if df is not None and len(df) > 0:
        if not os.path.isfile(pf_output):
            df.to_csv(pf_output, index=False)
        else:
            df.to_csv(pf_output, mode="a", index=False, header=False)

def main(env, args):
    # type: (Environment, argparse.Namespace) -> None

    gil = GenomeInfoList.init_from_file(args.pf_genome_list)

    # sequences = dict()          # type: Dict[str, str]
    # labels = Labels()

    remove_p(args.pf_output)            # start fresh
    num_genomes = len(gil)

    for counter, genome_info in enumerate(gil):
        logger.info("{}/{}: {}".format(counter, num_genomes, genome_info.name))
        list_sequence_label_pair = get_single_candidate_genes(
            env,
            genome_info,
            upstream_length_nt=args.upstream_length_nt,
            acs_kws={
                "num_starts_upstream": args.num_starts_upstream,
                "num_starts_downstream": args.num_starts_downstream,
                "upstream_length_nt": args.add_candidates_in_upstream_length,
                "downstream_length_nt": args.add_candidates_in_downstream_length,
            }

        )

        list_entries = list()

        for sl_pair in list_sequence_label_pair:
            seq, lab = sl_pair

            # sanity check: seqname is unique
            # if lab.seqname() in sequences.keys():
            #     raise RuntimeError("All sequences should have separate FASTA tags")

            # sequences[lab.seqname()] = seq
            # labels.add(lab)

            list_entries.append({
                "accession": lab.seqname(),
                "left": lab.left() + 1,
                "right": lab.right() + 1,
                "strand": lab.strand(),
                "sequence": seq,
                "num-upstream": args.num_starts_upstream,
                "num-downstream": args.num_starts_downstream,

            })

        df = pd.DataFrame(list_entries)
        append_data_frame_to_csv(df, args.pf_output)
    # print labels and sequences to file
    # mkdir_p(args.pd_output)

    # write_fasta_hash_to_file(sequences, os.path.join(args.pd_output, "sequence.fasta"))
    # write_labels_to_file(labels, os.path.join(args.pd_output, "ncbi.gff"))




if __name__ == "__main__":
    main(my_env, parsed_args)
