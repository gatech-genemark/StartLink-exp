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
from sbsp_general import Environment
from sbsp_general.labels import Label, Coordinates, Labels
from sbsp_io.labels import read_labels_from_file
from sbsp_io.sequences import read_fasta_into_hash

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #

parser = argparse.ArgumentParser("Description of driver.")

parser.add_argument('--pf-genome-list', required=True, help="List of genome names")
parser.add_argument('--num-of-starts-upstream', type=int, required=True, help="Number of starts to add upstream of true")
parser.add_argument('--num-of-starts-downstream', type=int, required=True, help="Number of starts to add downstream of true")

parser.add_argument('--pf-sequences-output', required=True)
parser.add_argument('--pf-labels-output', required=True)

parser.add_argument('--mean-inter-start-distance-aa', type=int, required=False, default=5, help="Distance between starts")

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


def get_entry_for_label(sequences, label, tag):
    # type: (Dict[str, Seq], Label, int) -> Union[None, Dict[str, Any]]

    result = dict()

    upstream_aa = 60

    num_candidates = 1
    seq = sequences[label.seqname()]

    left = label.coordinates().left
    right = label.coordinates().right
    strand = label.strand()
    min_gene_length = 48

    if strand == "+":

        if not (is_valid_start(seq[left:left+3], strand) and is_valid_stop(seq[right-2:right+1], strand)):
            return None

        # upstream
        curr_pos = left - 3
        while curr_pos >= 0:
            codon = seq[curr_pos:curr_pos+3]
            if is_valid_start(codon, strand):
                num_candidates += 1
            if is_valid_stop(codon, strand):
                break
            curr_pos -= 3

        # downstream
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
            frag_left = left - (upstream_aa*3)
            if frag_left >= 0:
                frag = seq[frag_left:right+1]
            else:
                return None
        else:
            frag_right = right + (upstream_aa*3)
            if frag_right >= len(seq):
                return None
            else:
                frag = seq[left:frag_right+1]
                frag = frag.reverse_complement()

        if len(frag) == 0:
            return None

        frag = clean_up_stops(frag)
        frag = add_starts(frag, upstream_aa*3)

        result = {
            "sequence": frag.strip(),
            "label": Label(
                Coordinates(upstream_aa*3, len(frag)-1, "+"),
                seqname=label.seqname() + "_tag{}".format(tag),
                attributes=label._attributes
            )
        }

        return result

    return None

def get_sequences_for_single_genome(env, gcfid):
    # type: (Environment, str) -> List[Dict[str, Any]]

    pd_gcfid = os.path.join(env["pd-data"], gcfid)

    pf_sequences = os.path.join(pd_gcfid, "sequence.fasta")
    pf_ncbi = os.path.join(pd_gcfid, "ncbi.gff")

    sequences = read_fasta_into_hash(pf_sequences)
    labels = read_labels_from_file(pf_ncbi)

    result = list()

    counter = 0
    for lab in labels:

        if lab.is_partial():
            continue

        if lab.is_frameshifted():
            continue

        # not hypothetical
        if lab.get_attribute_value("product") is not None and "hypothetical" in lab.get_attribute_value("product"):
            continue

        entry = get_entry_for_label(sequences, lab, tag=counter)

        if entry is not None:
            result.append(entry)
            counter += 1

    return result



def get_sequences_per_genome(env, df_genomes):
    # type: (Environment, pd.DataFrame) -> List[Dict[str, Any]]

    results = list()
    for index, row in df_genomes.iterrows():

        tmp_results = get_sequences_for_single_genome(env, row.iloc[0])

        results += tmp_results

    return results



def main(env, args):
    # type: (Environment, argparse.Namespace) -> None

    from sbsp_io.general import read_rows_to_list
    my_list = read_rows_to_list(args.pf_genome_list)
    my_list = [l.strip().split()[0] for l in my_list]


    df_genomes = pd.DataFrame({"gcfid": my_list})

    info = get_sequences_per_genome(env, df_genomes)

    accession_to_seq = {
        d["label"].seqname(): d["sequence"] for d in info
    }

    labels = Labels(
        [d["label"] for d in info]
    )

    from sbsp_io.labels import write_labels_to_file
    write_labels_to_file(labels, args.pf_labels_output)

    from sbsp_io.sequences import write_fasta_hash_to_file
    write_fasta_hash_to_file(accession_to_seq, args.pf_sequences_output)


if __name__ == "__main__":
    main(my_env, parsed_args)
