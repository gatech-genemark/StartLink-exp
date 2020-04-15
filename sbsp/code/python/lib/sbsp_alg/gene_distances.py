
import os
import math
import logging
import random
import string
from typing import *
import shutil

from Bio.Phylo.PAML import codeml

from sbsp_general.general import get_value
from sbsp_io.general import generate_random_non_existing_filename, write_string_to_file, mkdir_p

logger = logging.getLogger(__name__)


def write_to_temporary_alignment_file(pf_tmp, list_sequences):
    # type: (str, List[str]) -> str

    if len(list_sequences) == 0:
        raise ValueError("No sequences to write to file")

    for s in list_sequences:
        if len(s) != len(list_sequences[0]):
            raise ValueError("Sequences should have the same length")

    text = "{} {}".format(2, len(list_sequences[0]))
    for i in range(1, len(list_sequences)+1):

        text += "\nsequence_{}  {}".format(i, list_sequences[i-1])

    text += "\n"

    write_string_to_file(text, pf_tmp)

    return pf_tmp


def _run_codeml(seq_a, seq_b, **kwargs):
    # type: (str, str, Dict[str, Any]) -> Dict[str, Any]

    pd_work = get_value(kwargs, "pd_work", ".", default_if_none=True)
    pf_ctl = get_value(kwargs, "pf_ctl", None)

    if pf_ctl is None:
        raise ValueError("Cannot compute distance without CTL file for CodeML")

    if not os.path.isfile(pf_ctl):
        raise ValueError("File doesn't exist: {}".format(pf_ctl))

    random_name = generate_random_non_existing_filename(pd_work)

    pd_codeml_run = os.path.join(pd_work, random_name)
    mkdir_p(pd_codeml_run)

    shutil.copyfile(pf_ctl, os.path.join(pd_codeml_run, "codeml.ctl"))

    pf_sequences = os.path.join(pd_codeml_run, "in.phy")
    write_to_temporary_alignment_file(pf_sequences, [seq_a, seq_b])

    write_string_to_file("(1)\n", os.path.join(pd_codeml_run, "in.tre"))

    # run code ml
    scorer = codeml.Codeml(tree=os.path.join(pd_codeml_run, "in.tre"), alignment=pf_sequences, out_file=os.path.join(pd_codeml_run, "out.txt"), working_dir=pd_codeml_run)

    try:
        results = scorer.run(ctl_file="codeml.ctl", verbose=False)
    except Exception:
        results = {}

    shutil.rmtree(pd_codeml_run)

    return results

def compute_distance_ds(seq_a, seq_b, **kwargs):
    # type: (str, str, Dict[str, Any]) -> float

    results = _run_codeml(seq_a, seq_b, **kwargs)
    try:
        return results["pairwise"]["sequence_1"]["sequence_2"]["dS"]
    except KeyError:
        return 0


def compute_distance_dn(seq_a, seq_b, **kwargs):
    # type: (str, str, Dict[str, Any]) -> float

    results = _run_codeml(seq_a, seq_b, **kwargs)
    try:
        return results["pairwise"]["sequence_1"]["sequence_2"]["dN"]
    except KeyError:
        return 0


def compute_distance_mismatch_aa(seq_a, seq_b, **kwargs):
    # type: (str, str, Dict[str, Any]) -> float

    if len(seq_a) != len(seq_b) or len(seq_a) == 0:
        raise ValueError("Sequences should have the same, non-zero length: {} vs {}".format(len(seq_a), len(seq_b)))

    matching_status = [
        seq_a[i] == seq_b[i] for i in range(len(seq_a)) if seq_a[i] != "-" and seq_b[i] != "-"
    ]

    return matching_status.count(False) / float(len(matching_status))


def compute_synonymous_fraction(aa_seq_a, aa_seq_b, nt_seq_a, nt_seq_b, **kwargs):
    # type: (str, str, str, str, Dict[str, Any]) -> float

    if len(aa_seq_a) != len(aa_seq_b) or len(aa_seq_a) == 0:
        raise ValueError("Sequences should have the same, non-zero length: {} vs {}".format(
            len(aa_seq_a), len(aa_seq_b)))

    num_aa = len(aa_seq_a)

    total_aa_no_gap = 0
    synonymous = 0
    for i_aa in range(num_aa):
        if aa_seq_a[i_aa] == "-" or aa_seq_b[i_aa] == "-":
            continue

        total_aa_no_gap += 1

        i_nt = i_aa * 3

        # check if synonymous or not
        if nt_seq_a[i_nt:i_nt+3] != nt_seq_b[i_nt:i_nt+3]:
            if aa_seq_a[i_aa] == aa_seq_b[i_aa]:
                synonymous += 1

    if total_aa_no_gap == 0:
        return 0

    return synonymous / float(total_aa_no_gap)

def compute_non_synonymous_fraction(aa_seq_a, aa_seq_b, nt_seq_a, nt_seq_b, **kwargs):
    # type: (str, str, str, str, Dict[str, Any]) -> float

    if len(aa_seq_a) != len(aa_seq_b) or len(aa_seq_a) == 0:
        raise ValueError("Sequences should have the same, non-zero length: {} vs {}".format(
            len(aa_seq_a), len(aa_seq_b)))

    num_aa = len(aa_seq_a)

    total_aa_no_gap = 0
    non_synonymous = 0
    for i_aa in range(num_aa):
        if aa_seq_a[i_aa] == "-" or aa_seq_b[i_aa] == "-":
            continue

        total_aa_no_gap += 1

        i_nt = i_aa * 3

        # check if non_synonymous or not
        if nt_seq_a[i_nt:i_nt+3] != nt_seq_b[i_nt:i_nt+3]:
            if aa_seq_a[i_aa] != aa_seq_b[i_aa]:
                non_synonymous += 1

    if total_aa_no_gap == 0:
        return 0

    return non_synonymous / float(total_aa_no_gap)


def compute_non_synonymous_poisson(aa_seq_a, aa_seq_b, nt_seq_a, nt_seq_b, **kwargs):
    # type: (str, str, str, str, Dict[str, Any]) -> float

    fraction = compute_non_synonymous_fraction(aa_seq_a, aa_seq_b, nt_seq_a, nt_seq_b, **kwargs)

    if fraction == 0:
        return 0

    return -math.log(1 - fraction)

def compute_synonymous_poisson(aa_seq_a, aa_seq_b, nt_seq_a, nt_seq_b, **kwargs):
    # type: (str, str, str, str, Dict[str, Any]) -> float

    fraction = compute_synonymous_fraction(aa_seq_a, aa_seq_b, nt_seq_a, nt_seq_b, **kwargs)

    if fraction == 0:
        return 0

    return -math.log(1 - fraction)







