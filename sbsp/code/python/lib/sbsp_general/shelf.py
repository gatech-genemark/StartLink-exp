import logging
from typing import *

from sbsp_general.general import os_join, get_value

logger = logging.getLogger(__name__)


def create_q_key_3p(accession, left, right, strand):
    # type: (str, int, int, str) -> str
    if strand == "+":
        return "{};{};{};{}".format(accession, "", right, strand)
    else:
        return "{};{};{};{}".format(accession, left, "", strand)


def create_q_key_5p_3p(accession, left, right, strand):
    # type: (str, int, int, str) -> str
    if strand == "+":
        return "{};{};{};{}".format(accession, left, right, strand)
    else:
        return "{};{};{};{}".format(accession, left, right, strand)


def add_q_key_3p_to_df(df, column):
    # type: (pd.DataFrame, str) -> None

    df[column] = df.apply(
        lambda row: create_q_key_3p(
            row["q-accession"], int(row["q-left"]), int(row["q-right"]), row["q-strand"]
        ),
        axis=1
    )


def map_key_3p_to_label(labels):
    # type: (Labels) -> Dict[str, Label]

    return {
        create_q_key_3p(l.seqname(), l.left(), l.right(), l.strand()): l for l in labels
    }


def map_key_3p_to_df_group(df_sbsp_details):
    # type; (pd.DataFrame) -> Dict[str, pd.DataFrame]

    return {
        key: df_group for key, df_group in df_sbsp_details.groupby("q-key-3p", as_index=False)
    }


def labels_match_5p_3p(label_a, label_b):
    # type: (Label, Label) -> bool

    key_a = create_q_key_5p_3p(label_a.seqname(), label_a.left(), label_a.right(), label_a.strand())
    key_b = create_q_key_5p_3p(label_b.seqname(), label_b.left(), label_b.right(), label_b.strand())

    return key_a == key_b


def next_name(pd_work, **kwargs):
    # type: (str, Dict[str, Any]) -> str

    ext = get_value(kwargs, "ext", "pdf")
    if "counter" not in next_name.__dict__: next_name.counter = -1
    next_name.counter += 1
    return os_join(pd_work, "{}.{}".format(next_name.counter, ext))


def compute_gc_from_file(pf_sequence):
    # type: (str) -> float

    from sbsp_io.sequences import read_fasta_into_hash
    sequences = read_fasta_into_hash(pf_sequence)


    counts = {"A": 0, "C": 0, "G": 0, "T": 0}

    for seq in sequences.values():
        for s in seq:
            if s.upper() in {"A", "C", "G", "T"}:
                counts[s] += 1

    total = sum(counts.values())
    count_gc = counts["G"] + counts["C"]

    if total == 0:
        return 0.0

    return count_gc / float(total)
