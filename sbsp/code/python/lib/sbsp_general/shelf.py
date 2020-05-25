import logging
from typing import *
from collections import Counter
import numpy as np

from sbsp_alg.sbsp_steps import run_msa_on_sequences
from sbsp_general.general import os_join, get_value
from sbsp_options.sbsp import SBSPOptions

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


def compute_gc_from_sequences(sequences):
    # type: (Dict[str, Seq]) -> float

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

def compute_gc_from_file(pf_sequence):
    # type: (str) -> float

    from sbsp_io.sequences import read_fasta_into_hash
    sequences = read_fasta_into_hash(pf_sequence)

    return compute_gc_from_sequences(sequences)





def bin_by_gc(df, step=1):
    # type: (pd.DataFrame, int) -> List[Tuple[float, float, pd.DataFrame]]

    gc_ranges = range(20, 80, step)
    result = list()
    a = 0
    for b in gc_ranges:
        result.append(
            (a, b, df[(df["GC"] >= a) & (df["GC"]  < b)])
        )
        a = b
    return result


def get_consensus_sequence(dict_mat):
    # type: (Dict[str, List[float]]) -> str

    num_positions = len(next(iter(dict_mat.values())))
    out = ""
    for n in range(num_positions):
        best_letter = None
        best_val = None

        for letter in dict_mat.keys():
            if best_letter is None:
                best_letter = letter
                best_val = dict_mat[letter][n]
            else:
                if dict_mat[letter][n] > best_val:
                    best_letter = letter
                    best_val = dict_mat[letter][n]

        out += best_letter

    return out


def gather_consensus_sequences(env, df, col):
    # type: (Environment, pd.DataFrame, str) -> List[str]

    sequences = list()

    for idx in df.index:
        d = df.at[idx, col]     # type: Dict[str, List[float]]

        num_positions = len(next(iter(d.values())))
        out = ""
        for n in range(num_positions):
            best_letter = None
            best_val = None

            for letter in d.keys():
                if best_letter is None:
                    best_letter = letter
                    best_val = d[letter][n]
                else:
                    if d[letter][n] > best_val:
                        best_letter = letter
                        best_val = d[letter][n]

            out += best_letter
        sequences.append(out)

    return sequences


def create_extended_numpy_for_column_and_shifts(df, col, update_shifts, new_width):
    # type: (pd.DataFrame, str, List[int], int) -> np.ndarray
    df = df[~df[col].isna()]  # we only need non-NA
    example = df.at[df.index[0], col]

    n = len(df)  # number of examples
    w = new_width
    b = len(example)  # number of bases (letters)

    mat = np.zeros((n, w, b), dtype=float)

    # fill the array
    for n_pos, idx in enumerate(df.index):
        dict_arr = df.at[idx, col]

        # for each base
        for b_pos, letter in enumerate(sorted(dict_arr.keys())):
            for w_pos, value in enumerate(dict_arr[letter]):
                shifted_w_pos = w_pos + update_shifts[n_pos]
                mat[n_pos, shifted_w_pos, b_pos] = value
    return mat


def sort_sequences_by_first_non_gap_and_consensus(list_seqs):
    # type: (List[str]) -> List[str]

    def first_non_gap(l_seq):
        # type: (str) -> int
        for p in range(len(l_seq)):
            if l_seq[p] != "-":
                return p

        raise ValueError("Sequence is all gaps")

    pos_to_list_seqs = dict()
    for l in list_seqs:
        p = first_non_gap(l)
        if p not in pos_to_list_seqs.keys():
            pos_to_list_seqs[p] = list()

        pos_to_list_seqs[p].append(l)


    # reappend into single list and sort per position
    output = list()
    output_counts = list()
    for p in sorted(pos_to_list_seqs.keys()):
        # get counts per item
        counter = Counter(pos_to_list_seqs[p])
        sorted_counter = sorted([(x, counter[x]) for x in counter], key= lambda v: v[0])

        output += [x[0] for x in sorted_counter]
        output_counts += [x[1] for x in sorted_counter]

    return output, output_counts


def print_reduced_msa(msa_t, sort_by_starting_position=False, n=None):
    # type: (MSAType, bool) -> str

    list_sequences = [x.seq._data for x in msa_t.list_alignment_sequences]

    if sort_by_starting_position:
        list_sequences, counts = sort_sequences_by_first_non_gap_and_consensus(list_sequences)

    out = ""
    counter = 0
    for s, c in zip(list_sequences, counts):
        print(f"{s}\t{c}")
        out += "{}    {}\n".format(s, c)

        if n is not None and counter >= n:
            break
        counter += 1

    return out


def create_numpy_for_column_with_extended_motif(env, df, col, other=dict()):
    # type: (Environment, pd.DataFrame, str) -> np.ndarray

    example = df.at[df.index[0], col]

    # run alignment
    consensus_seqs = gather_consensus_sequences(env, df, col)
    msa_t = run_msa_on_sequences(env, consensus_seqs, SBSPOptions(env), outputorder="input-order")
    n = len(df)  # number of examples
    w = msa_t.alignment_length()
    b = len(example)  # number of bases (letters)
    other["msa_t"] = msa_t

    # get position of shift
    shifts = list()
    for s in msa_t.list_alignment_sequences:
        p = 0
        for pos in range(len(s)):
            if s[pos] != "-":
                p = pos
                break
        shifts.append(p)

    msa_t = run_msa_on_sequences(env, consensus_seqs, SBSPOptions(env), outputorder="tree-order")

    print_reduced_msa(msa_t, sort_by_starting_position=True)

    return create_extended_numpy_for_column_and_shifts(df, col, shifts, w), shifts


def get_position_distributions_by_shift(df, col, shifts):
    # type: (pd.DataFrame, str, List[int]) -> Dict[int, List[Dict[int,str]]]

    result = dict()
    for n in range(len(df.index)):
        idx = df.index[n]
        s = shifts[n]

        if s not in result:
            result[s] = list()

        result[s].append(df.at[idx, col])

    return result