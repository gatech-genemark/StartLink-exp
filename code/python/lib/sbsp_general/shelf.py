import logging
import math
import os
from typing import *
from collections import Counter
import numpy as np
import sbsp_viz
import seaborn
from matplotlib import pyplot as plt

from sbsp_alg.shelf import run_msa_on_sequences
from sbsp_container.gms2_mod import GMS2Mod
from sbsp_general.GMS2Noncoding import GMS2Noncoding
from sbsp_general.MotifModel import MotifModel
from sbsp_general.general import os_join, get_value, run_shell_cmd
from sbsp_io.general import remove_p, convert_multi_fasta_to_single
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


def relative_entropy(motif, background, component=None):
    # type: (MotifModel, GMS2Noncoding, str) -> float

    df_motif = motif.pwm_to_df()
    arr_bgd = background.pwm_to_array(0)

    result = 0.0

    if component in {"motif", "both", None}:
        for idx in df_motif.index:
            for i, l in enumerate(sorted(df_motif.columns)):
                result += df_motif.at[idx, l] * math.log2(df_motif.at[idx, l] / arr_bgd[i])

    if component in {"spacer", "both", None} and motif._spacer is not None:
        sp = motif._spacer
        sp_length = len(sp)
        for i in range(sp_length):
            result += sp[i] * math.log2(sp[i] / (1.0 / sp_length))

    return result


def run_gms2_prediction_with_model(pf_sequence, pf_new_mod, pf_new_pred):
    # type: (str, str, str) -> None

    from sbsp_general import ENV
    bin_external = ENV["pd-bin-external"]
    prog = f"{bin_external}/gms2/gmhmmp2"
    mgm_mod = f"{bin_external}/gms2/mgm_11.mod"
    cmd = f"{prog} -m {pf_new_mod} -M {mgm_mod} -s {pf_sequence} -o {pf_new_pred} --format gff"
    run_shell_cmd(cmd)


def train_gms2_model(env, pf_new_seq, pf_new_labels, **kwargs):
    group = get_value(kwargs, "group", "A", default_if_none=True)
    pf_mod = os_join(env["pd-work"], "a.mod")
    cmd = f"cd {env['pd-work']}; "
    cmd += f"{env['pd-bin-external']}/gms2/biogem gms2-training -s {pf_new_seq} -l {pf_new_labels} -m {pf_mod} --order-coding 5 --order-noncoding 2 --only-train-on-native 1 --genetic-code 11 --order-start-context 2 --fgio-dist-thr 25 --genome-group {group} --ga-upstr-len-rbs 20 --align right --ga-width-rbs 6"
    run_shell_cmd(
        cmd
    )
    mod = GMS2Mod.init_from_file(pf_mod)
    # remove_p(pf_mod)

    return mod


def train_and_create_models(env, pf_labels, pf_sequences, **kwargs):
    # type: (Environment, str, str) -> GMS2Mod
    pf_new_seq, pf_new_labels = convert_multi_fasta_to_single(env, pf_sequences, pf_labels)

    mod = train_gms2_model(env, pf_new_seq, pf_new_labels, **kwargs)
    remove_p(pf_new_labels)
    remove_p(pf_new_seq)

    return mod


def append_to_file(a_str, pf_output):
    # type: (str, str) -> None

    with open(pf_output, "a") as f:
        f.write(a_str)


def add_toolp_rbs_to_gms2_model(env, pf_sequence, pf_toolp, pf_gms2_mod, pf_new_mod, **kwargs):
    # type: (Environment, str, str, str, str) -> None

    group = get_value(kwargs, "group", None)

    # run toolp and create model file
    mod = train_and_create_models(
        env,
        pf_labels=pf_toolp,
        pf_sequences=pf_sequence,
        group=group
    )
    rbs_toolp = mod.items["RBS_MAT"]      # type: Dict[str, List[float]]
    spacer = mod.items["RBS_POS_DISTR"]

    cmd = ""

    # remove RBS_MAT and RBS_POS_DISTR from new model
    # cmd += " awk '{if ($1 == \"$RBS_MAT\") NR += 4 ; else print }' " + "{} > {}".format(pf_gms2_mod, pf_new_mod)
    cmd += "awk 'BEGIN{sut=0} {if (sut == 1) {l=substr($1,1,1);  if (l != \"$\") next ; else {sut=0; print}} "
    cmd += "else if ($1 == \"$RBS_MAT\" || $1 == \"$RBS_POS_DISTR\") sut = 1; else print }' "
    cmd += "{} > {}".format(pf_gms2_mod, pf_new_mod)
    run_shell_cmd(cmd)

    # write toolp RBS_MAT to new model file
    rbs_as_str = "\n\n$RBS_MAT\n"
    for i in sorted(rbs_toolp.keys()):
        rbs_as_str += str(i) + " " + " ".join([str(x) for x in rbs_toolp[i]]) + "\n"
    rbs_as_str += "\n\n"

    rbs_as_str += "$RBS_POS_DISTR\n"
    for i in sorted(spacer.keys()):
        rbs_as_str += str(i) + " " + str(spacer[i]) + "\n"
    rbs_as_str += "\n\n"


    append_to_file(
        rbs_as_str, pf_new_mod
    )

    return


def plot_candidate_codons(env, df, codons, cmap=None):
    # type: (Environment, pd.DataFrame, List[str]) -> None

    fig, ax = plt.subplots()
    from sbsp_viz.colormap import ColorMap as CM

    for c in sorted(codons):
        seaborn.regplot(df["GC"].astype(float).values, df[c].astype(float).values, label=c,
                        lowess=True, scatter_kws={"s": 5, "alpha": 0.1},
                        color=cmap[c]
                        )

    ax.set_ylim([-0.05,1.05])
    ax.set_ylabel("Probability")
    ax.set_xlabel("GC")
    leg = ax.legend()
    for lh in leg.legendHandles:
        lh.set_alpha(1)

    plt.show()

    # bacteria vs archaea
    fig, axes = plt.subplots(1, 2, sharex="all", sharey="all")

    for t, ax in zip(["Bacteria", "Archaea"], axes.ravel()):
        df_tmp = df[df["Type"] == t]
        for c in sorted(codons):
            seaborn.regplot(df_tmp["GC"].astype(float).values, df_tmp[c].astype(float).values, label=c,
                            lowess=True, scatter_kws={"s": 5, "alpha": 0.1}, ax=ax, color=cmap[c])

        ax.set_ylim([-0.05, 1.05])
        ax.set_ylabel("Probability")
        ax.set_xlabel("GC")
        ax.set_title(t)
        leg = ax.legend()
        for lh in leg.legendHandles:
            lh.set_alpha(1)

    plt.show()

    # group
    fig, axes = plt.subplots(2, 2, sharex="all", sharey="all")

    for t, ax in zip(list("ABCD"), axes.ravel()):
        df_tmp = df[df["GENOME_TYPE"] == t]
        for c in sorted(codons):
            seaborn.regplot(df_tmp["GC"].astype(float).values, df_tmp[c].astype(float).values, label=c,
                            lowess=True, scatter_kws={"s": 5, "alpha": 0.1}, ax=ax, color=cmap[c])

        ax.set_ylim([-0.05, 1.05])
        ax.set_ylabel("Probability")
        ax.set_xlabel("GC")
        ax.set_title(t)
        leg = ax.legend()
        for lh in leg.legendHandles:
            lh.set_alpha(1)

    plt.show()


def plot_candidate_starts(env, df):
    # type: (Environment, pd.DataFrame) -> None
    from sbsp_viz.colormap import ColorMap as CM
    plot_candidate_codons(env, df, ["ATG", "GTG", "TTG"],
                          CM.get_map("starts"))


def plot_candidate_stops(env, df):
    # type: (Environment, pd.DataFrame) -> None
    from sbsp_viz.colormap import ColorMap as CM
    plot_candidate_codons(env, df, ["TAA", "TAG", "TGA"],
                          CM.get_map("stops"))


def fix_names(r):
    # type: (pd.Series) -> str
    return "{}. {}".format(
        r["Genome"][0], r["Genome"].split("_")[1]
    )


def append_data_frame_to_csv(df, pf_output):
    # type: (pd.DataFrame, str) -> None
    if df is not None and len(df) > 0:
        try:
            if not os.path.isfile(pf_output):
                df.to_csv(pf_output, index=False)
            else:
                df.to_csv(pf_output, mode="a", index=False, header=False)
        except FileNotFoundError:
            raise ValueError(f"Could not write to file {pf_output}")
