# Karl Gemayel
# Georgia Institute of Technology
#
# Created: 3/17/20
import math

import umap
import umap.plot

import logging
import argparse
import numpy as np
import pandas as pd
from typing import *

# noinspection All
from Bio.Align import AlignInfo, MultipleSeqAlignment

# noinspection PyUnresolvedReferences
import sbsp_log  # runs init in sbsp_log and configures logger

# Custom imports
from sbsp_alg.sbsp_steps import run_msa_on_sequences
from sbsp_general import Environment
from sbsp_general.msa_2 import MSAType
from sbsp_general.shelf import next_name, bin_by_gc, get_consensus_sequence, gather_consensus_sequences, \
    print_reduced_msa, create_numpy_for_column_with_extended_motif, get_position_distributions_by_shift
from sbsp_io.objects import load_obj
import seaborn
import matplotlib.pyplot as plt
import logomaker as lm
from matplotlib.font_manager import FontProperties



# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from sbsp_options.sbsp import SBSPOptions
from sbsp_viz.general import save_figure, FigureOptions
from sbsp_viz.shelf import loess_with_stde, create_mappable_for_colorbar

parser = argparse.ArgumentParser("Description of driver.")

parser.add_argument('--pf-input-arc', required=True, help="Input file")
parser.add_argument('--pf-input-bac', required=True, help="Input file")

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



def create_numpy_for_column(df, col):
    # type: (pd.DataFrame, str) -> np.ndarray

    df = df[~df[col].isna()]            # we only need non-NA
    example = df.at[df.index[0], col]

    n = len(df)                         # number of examples
    w = len(next(iter(example.values())))        # width (numbere of positions)
    b = len(example)                    # number of bases (letters)

    mat = np.empty((n,w,b), dtype=float)

    # fill the array
    for n_pos, idx in enumerate(df.index):
        dict_arr = df.at[idx, col]


        # for each base
        for b_pos, letter in enumerate(sorted(dict_arr.keys())):
            for w_pos, value in enumerate(dict_arr[letter]):
                if w_pos == 6:
                    print("HI")
                mat[n_pos, w_pos, b_pos] = value

    return mat


def visualize_matrix_column(env, df, col):
    # type: (Environment, pd.DataFrame, str) -> None

    # first, remove all NA for column
    df = df[~df[col].isna()]           # we only need non-NA


    fp = FontProperties()
    fp.set_family("monospace")

    # create N x 6 x 4 matrix for RBS
    mat = create_numpy_for_column(df, col)
    mat = mat.reshape((mat.shape[0], mat.shape[1] * mat.shape[2]))

    # get interesting features to view data by
    gc = df["GC"]
    group = df["GENOME_TYPE"]

    for r in range(1):

        reducer = umap.UMAP(random_state=r)
        reducer = reducer.fit(mat)
        embedding = reducer.embedding_
        print(embedding.shape)

        # fig, ax = plt.subplots()
        #
        # plt.scatter(embedding[:, 0], embedding[:, 1], c=gc, marker="+")
        # plt.colorbar()
        # plt.show()
        # themes = ["fire", "viridis", "inferno", "blue", "red", "green", "darkblue", "darkred", "darkgreen"]
        # fig, axes = plt.subplots(3, 3)
        # for ax, theme in zip(axes.ravel(), themes):
        #     fig, ax = plt.subplots()
        #     umap.plot.points(reducer, values=gc, theme=theme, )
        #     plt.show()
        ax = umap.plot.points(reducer, values=gc, cmap="viridis")
        mappable = create_mappable_for_colorbar(gc, "viridis")
        plt.colorbar(mappable)
        plt.title(col)
        plt.tight_layout()
        save_figure(FigureOptions(
            save_fig=next_name(env["pd-work"])
        ))
        plt.show()

        umap.plot.points(reducer, labels=group.values, color_key_cmap="Paired")
        plt.title(col)
        plt.tight_layout()
        save_figure(FigureOptions(
            save_fig=next_name(env["pd-work"])
        ))
        plt.show()

        # umap.plot.points(reducer, labels=group.values, color_key_cmap="Dark2")
        # plt.title(col)
        # save_figure(FigureOptions(
        #     save_fig=next_name(env["pd-work"])
        # ))
        # plt.show()


        umap.plot.points(reducer, labels=df["Type"])
        plt.title(col)
        plt.tight_layout()
        save_figure(FigureOptions(
            save_fig=next_name(env["pd-work"])
        ))
        plt.show()

    # fig, ax = plt.subplots()
    # plt.scatter(embedding[:, 0], embedding[:, 1], c=group.values, cmap=cm.brg)
    # plt.show()


def build_consensus_from_consensus(env, df, col):
    # type: (Environment, pd.DataFrame, str) -> None
    df = df[~df[col].isna()] .copy()           # we only need non-NA

    consensus_seqs = gather_consensus_sequences(env, df, col)

    msa_t = run_msa_on_sequences(env, consensus_seqs, SBSPOptions(env, gapopen=10000), outputorder="tree-order")

    # print(msa_t.to_string())
    #
    # print(consensus_seqs)

    summary_align = AlignInfo.SummaryInfo(MultipleSeqAlignment(msa_t.list_alignment_sequences))
    con = summary_align.dumb_consensus()

    # print(con)
    # print(summary_align)
    seqs = [x.seq._data for x in  msa_t.list_alignment_sequences]
    counts_mat = lm.alignment_to_matrix(sequences=seqs, to_type='counts', characters_to_ignore='.-X')

    # Counts matrix -> Information matrix
    info_mat = lm.transform_matrix(counts_mat,
                                   from_type='counts',
                                   to_type='information')

    lm.Logo(info_mat)
    plt.show()

    from collections import Counter

    print("New set")
    counter = Counter(consensus_seqs)
    sorted_counter = counter.most_common()
    print("\n".join([str(x) for x in sorted_counter]))


def get_letter_probabilities_and_gc_for_position(df, col, p):
    # type: (pd.DataFrame, str, int) -> Dict[str, np.ndarray]

    example = df.at[df.index[0], col]
    n = len(df)

    result = {
        x: np.empty((n, 2), dtype=float) for x in example.keys()
    }

    for n_pos, idx in enumerate(df.index):
        d = df.at[idx, col]     # type: Dict[str, List[float]]
        for l in d.keys():
            result[l][n_pos, 0] = d[l][p]
            result[l][n_pos, 1] = df.at[idx, "GC"]

    return result


def plot_letter_over_gc(env, df, col):
    # type: (Environment, pd.DataFrame, str) -> None

    for p in range(6):

        # get data per letter for position
        dict_letter_prob_at_gc = get_letter_probabilities_and_gc_for_position(df, col, p)

        letters = dict_letter_prob_at_gc.keys()

        fig, axes = plt.subplots(2, math.ceil(len(letters)/2), sharex="all", sharey="all")
        for l, ax in zip(letters, axes.ravel()[:len(letters)]):
            ax.scatter(dict_letter_prob_at_gc[l][:, 1], dict_letter_prob_at_gc[l][:, 0],
                       marker="+")
            ax.set_title(f"{p, l}")
        plt.show()


def find_most_populated_block(msa_t, width):
    # type: (MSAType, int) -> int
    population_per_position = [0] * width

    begin = None     # indicates position of first column in block  (inclusive)
    end = None        # indicates position of last column in block  (inclusive)

    if width > msa_t.alignment_length():
        return 0

    def compute_population_at_pos(local_pos):
        # type: (int) -> int
        return sum(1 for i in range(msa_t.number_of_sequences()) if msa_t[i][local_pos] != "-")

    # get first block
    for p in range(width):
        population_per_position[p] = compute_population_at_pos(p)

    best_pop = sum(population_per_position)
    best_pos = 0

    curr_begin = 0
    for p in range(width, msa_t.alignment_length()):
        pop = compute_population_at_pos(p)
        population_per_position[curr_begin] = pop

        new_pop = sum(population_per_position)
        if new_pop > best_pop:
            best_pop = new_pop
            best_pos = p - width

        curr_begin += 1
        if curr_begin == len(population_per_position):
            curr_begin = 0

    return best_pos


def create_numpy_for_column_and_shifts(df, col, update_shifts):
    # type: (pd.DataFrame, str, List[int]) -> np.ndarray
    df = df[~df[col].isna()]  # we only need non-NA
    example = df.at[df.index[0], col]

    n = len(df)  # number of examples
    w = len(next(iter(example.values())))  # width (numbere of positions)
    b = len(example)  # number of bases (letters)

    mat = np.empty((n, w, b), dtype=float)

    # fill the array
    for n_pos, idx in enumerate(df.index):
        dict_arr = df.at[idx, col]

        # for each base
        for b_pos, letter in enumerate(sorted(dict_arr.keys())):
            for w_pos, value in enumerate(dict_arr[letter]):
                shifted_w_pos = w_pos + update_shifts[n_pos]
                if shifted_w_pos < 0 or shifted_w_pos >= w:
                    continue

                mat[n_pos, shifted_w_pos, b_pos] = value
    return mat


def create_numpy_for_column_and_realign(env, df, col):
    # type: (Environment, pd.DataFrame, str) -> np.ndarray

    example = df.at[df.index[0], col]

    n = len(df)  # number of examples
    w = len(next(iter(example.values())))  # width (numbere of positions)
    b = len(example)  # number of bases (letters)

    # run alignment
    consensus_seqs = gather_consensus_sequences(env, df, col)
    msa_t = run_msa_on_sequences(env, consensus_seqs, SBSPOptions(env), outputorder="input-order")

    # get position of shift
    shifts = list()
    for s in msa_t.list_alignment_sequences:
        p = 0
        for pos in range(len(s)):
            if s[pos] != "-":
                p = pos
                break
        shifts.append(p)

    # find position of block in alignment < w where most sequences exist
    pos_of_block = find_most_populated_block(msa_t, w)

    update_shifts = [x-pos_of_block for x in shifts]
    return create_numpy_for_column_and_shifts(df, col, update_shifts), update_shifts




def plot_letter_over_gc_with_alignment_per_bin(env, df, col):
    # type: (Environment, pd.DataFrame, str) -> None
    binned_dfs = bin_by_gc(df)
    binned_arrays = list()
    # for each bin
    for info in binned_dfs:
        lower, upper, df_gc = info

        if len(df_gc) <= 1:
            continue

        array, update_shifts = create_numpy_for_column_and_realign(env, df_gc, col)


        binned_arrays.append({
            "GC": df_gc["GC"],
            "motifs": array,
            "shifts": update_shifts
        })

    example = df.at[df.index[0], col]       # type: Dict[str, List[float]]
    w = len(next(iter(example.values())))  # width (numbere of positions)
    b = len(example)  # number of bases (letters)

    for w_pos in range(w):
        letters = example.keys()
        letter_to_idx = {x: x_pos for x_pos, x in enumerate(sorted(letters))}

        fig, axes = plt.subplots(2, math.ceil(len(letters) / 2), sharex="all", sharey="all")
        for l, ax in zip(letters, axes.ravel()[:len(letters)]):
            # go through GC bins
            all_gc = list()
            all_probs = list()

            for ba in binned_arrays:
                arr = ba["motifs"]
                gc = ba["GC"].values
                shifts = ba["shifts"]

                for index in range(len(shifts)):
                    if shifts[index] > 0:
                        if w_pos < shifts[index]:
                            continue
                    elif shifts[index] < 0:
                        if w_pos >= w + shifts[index]:
                            continue


                    # shifted_pos = w_pos - shifts[index]
                    # if shifted_pos < 0 or shifted_pos >= w:
                    #     continue

                    all_gc.append(gc[index])
                    if arr[index, w_pos, letter_to_idx[l]] < 0 or arr[index, w_pos, letter_to_idx[l]] > 1:
                        raise ValueError("Something's up")
                        break
                    all_probs.append(arr[index, w_pos, letter_to_idx[l]])

            ax.scatter(all_gc, all_probs, marker="+")
            ax.set_title(f"{l}")
            # plt.show()
        plt.show()


def plot_letter_over_gc_with_single_alignment(env, df, col):
    # type: (Environment, pd.DataFrame, str) -> None
    binned_dfs = bin_by_gc(df)
    binned_arrays = list()

    df = df[~df[col].isna()].copy()            # we only need non-NA
    df[f"CONSENSUS_{col}"] = df.apply(lambda r: get_consensus_sequence(r[col]), axis=1)
    df = df[df.groupby(f"CONSENSUS_{col}")[f"CONSENSUS_{col}"].transform(len) > 5]

    array, update_shifts = create_numpy_for_column_with_extended_motif(env, df, col)

    binned_arrays.append({
        "GC": df["GC"],
        "motifs": array,
        "shifts": update_shifts
        })

    example = df.at[df.index[0], col]       # type: Dict[str, List[float]]
    w = len(next(iter(example.values())))  # width (numbere of positions)
    b = len(example)  # number of bases (letters)

    for w_pos in range(array.shape[1]):
        letters = example.keys()
        letter_to_idx = {x: x_pos for x_pos, x in enumerate(sorted(letters))}

        fig, axes = plt.subplots(2, math.ceil(len(letters) / 2), sharex="all", sharey="all")
        ax_counter = 0
        for l, ax in zip(letters, axes.ravel()[:len(letters)]):
            # go through GC bins
            all_gc = list()
            all_probs = list()

            for ba in binned_arrays:
                arr = ba["motifs"]
                gc = ba["GC"].values
                shifts = ba["shifts"]

                for index in range(len(shifts)):

                    shifted_position = w_pos
                    # print(w_pos, shifted_position)

                    # shifted_pos = w_pos - shifts[index]
                    # if shifted_pos < 0 or shifted_pos >= w:
                    #     continue
                    if w_pos < shifts[index] or w_pos >= shifts[index] + 6:
                        continue

                    all_gc.append(gc[index])
                    if arr[index, shifted_position, letter_to_idx[l]] < 0 or arr[index, shifted_position, letter_to_idx[l]] > 1:
                        raise ValueError("Something's up")
                    all_probs.append(arr[index, shifted_position, letter_to_idx[l]])

            # ax.scatter(all_gc, all_probs, marker="+")
            # seaborn.regplot(all_gc, all_probs, ax=ax, lowess=True, scatter_kws={"s": 5, "alpha": 0.3})
            ax.set_title(f"{l, w_pos}")

            df = pd.DataFrame({"GC": all_gc, "Probability": all_probs})
            df.sort_values("GC", inplace=True)
            loess_with_stde(df, "GC", "Probability", ax, None)

            ax_counter += 1
            # plt.show()
        plt.show()


def plot_letter_over_position(env, df, col, title=""):
    # type: (Environment, pd.DataFrame, str, str) -> None

    collect = dict()
    array, update_shifts = create_numpy_for_column_with_extended_motif(env, df, col,
                                                                       collect)
    df_original = df
    binned_arrays = [{
        "GC": df["GC"],
        "motifs": array,
        "shifts": update_shifts
    }]

    example = df.at[df.index[0], col]  # type: Dict[str, List[float]]
    w = len(next(iter(example.values())))  # width (numbere of positions)
    b = len(example)  # number of bases (letters)



    letters = example.keys()
    letter_to_idx = {x: x_pos for x_pos, x in enumerate(sorted(letters))}

    # fig, axes = plt.subplots(2, math.ceil(len(letters) / 2), sharex="all", sharey="all")
    fig = plt.figure(figsize=(10,12))
    shape = (4,2)

    ax1 = plt.subplot2grid(shape, (0, 0))
    ax2 = plt.subplot2grid(shape, (0, 1))
    ax3 = plt.subplot2grid(shape, (1, 0))
    ax4 = plt.subplot2grid(shape, (1, 1))
    ax_logo = plt.subplot2grid(shape, (3, 0))
    ax_counts = plt.subplot2grid(shape, (2,0))
    ax_pos_dist = plt.subplot2grid(shape, (2, 1))
    ax_text = plt.subplot2grid(shape, (3,1))


    axes = [ax1, ax2, ax3, ax4]

    # for each letter
    # for l, ax in zip(letters, axes.ravel()[:len(letters)]):
    ylim=[-0.1, 1.1]
    for l, ax in zip(letters, axes):
        # for each position in motif
        # go through df and accumulate values
        all_gc = list()
        all_probs = list()
        for w_pos in range(array.shape[1]):

            for ba in binned_arrays:
                arr = ba["motifs"]
                gc = ba["GC"].values
                shifts = ba["shifts"]

                for index in range(len(shifts)):

                    shifted_position = w_pos
                    # print(w_pos, shifted_position)

                    # shifted_pos = w_pos - shifts[index]
                    # if shifted_pos < 0 or shifted_pos >= w:
                    #     continue
                    if w_pos < shifts[index] or w_pos >= shifts[index] + 6:
                        continue

                    all_gc.append(shifted_position)

                    if arr[index, shifted_position, letter_to_idx[l]] < 0 or arr[
                        index, shifted_position, letter_to_idx[l]] > 1:
                        raise ValueError("Something's up")
                    all_probs.append(arr[index, shifted_position, letter_to_idx[l]])

            # ax.scatter(all_gc, all_probs, marker="+")
            # seaborn.regplot(all_gc, all_probs, ax=ax, lowess=True, scatter_kws={"s": 5, "alpha": 0.3})
        ax.set_title(f"{l}")

        df = pd.DataFrame({"Position": all_gc, "Probability": all_probs})
        df.sort_values("Position", inplace=True)

        # seaborn.kdeplot(df["Position"], df["Probability"], cmap="Reds", ax=ax)


        df_mean = df.groupby("Position", as_index=False).mean()
        seaborn.boxplot("Position", "Probability", data=df, ax=ax, color="red", fliersize=0)
        seaborn.lineplot(df_mean["Position"], df_mean["Probability"], ax=ax, color="blue")
        ax.set_ylim(ylim)
        # loess_with_stde(df, "Position", "Probability", ax, None)

            # plt.show()

    # add logo
    ax=ax_logo
    msa_t = collect["msa_t"]
    seqs = [x.seq._data for x in msa_t.list_alignment_sequences]
    counts_mat = lm.alignment_to_matrix(sequences=seqs, to_type='counts', characters_to_ignore='.-X')

    # Counts matrix -> Information matrix
    info_mat = lm.transform_matrix(counts_mat,
                                   from_type='counts',
                                   to_type='information')

    lm.Logo(info_mat, ax=ax, color_scheme="classic")
    ax.set_ylim([0,2])


    # add distplot of starting positions
    ax = ax_counts
    # seaborn.distplot(update_shifts, ax=ax)
    counter = Counter(update_shifts)
    total = sum(counter.values())
    to_add = sorted(set(range(4)).difference(counter.keys()))
    normalized = [[x, 100*counter[x]/total] for x in counter] + [[x, 0] for x in to_add]
    normalized = np.array(normalized)
    seaborn.barplot(normalized[:, 0], normalized[:, 1], ax=ax, color="blue")
    ax.set_ylim([0, 100])
    ax.set_ylabel("Probability")
    ax.set_xlabel("Shift in consensus")

    ### Plot position distribution
    col_pos = col.replace("_MAT", "_POS_DISTR")
    ax = ax_pos_dist
    shift_to_pos_dist = get_position_distributions_by_shift(df_original, col_pos, update_shifts)
    for s in sorted(shift_to_pos_dist.keys()):
        list_pos_dist = shift_to_pos_dist[s]

        # average positions
        values = dict()
        for l in list_pos_dist:
            try:
                for i in l.keys():
                    if i not in values.keys():
                        values[i] = list()
                    values[i].append(l[i])
            except Exception:
                continue
        for i in values.keys():
            values[i] = np.mean(values[i])

        total = sum(values.values())
        for i in values.keys():
            values[i] /= total

        x = sorted(values.keys())
        y = [values[a] for a in x]

        seaborn.lineplot(x,y, label=s, ax=ax)

    ax.legend()


    # TEXT
    ax = ax_text
    from matplotlib.font_manager import FontProperties
    fp = FontProperties()
    fp.set_family("monospace")
    print("here")
    print(print_reduced_msa(msa_t, True, n=10))
    ax.text(0, 0, print_reduced_msa(msa_t, True, n=10),
            horizontalalignment='left',
            verticalalignment='center',
            fontproperties=fp)
    ax.set_xlim([-0.2, 0.4])
    ax.set_ylim([-0.4, 0.4])
    # ax.axis("off",)
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)

    plt.suptitle("Gc range: {}. Num Data points: {}".format(title, msa_t.number_of_sequences()))
    # save_figure(FigureOptions(save_fig=next_name(env["pd-work"])))
    plt.tight_layout()
    plt.subplots_adjust(top=0.9)

    plt.savefig(next_name(env["pd-work"]))
    plt.show()

def plot_letter_per_gc_over_position_with_alignment_per_gc(env, df, col):
    # type: (Environment, pd.DataFrame, str) -> None
    df = df[~df[col].isna()].copy()  # we only need non-NA
    binned_dfs = bin_by_gc(df, step=5)
    binned_arrays = list()
    # for each bin
    # for info in binned_dfs:
    #     lower, upper, df_gc = info
    #     df_gc = df_gc.copy()
    #     if len(df_gc) <= 1:
    #         continue
    #
    #     array, update_shifts = create_numpy_for_column_and_realign(env, df_gc, col)
    #
    #     binned_arrays.append({
    #         "GC": df_gc["GC"],
    #         "motifs": array,
    #         "shifts": update_shifts
    #     })
    #
    #     df_gc[f"CONSENSUS_{col}"] = df_gc.apply(lambda r: get_consensus_sequence(r[col]), axis=1)
    #     df_gc = df_gc[df_gc.groupby(f"CONSENSUS_{col}")[f"CONSENSUS_{col}"].transform(len) > 5]

    example = df.at[df.index[0], col]  # type: Dict[str, List[float]]
    w = len(next(iter(example.values())))  # width (numbere of positions)
    b = len(example)  # number of bases (letters)

    for info_gc in binned_dfs:
        lower, upper, df_gc = info_gc
        df_gc = df_gc.copy()

        if len(df_gc) <= 1:
            continue

        df_gc[f"CONSENSUS_{col}"] = df_gc.apply(lambda r: get_consensus_sequence(r[col]), axis=1)
        df_gc = df_gc[df_gc.groupby(f"CONSENSUS_{col}")[f"CONSENSUS_{col}"].transform(len) > 5]

        if len(df_gc) <= 1:
            continue

        print(f"GC in [{lower}, {upper}]")
        plot_letter_over_position(env, df_gc, col, f"[{lower}, {upper}]")


def plot_candidate_codons(env, df, codons):
    # type: (Environment, pd.DataFrame, List[str]) -> None

    fig, ax = plt.subplots()

    for c in sorted(codons):
        seaborn.regplot(df["GC"].astype(float).values, df[c].astype(float).values, label=c,
                        lowess=True, scatter_kws={"s": 5, "alpha": 0.1})

    ax.set_ylim([-0.05,1.05])
    ax.set_ylabel("Probability")
    ax.set_xlabel("GC")
    leg = ax.legend()
    for lh in leg.legendHandles:
        lh.set_alpha(1)

    plt.show()

def plot_candidate_starts(env, df):
    # type: (Environment, pd.DataFrame) -> None
    plot_candidate_codons(env, df, ["ATG", "GTG", "TTG"])

def plot_candidate_stops(env, df):
    # type: (Environment, pd.DataFrame) -> None
    plot_candidate_codons(env, df, ["TAA", "TAG", "TGA"])



def main(env, args):
    # type: (Environment, argparse.Namespace) -> None
    df_bac = load_obj(args.pf_input_bac)        # type: pd.DataFrame
    df_arc = load_obj(args.pf_input_arc)        # type: pd.DataFrame
    df_bac["Type"] = "Bacteria"
    df_arc["Type"] = "Archaea"

    df = pd.concat([df_bac, df_arc], sort=False)
    # df = df.sample(100)
    df["GENOME_TYPE"] = df["GENOME_TYPE"].apply(lambda x: x.strip().split("-")[1].upper())
    df.loc[df["GENOME_TYPE"] == "D2", "GENOME_TYPE"] = "D"

    df.reset_index(inplace=True)
    import matplotlib
    matplotlib.rcParams.update({
        # "pgf.texsystem": "pdflatex",
        'font.family': 'serif',
        'text.usetex': False,
        'pgf.rcfonts': False,
    })

    visualize_matrix_column(env, df, "RBS_MAT")
    visualize_matrix_column(env, df[(df["Type"] == "Bacteria") & (df["GENOME_TYPE"] == "C")], "PROMOTER_MAT")
    #
    # build_consensus_from_consensus(env, df[df["GENOME_TYPE"] == "A"] , "RBS_MAT")
    # build_consensus_from_consensus(env, df[df["GENOME_TYPE"] == "B"], "RBS_MAT")
    # build_consensus_from_consensus(env, df[df["GENOME_TYPE"] == "C"], "RBS_MAT")
    # build_consensus_from_consensus(env, df[df["GENOME_TYPE"] == "D"], "RBS_MAT")
    #
    #
    # build_consensus_from_consensus(env, df[df["GENOME_TYPE"] == "C"], "PROMOTER_MAT")
    # build_consensus_from_consensus(env, df[df["GENOME_TYPE"] == "D"], "PROMOTER_MAT")
    #
    #
    # plot_letter_over_gc(env, df, "RBS_MAT")
    # plot_letter_over_gc_with_alignment_per_bin(env, df[df["GENOME_TYPE"] == "A"], "RBS_MAT")
    # plot_letter_over_gc_with_single_alignment(env,  df[df["GENOME_TYPE"] != "D"], "RBS_MAT")
    # plot_letter_over_gc_with_single_alignment(env, df[df["GENOME_TYPE"] == "C"], "PROMOTER_MAT")

    # plot_letter_per_gc_over_position_with_alignment_per_gc(env,df[df["GENOME_TYPE"].isin({"A", "C"})], "RBS_MAT")
    # plot_letter_per_gc_over_position_with_alignment_per_gc(env, df[df["GENOME_TYPE"].isin({"C"})], "PROMOTER_MAT")
    #
    # plot_candidate_starts(env, df)
    # plot_candidate_stops(env, df)

    # for g in ["A", "B", "C"]:
    #     plot_candidate_starts(env, df[df["GENOME_TYPE"] == g])


if __name__ == "__main__":
    main(my_env, parsed_args)
