
import numpy as np
from Bio import pairwise2
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage
import operator

from sbsp_general.general import get_value



def select_alignment_with_smallest_number_of_gaps(alignments):

    if len(alignments) == 0:
        return None


    alignments_sorted = sorted(alignments, key=operator.itemgetter(0, 1))
    return alignments_sorted[0]

    smallest_num_gaps = float('inf')

    results = list()
    for a in alignments:
        curr_gaps = a[0].count("-") + a[1].count("-")

        if curr_gaps < smallest_num_gaps:
            smallest_num_gaps = curr_gaps
            results.append(a)

    result = None

    return result

def add_stop_codon_to_blosum(matrix, stopToAA=-4, stopToStop=1):
    # type: (dict, int, int) -> dict
    # get unique aa

    def add_letter_to_blossum(matrix, letter, to_other=-4, to_itself=1):
        # type: (dict, str, int, int) -> dict
        # get unique aa
        uniqueAA = set()
        for a in matrix:
            uniqueAA.add(a[0])
            uniqueAA.add(a[1])

        # for each, add penalty score
        for aa in uniqueAA:
            matrix[(aa, letter)] = to_other
            matrix[(letter, aa)] = to_itself

        matrix[(letter, letter)] = to_itself

    uniqueAA = set()
    for a in matrix:
        uniqueAA.add(a[0])
        uniqueAA.add(a[1])

    # for each, add penalty score
    for aa in uniqueAA:
        matrix[(aa, "*")] = stopToAA
        matrix[("*", aa)] = stopToAA

    matrix[("*", "*")] = stopToStop

    add_letter_to_blossum(matrix, "J")  # fixme: put it somewhere better

def count_differences_in_collapsed_alignment(seq1, seq2):

    if len(seq1) != len(seq2):
        raise ValueError("Sequences should have the same length")

    difference = 0          # number of non-gap mismatches
    length = 0              # collapsed length

    for i in range(len(seq1)):
        if seq1[i] != "-" and seq2[i] != "-":
            if seq1[i] != seq2[i]:
                difference += 1
            length += 1

    return {"difference": difference, "length": length}


def compute_ga_distance(seq1, seq2, matrix):

    # run Needleman-Wunsch
    nw = pairwise2.align.globaldx(seq1, seq2, matrix)

    if len(nw) == 0:
        return float('inf')

    diffs = count_differences_in_collapsed_alignment(nw[0][0], nw[0][1])
    distance = - np.log2(1 - diffs["difference"] / float(diffs["length"]))


    return distance


def compute_ga_distance_info(seq1, seq2, matrix):
    nw = pairwise2.align.globaldx(seq1, seq2, matrix)

    diffs = count_differences_in_collapsed_alignment(nw[0][0], nw[0][1])
    return diffs


def global_alignment(seq1, seq2, matrix):

    nw = pairwise2.align.globaldx(seq1, seq2, matrix)

    return nw[0][0], nw[0][1]

def global_alignment_nt(seq1, seq2):
    nw = pairwise2.align.globalxs(seq1, seq2, -4, -2)

    if len(nw) == 0:
        return ["", "", 0, 0, 0]

    # print nw
    return select_alignment_with_smallest_number_of_gaps(nw)

def global_alignment_aa(seq1, seq2, matrix):
    nw = pairwise2.align.globaldx(seq1, seq2, matrix)

    if len(nw) == 0:
        return ["", "", 0, 0, 0]

    return select_alignment_with_smallest_number_of_gaps(nw)

def global_alignment_aa_with_gap(seq1, seq2, matrix):
    nw = pairwise2.align.globalds(seq1, seq2, matrix, -4, -2)

    if len(nw) == 0:
        return ["", "", 0, 0, 0]

    return select_alignment_with_smallest_number_of_gaps(nw)


def local_alignment_aa(seq1, seq2, matrix, gap_open=-4, gap_extend=-2):
    nw = pairwise2.align.localds(seq1, seq2, matrix, gap_open, gap_extend)
    # nw = pairwise2.align.localdx(seq1, seq2, matrix)

    if len(nw) == 0:
        return ["", "", 0, 0, 0]

    return select_alignment_with_smallest_number_of_gaps(nw)




def compute_distances_for_list_of_pairs(genome_pairs, func, matrix, genome_ortholog_pair, sequences):
    distances = dict()
    for stp in genome_pairs:
        source_genome = stp["source"]
        target_genome = stp["target"]

        source_fasta_id = genome_ortholog_pair[source_genome]
        target_fasta_id = genome_ortholog_pair[target_genome]

        source_sequence = sequences[source_fasta_id]
        target_sequence = sequences[target_fasta_id]

        distance = func(source_sequence, target_sequence, matrix)

        if source_genome not in distances:
            distances[source_genome] = dict()
        if target_genome not in distances:
            distances[target_genome] = dict()

        distances[source_genome][target_genome] = distance
        distances[target_genome][source_genome] = distance
    return distances


def seriation(Z, N, cur_index):
    '''
        input:
            - Z is a hierarchical tree (dendrogram)
            - N is the number of points given to the clustering process
            - cur_index is the position in the tree for the recursive traversal
        output:
            - order implied by the hierarchical tree Z

        seriation computes the order implied by a hierarchical tree (dendrogram)
    '''
    if cur_index < N:
        return [cur_index]
    else:
        left = int(Z[cur_index - N, 0])
        right = int(Z[cur_index - N, 1])
        return (seriation(Z, N, left) + seriation(Z, N, right))



def compute_serial_matrix(dist_mat, method="average"):
    '''
        input:
            - dist_mat is a distance matrix
            - method = ["ward","single","average","complete"]
        output:
            - seriated_dist is the input dist_mat,
              but with re-ordered rows and columns
              according to the seriation, i.e. the
              order implied by the hierarchical tree
            - res_order is the order implied by
              the hierarhical tree
            - res_linkage is the hierarhical tree (dendrogram)

        compute_serial_matrix transforms a distance matrix into
        a sorted distance matrix according to the order implied
        by the hierarchical tree (dendrogram)
    '''
    N = len(dist_mat)
    flat_dist_mat = squareform(dist_mat)
    res_linkage = linkage(flat_dist_mat, method=method)
    res_order = seriation(res_linkage, N, N + N - 2)
    seriated_dist = np.zeros((N, N))
    a, b = np.triu_indices(N, k=1)
    seriated_dist[a, b] = dist_mat[[res_order[i] for i in a], [res_order[j] for j in b]]
    seriated_dist[b, a] = seriated_dist[a, b]

    return seriated_dist, res_order, res_linkage



def k2p_distance(seq_a, seq_b, **kwargs):
    # type: (str, str) -> float

    if len(seq_a) != len(seq_b):
        raise ValueError("Sequence sizes are not the same: {} != {}".format(len(seq_a), len(seq_b)))


    kimura_on_3rd = get_value(kwargs, "kimura_on_3rd", False)

    from math import log, sqrt

    alignment_length = len(seq_a)

    distance = 0

    # compute transversions and transitions
    transitions = {"AG", "GA", "CT", "TC"}
    transversions = {"AC", "CA", "AT", "TA",
                     "GC", "CG", "GT", "TG"}

    ts_count = 0            # transitions
    tv_count = 0            # transversions

    ungapped_length = 0

    begin = 0
    end = alignment_length
    step = 1

    if kimura_on_3rd:
        begin = 2
        step = 3

    for i in range(begin, end, step):

        # ignore gaps
        if seq_a[i] == "-" or seq_b[i] == "-":
            continue

        ungapped_length += 1
        combined = seq_a[i] + seq_b[i]

        if combined in transitions:
            ts_count += 1
        elif combined in transversions:
            tv_count += 1

    if ungapped_length == 0:
        return 0

    p = float(ts_count) / ungapped_length
    q = float(tv_count) / ungapped_length

    try:
        distance = -0.5 * log((1 - 2*p - q) * sqrt(1 - 2*q))
    except ValueError:
        raise ValueError("Can't take log of negative value")

    return distance
