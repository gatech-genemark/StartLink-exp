import copy
import logging
from typing import *

import matplotlib.pyplot as plt
from matplotlib_venn import venn3

from sbsp_general.labels import create_key_3prime_from_label, create_gene_key_from_label, Labels
from sbsp_viz.general import FigureOptions

logger = logging.getLogger(__name__)


def get_set_3prime_keys(labels):
    # type: (Labels) -> Set[str]
    return set(create_key_3prime_from_label(l) for l in labels)


def get_set_gene_keys(labels):
    # type: (Labels) -> Set[str]
    return set(create_gene_key_from_label(l) for l in labels)


def get_intersection(list_set):
    # type: (List[Set]) -> Set
    if len(list_set) == 0:
        return set()

    acc = copy.copy(list_set[0])

    for i in range(1, len(list_set)):
        acc.intersection_update(list_set[i])

    return acc


def reduce_labels_to_genes_in_all(list_labels):
    # type: (List[Labels]) -> List[Labels]

    common_3prime_keys = get_intersection([get_set_3prime_keys(ls) for ls in list_labels])

    # extract labels with common 3prime
    list_labels_common_3prime = [
        l.get_multiple_by_3prime_keys(common_3prime_keys) for l in list_labels
    ]

    return list_labels_common_3prime

def numbers_for_3d_venn(labels_a, labels_b, labels_c):
    # type: (Labels, Labels, Labels) -> Dict[str, int]

    keys_a = get_set_gene_keys(labels_a)
    keys_b = get_set_gene_keys(labels_b)
    keys_c = get_set_gene_keys(labels_c)

    keys_ab = keys_a.intersection(keys_b)
    keys_bc = keys_b.intersection(keys_c)
    keys_ac = keys_a.intersection(keys_c)
    keys_abc = keys_ab.intersection(keys_ac).intersection(keys_bc)

    n_abc = len(keys_abc)
    n_ab = len(keys_ab) - n_abc
    n_ac = len(keys_ac) - n_abc
    n_bc = len(keys_bc) - n_abc

    n_a = len(keys_a) - (n_ab + n_ac + n_abc)
    n_b = len(keys_b) - (n_ab + n_bc + n_abc)
    n_c = len(keys_c) - (n_ac + n_bc + n_abc)


    # Make a Basic Venn
    label_value_pair = {
        "100": n_a,
        "010": n_b,
        "001": n_c,
        "110": n_ab,
        "101": n_ac,
        "011": n_bc,
        "111": n_abc
    }

    return label_value_pair


def venn_diagram_5prime(labels_a, labels_b, labels_c, figure_options=None):
    # type: (Labels, Labels, Labels, FigureOptions) -> None

    # first, reduce each set to common genes
    list_labels_common_3prime = reduce_labels_to_genes_in_all([labels_a, labels_b, labels_c])
    label_value_pair = numbers_for_3d_venn(*list_labels_common_3prime)


    fig, ax = plt.subplots()

    # venn3([set(get_set_gene_keys(labels)) for labels in list_labels_common_3prime],
    #       set_labels=[labels.name for labels in list_labels_common_3prime])

    # create equal sized circles
    v = venn3([1, 1, 1, 1, 1, 1, 1],
          set_labels=[labels.name for labels in list_labels_common_3prime])

    for key, value in label_value_pair.items():
        v.get_label_by_id(key).set_text(value)



    # Add title and annotation
    FigureOptions.set_properties_for_axis(ax, figure_options)

    if figure_options is not None and figure_options.save_fig is not None:
        plt.savefig(figure_options.save_fig, bbox_inches='tight')

    # Show it
    plt.show()
