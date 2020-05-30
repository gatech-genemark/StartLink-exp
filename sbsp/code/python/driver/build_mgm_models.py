# Karl Gemayel
# Georgia Institute of Technology
#
# Created: 3/17/20
import logging
import argparse
import numpy as np
import pandas as pd
from typing import *

# noinspection All
import pathmagic


# noinspection PyUnresolvedReferences
import sbsp_log  # runs init in sbsp_log and configures logger

# Custom imports
from sbsp_general import Environment

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from sbsp_general.MGMMotifModel import MGMMotifModel
from sbsp_general.general import get_value
from sbsp_general.mgm_motif_model_all_gc import MGMMotifModelAllGC
from sbsp_general.shelf import bin_by_gc, get_consensus_sequence, create_numpy_for_column_with_extended_motif, \
    get_position_distributions_by_shift
from sbsp_io.objects import load_obj, save_obj
from sbsp_viz.mgm_motif_model import MGMMotifModelVisualizer

parser = argparse.ArgumentParser("Build MGM start models.")

parser.add_argument('--pf-input-arc', required=True, help="Input file")
parser.add_argument('--pf-input-bac', required=True, help="Input file")

parser.add_argument('--pf-output', required=True)

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

def fix_genome_type(df):
    # type: (pd.DataFrame) -> None
    df["GENOME_TYPE"] = df["GENOME_TYPE"].apply(lambda x: x.strip().split("-")[1].upper())
    df.loc[df["GENOME_TYPE"] == "D2", "GENOME_TYPE"] = "D"

def read_archaea_bacteria_inputs(pf_arc, pf_bac):
    # type: (str, str) -> pd.DataFrame
    df_bac = load_obj(pf_bac)  # type: pd.DataFrame
    df_arc = load_obj(pf_arc)  # type: pd.DataFrame
    df_bac["Type"] = "Bacteria"
    df_arc["Type"] = "Archaea"

    df = pd.concat([df_bac, df_arc], sort=False)
    df.reset_index(inplace=True)
    fix_genome_type(df)

    return df

def mat_to_dict(mat):
    # type: (np.ndarray) -> Dict[str, List[float]]

    index_to_letter = {
        i: x for i, x in enumerate(list("ACGT"))
    }

    result = dict()

    P, L = mat.shape
    for l in range(L):
        l_str = index_to_letter[l]
        result[l_str] = list(mat[:, l])

    return result

def build_mgm_motif_model_for_gc(env, df, col, **kwargs):
    # type: (Environment, pd.DataFrame, str, Dict[str, Any]) -> MGMMotifModel

    min_consensus_occurence = get_value(kwargs, "min_consensus_occurence", 5)
    title = get_value(kwargs, "title", "")

    # filter out consensus sequences that don't appear very frequently
    df = df[df.groupby(f"CONSENSUS_{col}")[f"CONSENSUS_{col}"].transform(len) > min_consensus_occurence]

    if len(df) <= 1:
        return

    # run alignment of consensus sequences and get back shifts
    collect = dict()
    array, update_shifts = create_numpy_for_column_with_extended_motif(env, df, col,
                                                                       collect)



    # two widths:
    original_width = len(df.iloc[0][f"CONSENSUS_{col}"])
    extended_width = array.shape[1]

    # array has shape: N x P x L
    # N: Number of motifs
    # P: Width of extended motif
    # L: number of letters
    N, P, L = array.shape

    extended_motif = np.sum(array, 0)
    extended_motif = np.divide(
        extended_motif,
        extended_motif.sum(1).reshape(P, 1)
    )

    extended_motif_dict = mat_to_dict(extended_motif)

    # get prior probabilities on shift position
    counter = Counter(update_shifts)
    total = sum(counter.values())
    to_add = sorted(set(range(max(counter.keys()))).difference(counter.keys()))
    normalized = [[x, 100 * counter[x] / total] for x in counter] + [[x, 0] for x in to_add]
    normalized = np.array(normalized)
    shifts_dict = {normalized[x, 0]: normalized[x, 1] for x in range(normalized.shape[0])}

    # get position distributions
    col_pos = col.replace("_MAT", "_POS_DISTR")
    shift_to_pos_dist = get_position_distributions_by_shift(df, col_pos, update_shifts)

    position_distributions_by_shift = dict()        # type: Dict[int, Dict[int, float]]
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

        position_distributions_by_shift[s] = {
            a: b for a, b in zip(x,y)
        }

    # compile into single model

    mgm_mm = MGMMotifModel(shifts_dict, extended_motif_dict, original_width, position_distributions_by_shift)

    MGMMotifModelVisualizer.visualize(mgm_mm, title=title)
    return mgm_mm


def build_mgm_motif_models_for_all_gc(env, df, name, **kwargs):
    # type: (Environment, pd.DataFrame, str, Dict[str, Any]) -> MGMMotifModelAllGC
    df = df[~df[name].isna()].copy()  # we only need non-NA

    bin_size = get_value(kwargs, "bin_size", 5, default_if_none=True)

    # get consensus sequences for all motifs
    df[f"CONSENSUS_{name}"] = df.apply(lambda r: get_consensus_sequence(r[name]), axis=1)

    # bin dataframes by GC
    binned_dfs = bin_by_gc(df, step=bin_size)

    # for each binned dataframe, build specific model
    list_mgm_models = list() # type: List[Tuple[float, float, MGMMotifModel]]
    for info in binned_dfs:
        lower, upper, df_gc = info



        mgm_mm = None
        if len(df_gc) > 1:
            mgm_mm = build_mgm_motif_model_for_gc(env, df_gc, name, title=f"[{lower},{upper}]", **kwargs)

        if lower == 30 and upper == 35 and mgm_mm is None:
            print('hi')

        if mgm_mm is None:
            # use previous model
            if len(list_mgm_models) > 0:
                list_mgm_models.append([lower, upper, mgm_mm])
        else:
            list_mgm_models.append([lower, upper, mgm_mm])

    return MGMMotifModelAllGC(list_mgm_models)


def build_mgm_models(env, df, pf_output):
    # type: (Environment, pd.DataFrame, str) -> None

    # build the following models
    # {
    #   Bacteria: {
    #       RBS: {
    #           A,C: models_gc_rbs_ac
    #           B: models_gc_rbs_b
    #       }, PROMOTER: {
    #           C: models_gc_promoter_c
    #   },
    #   Archaea: {
    #       RBS: {
    #           A,D: models_gc_rbs_ad
    #       },
    #       PROMOTER: {
    #           D: models_gc_promoter_d
    #       }
    #   },

    type_model_group = {
        "Bacteria": {
            "RBS": {
                "AC", "B"
            },
            "PROMOTER": {
                "C"
            }
        },
        "Archaea": {
            "RBS": {
                "AD"
            },
            "PROMOTER": {
                "D"
            }
        }
    }

    name_to_models = dict()     # type: Dict[str, Dict[str, Dict[str, MGMMotifModelAllGC]]]
    for species_type in type_model_group.keys():
        name_to_models[species_type] = dict()       # type: Dict[str, Dict[str, MGMMotifModelAllGC]]
        for name in type_model_group[species_type].keys():
            name_to_models[species_type][name] = dict()
            for group in type_model_group[species_type][name]:
                if group != "C":
                    continue
                name_to_models[species_type][name][group] = build_mgm_motif_models_for_all_gc(
                    env, df[(df["Type"] == species_type) & (df["GENOME_TYPE"].isin(set(group)))], name + "_MAT"
                )

    save_obj(name_to_models, pf_output)


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None
    df = read_archaea_bacteria_inputs(args.pf_input_arc, args.pf_input_bac)

    build_mgm_models(env, df, args.pf_output)





if __name__ == "__main__":
    main(my_env, parsed_args)
