# Karl Gemayel
# Georgia Institute of Technology
#
# Created: 3/17/20
import logging
import argparse
import math

import pandas as pd
import numpy as np
from typing import *
import matplotlib.pyplot as plt
import logomaker as lm

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import sbsp_log  # runs init in sbsp_log and configures logger

# Custom imports
from sbsp_general import Environment

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from sbsp_general.GMS2Noncoding import GMS2Noncoding
from sbsp_general.shelf import next_name
from sbsp_io.objects import load_obj

parser = argparse.ArgumentParser("Description of driver.")

parser.add_argument('--pf-data', required=True)
parser.add_argument('--motif-type', required=True, choices=["RBS", "PROMOTER"])
parser.add_argument('--group', choices=["group-a", "group-b", "group-c"], nargs="+", required=True)

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


def motif_dict_to_df(motif_dict):
    # type: (Dict[str, List[float]]) -> pd.DataFrame

    header = sorted(motif_dict.keys())

    list_entries = list()
    motif_length = len(next(iter(motif_dict.values())))

    for i in range(motif_length):
        entry = {
            x: motif_dict[x][i] for x in header
        }

        list_entries.append(entry)

    return pd.DataFrame(list_entries)

def get_models_by_gc(df, gc_values, motif_type):
    # type: (pd.DataFrame, List[float], str) -> List[[pd.DataFrame,np.ndarray]]

    df.sort_values("GC", inplace=True)

    cpi = 0
    result = list()
    for i, gc in enumerate(gc_values):

        while cpi < len(df) and df.at[df.index[cpi], "GC"] < gc:
            cpi += 1
        cpi += 2

        if cpi >= len(df):
            break

        if i < len(gc_values) - 1 and df.at[df.index[cpi], "GC"] >= gc_values[i+1]:
            result.append(None)
            continue

        print(df.at[df.index[cpi], "GC"], df.at[df.index[cpi], "Genome"])
        result.append([motif_dict_to_df(df.at[df.index[cpi], f"{motif_type}_MAT"]),
                      GMS2Noncoding(df.at[df.index[cpi], "NON_MAT"]).pwm_to_array(0)])
                        # [0.25]*4])

    return result


def background_from_gc(gc):
    # type: (float) -> List[float]
    if gc > 1:
        gc /= 100.0

    at = 1-gc
    g = c = gc / 2.0
    a = at / 2.0
    t = 1 - a - g - c
    return [a, c, g, t]

def main(env, args):
    # type: (Environment, argparse.Namespace) -> None
    df_bac = load_obj(args.pf_data).reset_index()        # type: pd.DataFrame
    df_bac = df_bac[df_bac["GENOME_TYPE"].isin(args.group)]
    min_gc = 20
    max_gc = 70

    if args.motif_type == "PROMOTER":
        df_bac = df_bac[df_bac["GC"] >= 40].copy()

    gc_values = np.arange(min_gc, max_gc, 2)
    models = get_models_by_gc(df_bac, gc_values, motif_type=args.motif_type)

    num_plots = len(models)
    num_rows = int(math.sqrt(num_plots))
    num_cols = math.ceil(num_plots / float(num_rows))

    fig, axes = plt.subplots(num_rows, num_cols, sharex="all", sharey="all", figsize=(12,10))

    model_index = 0
    for r in range(num_rows):
        for c in range(num_cols):
            if model_index >= len(models):
                break

            if models[model_index] is None:
                model_index += 1
                continue

            bgd = [0.25]*4
            bgd = background_from_gc(gc_values[model_index])

            newmod = lm.transform_matrix(
                models[model_index][0], to_type="information", from_type="probability",
                background=models[model_index][1]
            )
            # from copy import copy
            # newmod = copy(models[model_index][0])
            # for idx in newmod.index:
            #     # see https://bioconductor.org/packages/release/bioc/vignettes/universalmotif/inst/doc/IntroductionToSequenceMotifs.pdf
            #
            #     uncertainty = sum(
            #         [newmod.at[idx, l] * math.log2(newmod.at[idx, l]) for l in newmod.columns]
            #     )
            #     fIC = math.log2(4) - uncertainty
            #     for i, l in enumerate(sorted(newmod.columns)):
            #         newmod.at[idx, l] = max(1 * newmod.at[idx, l] * math.log2(newmod.at[idx, l] / models[model_index][1][i]), 0)
            lm.Logo(newmod, ax=axes[r][c])

            axes[r][c].set_ylim(0, 2)
            axes[r][c].set_title(int(gc_values[model_index]))
            # fig.show()
            model_index += 1


    plt.tight_layout()
    plt.savefig(next_name(env["pd-work"]))
    plt.show()





if __name__ == "__main__":
    main(my_env, parsed_args)
