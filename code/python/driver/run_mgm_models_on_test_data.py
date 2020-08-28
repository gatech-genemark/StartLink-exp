# Karl Gemayel
# Georgia Institute of Technology
#
# Created: 3/17/20
import logging
import argparse
import multiprocessing
from multiprocessing import Process

import numpy as np
import pandas as pd
from typing import *

# noinspection All
from tqdm import tqdm

import pathmagic

# noinspection PyUnresolvedReferences
import sbsp_log  # runs init in sbsp_log and configures logger

# Custom imports
from sbsp_general.general import get_value
from sbsp_io.objects import load_obj
from sbsp_general import Environment
from sbsp_general.mgm_motif_model_all_gc import MGMMotifModelAllGC


# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #


parser = argparse.ArgumentParser("Description of driver.")

parser.add_argument('--pf-test', required=True)
parser.add_argument('--pf-mgm-models', required=True, help="Pickled file of MGM Motif models")
parser.add_argument('--pf-output', required=True, help="Output file with all scores")
parser.add_argument('--species-type', required=True, choices=["Archaea", "Bacteria"])

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


def extract_motif(r, pos_col):
    ups = r["upstream_nt"]
    upstream_length = len(ups)
    mw = 6
    pos = r[pos_col]
    left = int(upstream_length - pos - mw)
    return ups[left: left + mw]


def run_mgm_models_on_test_data(env, mgm_models, df_test, species_type, pf_output, **kwargs):
    # type: (Environment, Dict[str, Dict[str, Dict[str, MGMMotifModelAllGC]]], pd.DataFrame, str, str, Dict[str, Any]) -> None
    for idx in tqdm(df_test.index, total=len(df_test)):

        results = list()

        for name in mgm_models[species_type].keys():
            for group, model in mgm_models[species_type][name].items():
                tag = f"MGM_GC_{name}_{group}"

                # print(tag, df_test.at[idx, "Genome GC"])
                model_gc = model.get_model_by_gc(df_test.at[idx, "Genome GC"])
                result = model_gc.find_best_position_and_score(
                    df_test.at[idx, "upstream_nt"]
                )

                df_test.at[idx, tag + "_score_noprior"] = result[2]
                df_test.at[idx, tag + "_score"] = result[1]
                df_test.at[idx, tag + "_position"] = len(df_test.at[idx, "upstream_nt"]) - result[0] - model_gc.motif_width()

                results.append([
                    name, group, result[1], result[2], df_test.at[idx, tag + "_position"]
                ])
                # if df_test.at[idx, "motif"] == "AGGAGG" and abs(df_test.at[idx, tag + "_position"] - df_test.at[idx, "RBS_both_position"]) > 5:
                #     print("AGGAGG", extract_motif(df_test.loc[idx], tag + "_position"), df_test.at[idx, tag+"_position"])
                #
                #     model.get_model_by_gc(df_test.at[idx, "Genome GC"]).find_best_position_and_score(
                #         df_test.at[idx, "upstream_nt"]
                #     )

        # get best noprior score across models
        best = max(
            results,
            key=lambda x: x[3])

        df_test.at[idx, "mgm_best_position"] = best[4]
        df_test.at[idx, "mgm_best_score"] = best[3]
        df_test.at[idx, "mgm_best_name"] = best[0]
        df_test.at[idx, "mgm_best_group"] = best[1]


    return df_test


def parallelize_dataframe(df, func):
    num_cores = multiprocessing.cpu_count()-1  #leave one free to not freeze machine
    num_partitions = num_cores #number of partitions to split dataframe
    df_split = np.array_split(df, num_partitions)
    pool = multiprocessing.Pool(num_cores)
    df = pd.concat(pool.map(func, df_split))
    pool.close()
    pool.join()
    return df


def _par_helper(func, func_kwargs, proc_num, output):
    # type: (Callable, Dict[str, Any], int, List) -> None
    output[proc_num] = func(**func_kwargs)


def parallelize_dataframe_by_chunks(df, func, data_name, func_kwargs, **kwargs):
    # type: (pd.DataFrame, Callable, str, Dict[str, Any], Dict[str, Any]) -> pd.DataFrame

    num_processors = get_value(kwargs, "num_processors", multiprocessing.cpu_count()-1)

    df_split = np.array_split(df, num_processors)
    num_partitions = len(df_split)
    logger.info(f"Running on {num_partitions} partitions: " + " ".join([str(len(x)) for x in df_split]))

    active_processes = list()

    # run separate process on each split
    manager = multiprocessing.Manager()
    return_dict = manager.dict()

    # create N jobs
    for n in range(num_partitions):
        # p = Process(target=func,
        #             kwargs={data_name: df_split[n], **func_kwargs}
        #             )
        p = Process(target=_par_helper,
                    args=(func, {data_name: df_split[n], **func_kwargs}, n, return_dict))

        p.start()
        logger.info(f"Started {n}")
        active_processes.append(p)


    for p in active_processes:
        p.join()

    # print(return_dict)
    return pd.concat(return_dict.values(), sort=False, ignore_index=True)


def run_mgm_models_on_test_data_parralel(env, mgm_models, df_test, species_type, pf_output, **kwargs):
    # type: (Environment, Dict[str, Dict[str, Dict[str, MGMMotifModelAllGC]]], pd.DataFrame, str, str, Dict[str, Any]) -> None


    parallelize_dataframe(df_test, run_mgm_models_on_test_data(env, mgm_models, ))
    for idx in tqdm(df_test.index, total=len(df_test)):

        results = list()

        for name in mgm_models[species_type].keys():
            for group, model in mgm_models[species_type][name].items():
                tag = f"MGM_GC_{name}_{group}"

                # print(tag, df_test.at[idx, "Genome GC"])
                model_gc = model.get_model_by_gc(df_test.at[idx, "Genome GC"])
                result = model_gc.find_best_position_and_score(
                    df_test.at[idx, "upstream_nt"]
                )

                df_test.at[idx, tag + "_score_noprior"] = result[2]
                df_test.at[idx, tag + "_score"] = result[1]
                df_test.at[idx, tag + "_position"] = len(df_test.at[idx, "upstream_nt"]) - result[0] - model_gc.motif_width()

                results.append([
                    name, group, result[1], result[2], df_test.at[idx, tag + "_position"]
                ])
                # if df_test.at[idx, "motif"] == "AGGAGG" and abs(df_test.at[idx, tag + "_position"] - df_test.at[idx, "RBS_both_position"]) > 5:
                #     print("AGGAGG", extract_motif(df_test.loc[idx], tag + "_position"), df_test.at[idx, tag+"_position"])
                #
                #     model.get_model_by_gc(df_test.at[idx, "Genome GC"]).find_best_position_and_score(
                #         df_test.at[idx, "upstream_nt"]
                #     )

        # get best noprior score across models
        best = max(
            results,
            key=lambda x: x[3])

        df_test.at[idx, "best_position"] = best[4]
        df_test.at[idx, "best_score"] = best[3]
        df_test.at[idx, "best_name"] = best[0]
        df_test.at[idx, "best_group"] = best[1]



def main(env, args):
    # type: (Environment, argparse.Namespace) -> None
    mgm_models = load_obj(args.pf_mgm_models)       # type: Dict[str, Dict[str, Dict[str, MGMMotifModelAllGC]]]
    df_test = pd.read_csv(args.pf_test)                # type: pd.DataFrame
    # df_test = df_test.head(500).copy()
    run_mgm_models_on_test_data(env, mgm_models, df_test, args.species_type, args.pf_output)
    # df_test = parallelize_dataframe_by_chunks(df_test, run_mgm_models_on_test_data, "df_test", {
    #     "env": env, "mgm_models": mgm_models, "species_type": args.species_type, "pf_output": args.pf_output
    # })


    # return
    df_test.to_csv(args.pf_output, index=False)

if __name__ == "__main__":
    main(my_env, parsed_args)
