# Karl Gemayel
# Georgia Institute of Technology
#
# Created: 3/17/20
import logging
import argparse
import pandas as pd
from typing import *

# noinspection All
from tqdm import tqdm

import pathmagic

# noinspection PyUnresolvedReferences
import sbsp_log  # runs init in sbsp_log and configures logger

# Custom imports
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

def run_mgm_models_on_test_data(env, mgm_models, df_test, species_type, pf_output, **kwargs):
    # type: (Environment, Dict[str, Dict[str, Dict[str, MGMMotifModelAllGC]]], pd.DataFrame, str, str, Dict[str, Any]) -> None

    for idx in tqdm(df_test.index, total=len(df_test)):

        for name in mgm_models[species_type].keys():
            for group, model in mgm_models[species_type][name].items():
                tag = f"MGM_GC_{name}_{group}"

                # print(tag, df_test.at[idx, "Genome GC"])
                result = model.get_model_by_gc(df_test.at[idx, "Genome GC"]).find_best_position_and_score(
                    df_test.at[idx, "upstream_nt"]
                )

                df_test.at[idx, tag + "_score"] = result[1]
                df_test.at[idx, tag + "_position"] = result[0]


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None
    mgm_models = load_obj(args.pf_mgm_models)       # type: Dict[str, Dict[str, Dict[str, MGMMotifModelAllGC]]]
    df_test = pd.read_csv(args.pf_test)                # type: pd.DataFrame
    # df_test = df_test.sample(500)
    run_mgm_models_on_test_data(env, mgm_models, df_test, args.species_type, args.pf_output)
    df_test.to_csv(args.pf_output, index=False)

if __name__ == "__main__":
    main(my_env, parsed_args)
