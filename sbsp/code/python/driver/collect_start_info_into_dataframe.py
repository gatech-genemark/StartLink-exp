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
from sbsp_container.genome_list import GenomeInfoList, GenomeInfo
from sbsp_container.gms2_mod import GMS2Mod
from sbsp_general import Environment
import sbsp_argparse.parallelization
from sbsp_options.parallelization import ParallelizationOptions
from sbsp_parallelization.pbs import PBS
from sbsp_pbs_data.mergers import merge_identity
from sbsp_pbs_data.splitters import split_genome_info_list



# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from sbsp_general.general import os_join
from sbsp_general.shelf import compute_gc_from_file
from sbsp_io.objects import save_obj

parser = argparse.ArgumentParser("Description of driver.")

parser.add_argument('--pf-genome-list', required=True, help="Genome list")
parser.add_argument('--pf-output', required=True, help="Output file")

sbsp_argparse.parallelization.add_parallelization_options(parser)


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


def collect_start_info_from_gi(env, gi):
    # type: (Environment, GenomeInfo) -> Dict[str, Any]
    pd_genome = os_join(env["pd-data"], gi.name)
    pf_sequence = os_join(pd_genome, "sequence.fasta")

    gc = compute_gc_from_file(pf_sequence)

    pd_genome_run = os_join(env["pd-runs"], gi.name)
    pd_gms2 = os_join(pd_genome_run, "gms2")
    pf_mod = os_join(pd_gms2, "GMS2.mod")

    mod = GMS2Mod.init_from_file(pf_mod)

    return {
        "Genome": gi.name,
        "GC": 100*gc,
        **{
            x: mod.items[x] for x in {
                "GENOME_TYPE", "RBS_MAT", "RBS_MAT", "PROMOTER_MAT", "PROMOTER_WIDTH", "RBS_WIDTH",
                "RBS_POS_DISTR", "PROMOTER_POS_DISTR", "ATG", "GTG", "TTG", "TAA", "TGA", "TAG",
                "NON_MAT"
            } if x in mod.items.keys()
        }
    }

def collect_start_info_from_gil(env, gil, pf_output=None):
    # type: (Environment, GenomeInfoList) -> pd.DataFrame

    list_entries = list()
    for gi in tqdm(gil, total=len(gil)):
        entry = collect_start_info_from_gi(env, gi)
        list_entries.append(entry)

    return pd.DataFrame(list_entries)



def main(env, args):
    # type: (Environment, argparse.Namespace) -> None

    gil = GenomeInfoList.init_from_file(args.pf_genome_list)
    prl_options = ParallelizationOptions.init_from_dict(env, vars(args))


    if prl_options is not None:
        pbs = PBS(env, prl_options,
          splitter=split_genome_info_list,
          merger=merge_identity
          )

        output = pbs.run(
            data={"gil": gil},
            func=collect_start_info_from_gil,
            func_kwargs={
                "env": env,
            }
        )

        df = pd.concat(output, sort=False)
    else:
        df = collect_start_info_from_gil(env, gil)



    save_obj(df, args.pf_output)



if __name__ == "__main__":
    main(my_env, parsed_args)
