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

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from sbsp_general.general import os_join
from sbsp_general.labels import Labels
from sbsp_io.labels import read_labels_from_file
from sbsp_io.sequences import read_fasta_into_hash

parser = argparse.ArgumentParser("Description of driver.")

parser.add_argument('--pf-genome-list', required=True, help="Genome information list")
parser.add_argument('--pf-output', required=True, help="Output file")

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

def extract_upstream_sequence(label, sequences):
    # type: (Label, Dict

def create_mgm_test_data_for_genome(env, gi, **kwargs):
    # type: (Environment, GenomeInfo, Dict[str, Any]) -> pd.DataFrame
    pd_genome = os_join(env["pd-data"], gi.name)
    pf_sequence = os_join(pd_genome, "sequence.fasta")

    pd_genome_run = os_join(env["pd-runs"], gi.name)
    pf_gms2 = os_join(pd_genome_run, "gms2", "gms2.lst")
    pf_mod = os_join(pd_genome_run, "gms2", "GMS2.mod")

    labels = read_labels_from_lst_file(pf_gms2)     # type: Labels
    sequences = read_fasta_into_hash(pf_sequence)
    mod = GMS2Mod.init_from_file(pf_mod)

    # extract upstream regions
    list_entries = list()       # type: List[Dict[str, Any]]

    for l in labels:
        motif_type = l.get_attribute_value("motif-type")
        if motif_type == 1:     # RBS
            motif = mod.items["RBS_MAT"]
            motif_pos = mod.items["RBS_POS_DIST"]
        elif motif_type == 2:   # PROMOTER
            motif = mod.items["PROMOTER_MAT"]
            motif_pos = mod.items["PROMOTER_POS_DIST"]
        else:
            motif = None
            motif_pos = None

        # extract upstream sequence



    return pd.DataFrame(list_entries)



def create_mgm_test_data(env, gil, pf_output=None, **kwargs):
    # type: (Environment, GenomeInfoList, str, Dict[str, Any]) -> pd.DataFrame

    list_df = list()
    for gi in tqdm(gil, total=len(gil)):
        df_gi = create_mgm_test_data_for_genome(env, gi, **kwargs)
        list_df.append(df_gi)

    df = pd.concat(list_df, sort=False)
    return df

def main(env, args):
    # type: (Environment, argparse.Namespace) -> None
    gil = GenomeInfoList.init_from_file(args.pf_genome_list)
    df = create_mgm_test_data(env, gil)

    df.to_csv(args.pf_output)


if __name__ == "__main__":
    main(my_env, parsed_args)
