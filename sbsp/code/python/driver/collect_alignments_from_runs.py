# Karl Gemayel
# Georgia Institute of Technology
#
# Created: 3/17/20
import logging
import argparse
import shutil

import pandas as pd
from typing import *

# noinspection All
from tqdm import tqdm

import pathmagic

# noinspection PyUnresolvedReferences
import sbsp_log  # runs init in sbsp_log and configures logger

# Custom imports
from sbsp_container.genome_list import GenomeInfoList, GenomeInfo
from sbsp_general import Environment

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from sbsp_general.general import os_join
from sbsp_general.labels_comparison_detailed import LabelsComparisonDetailed
from sbsp_general.shelf import create_q_key_3p, add_q_key_3p_to_df
from sbsp_io.general import mkdir_p
from sbsp_io.labels import read_labels_from_file

parser = argparse.ArgumentParser("Description of driver.")

parser.add_argument('--pf-genome-list', required=True, help="List of genomes")

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


def collect_alignments_for_genome(env, gi):
    # type: (Environment, GenomeInfo) -> None
    pd_genome = os_join(env["pd-work"], gi.name)

    mkdir_p(pd_genome)

    pd_run = os_join(env["pd-runs"], gi.name)

    # load labels and data files
    pf_sbsp = os_join(pd_run, "sbsp", "accuracy", f"{gi.name}.gff")
    pf_gms2 = os_join(pd_run, "gms2", "gms2.gff")
    pf_ncbi = os_join(pd_run, "ncbi", "ncbi.gff")
    pf_sbsp_details = os_join(pd_run, "sbsp", "output.csv")

    common_options = {
        "ignore_frameshifted": True, "ignore_partial": True, "shift": 0
    }

    try:

        labels_sbsp = read_labels_from_file(pf_sbsp, name="SBSP", **common_options)
        labels_gms2 = read_labels_from_file(pf_gms2, name="GMS2", **common_options)
        labels_ncbi = read_labels_from_file(pf_ncbi, name="NCBI", **common_options)
        df_details = pd.read_csv(pf_sbsp_details)
        add_q_key_3p_to_df(df_details, "q-3prime")
    except FileNotFoundError:
        return

    # get genes where GMS2=SBSP
    lcd_full = LabelsComparisonDetailed(labels_gms2, labels_sbsp,
                                        name_a="gms2", name_b="sbsp")

    labels_gms2_eq_sbsp = lcd_full.match_3p_5p("a")

    # get labels where gms2_eq_sbsp doesn't match NCBI
    lcd2 = LabelsComparisonDetailed(labels_gms2_eq_sbsp, labels_ncbi, name_a="gms2_eq_sbsp", name_b="ncbi")
    labels_gms2_eq_sbsp_not_ncbi = lcd2.match_3p_not_5p("a")

    # get msa files for all these labels
    set_3prime_keys = {create_q_key_3p(l.seqname(), l.left(), l.right(), l.strand()) for l in labels_gms2_eq_sbsp_not_ncbi}

    df_gms2_eq_sbsp_not_ncbi = df_details[df_details["q-3prime"].isin(set_3prime_keys)]

    set_pf_msa_out = set(df_gms2_eq_sbsp_not_ncbi["pf-msa-output"])

    for pf_msa_out in set_pf_msa_out:
        shutil.copy(pf_msa_out, pd_genome)







def main(env, args):
    # type: (Environment, argparse.Namespace) -> None
    gil = GenomeInfoList.init_from_file(args.pf_genome_list)

    for gi in tqdm(gil, total=len(gil)):
        collect_alignments_for_genome(env, gi)



if __name__ == "__main__":
    main(my_env, parsed_args)
