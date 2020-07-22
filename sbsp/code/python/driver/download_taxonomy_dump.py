# Karl Gemayel
# Georgia Institute of Technology
#
# Created: 3/17/20
import logging
import argparse
import urllib.request
from typing import *

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import sbsp_log  # runs init in sbsp_log and configures logger

# Custom imports
from sbsp_container.assembly_summary import AssemblySummary
from sbsp_container.taxonomy_tree import TaxonomyTree
from sbsp_general import Environment

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from sbsp_general.data_download import download_data_by_ancestor
from sbsp_general.general import os_join, run_shell_cmd
from sbsp_io.assembly_summary import read_assembly_summary_into_dataframe
from sbsp_io.general import mkdir_p

parser = argparse.ArgumentParser("Download taxonomy information from NCBI.")

parser.add_argument('--pd-output', required=True)
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


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None

    # link to taxonomy dump
    lp_taxonomy = f"https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip"

    pd_output = args.pd_output

    mkdir_p(pd_output)
    pf_output = os_join(pd_output, "taxdump.zip")

    logger.info(f"Downloading file: {lp_taxonomy}")
    urllib.request.urlretrieve(lp_taxonomy, pf_output)


    logger.info("Download complete. Unzipping")
    run_shell_cmd(f"cd {pd_output}; unzip {pf_output}")


if __name__ == "__main__":
    main(my_env, parsed_args)
