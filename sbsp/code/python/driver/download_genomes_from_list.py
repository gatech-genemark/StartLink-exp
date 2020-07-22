# Karl Gemayel
# Georgia Institute of Technology
#
# Created: 3/17/20
import logging
import argparse
from typing import *

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import sbsp_log  # runs init in sbsp_log and configures logger

# Custom imports
from sbsp_container.assembly_summary import AssemblySummary
from sbsp_container.genome_list import GenomeInfoList, GenomeInfo
from sbsp_container.taxonomy_tree import TaxonomyTree
from sbsp_general import Environment

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from sbsp_general.data_download import download_data_by_ancestor, download_data_from_assembly_summary
from sbsp_io.assembly_summary import read_assembly_summary_into_dataframe

parser = argparse.ArgumentParser("Download genome sequence and label files from NCBI, for genomes from list.")

parser.add_argument('--pf-assembly-summary', required=True, help="Assembly summary file")
parser.add_argument('--pf-genome-list', required=True, help="Genome file")
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

    logger.info("Reading assembly file")
    df_assembly_summary = AssemblySummary.init_from_file(args.pf_assembly_summary)

    logger.info("Reading genome file")
    gil = GenomeInfoList.init_from_file(args.pf_genome_list)

    # only keep assemblies that match genomes from list
    wanted_gcfids = {gi.name: gi for gi in gil}
    df_assembly_summary["name"] = df_assembly_summary.apply(
        lambda r: f"{r['assembly_accession']}_{r['asm_name'].replace(' ' , '_')}", axis=1
    )




    df_assembly_summary = df_assembly_summary[df_assembly_summary["name"].isin(wanted_gcfids.keys())].copy()
    for i in df_assembly_summary.index:
        name = df_assembly_summary.at[i, "name"]
        gi = wanted_gcfids[name]        # type: GenomeInfo
        genetic_code = gi.attributes.get("genetic_code")
        df_assembly_summary.loc[i, "genetic_code"] = genetic_code

    logger.info(f"Request {len(gil)}. Found {len(df_assembly_summary)}")

    logger.info("Downloading genomes")
    download_data_from_assembly_summary(df_assembly_summary, env["pd-data"])


if __name__ == "__main__":
    main(my_env, parsed_args)
