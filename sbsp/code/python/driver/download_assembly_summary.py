# Karl Gemayel
# Georgia Institute of Technology
#
# Created: 3/17/20
import logging
import argparse
import pandas as pd
import urllib.request

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import sbsp_log  # runs init in sbsp_log and configures logger

# Custom imports
from sbsp_container.assembly_summary import AssemblySummary
from sbsp_general import Environment

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from sbsp_io.assembly_summary import read_assembly_summary_into_dataframe

parser = argparse.ArgumentParser("Download archaea and bacteria NCBI assembly summary files, and save them"
                                 " in a single tab-delimited tabular file.")

parser.add_argument("--database", required=True, choices=["refseq", "genbank"], help="Which database to use.")
parser.add_argument("--pf-output", required=True, help="Path to output file.")

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

    # link to assembly summary files
    lf_arc = f"https://ftp.ncbi.nlm.nih.gov/genomes/{args.database}/archaea/assembly_summary.txt"
    lf_bac = f"https://ftp.ncbi.nlm.nih.gov/genomes/{args.database}/bacteria/assembly_summary.txt"

    pf_out_arc = args.pf_output + "_arc"
    pf_out_bac = args.pf_output + "_bac"

    # download archaea
    logger.info(f"Downloading file: {lf_arc}")
    urllib.request.urlretrieve(lf_arc, pf_out_arc)

    # download bacteria
    logger.info(f"Downloading file: {lf_bac}")
    urllib.request.urlretrieve(lf_bac, pf_out_bac)

    logger.info("Download complete. Merging assemblies into single file")
    as_arc = AssemblySummary.init_from_file(pf_out_arc)
    as_bac = AssemblySummary.init_from_file(pf_out_bac)

    as_merged = pd.concat([as_arc, as_bac], sort=False)  # type: AssemblySummary
    as_merged.write(args.pf_output)




if __name__ == "__main__":
    main(my_env, parsed_args)
