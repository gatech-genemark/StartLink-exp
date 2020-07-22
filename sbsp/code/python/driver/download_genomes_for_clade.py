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
from sbsp_container.taxonomy_tree import TaxonomyTree
from sbsp_general import Environment

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from sbsp_general.data_download import download_data_by_ancestor
from sbsp_io.assembly_summary import read_assembly_summary_into_dataframe

parser = argparse.ArgumentParser("Download genome sequence and label files from NCBI, for a given clade.")

parser.add_argument('--pf-tree', required=True, help="Taxonomy tree object.")
parser.add_argument('--pf-assembly-summary', required=True, help="Assembly summary file")
parser.add_argument('--clade-id', required=True, help="Clade identifier (see clade-id-type)")
parser.add_argument('--clade-id-type', required=True, choices=["taxid", "name_txt"], help="Clade identifier type")

parser.add_argument('--pf-output-list', required=True, help="Path to output genome list")

parser.add_argument('--valid-assembly-levels', choices={"Complete Genome", "Scaffold", "Contig"}, nargs="+",
                    help="The assembly levels to choose from.")
parser.add_argument('--favor-assembly-level-order', action="store_true", default=False,
                    help="If set, Complete Genome is favored over Scaffold, favored over Contig")
parser.add_argument('--genomes-per-taxid', type=int, default=None,
                    help="The number of allowed genomes with the same taxonomy id number")

parser.add_argument('--dry-run', default=False, action="store_true", help="If set, dry run; nothing is downloaded")

parser.add_argument('--force-download', type=str, choices=["any", "annotation_changed", "no_download"], default=None)

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

    logger.info("Reading Taxonomy Tree")
    taxonomy_tree = TaxonomyTree.load(args.pf_tree)

    logger.info("Reading assembly file")
    df_assembly_summary = AssemblySummary.init_from_file(args.pf_assembly_summary)

    logger.info("Downloading genomes")
    download_data_by_ancestor(args.clade_id, args.clade_id_type, taxonomy_tree, df_assembly_summary,
                              env["pd-data"], pf_output_list=args.pf_output_list, dry_run=args.dry_run,
                              valid_assembly_levels=args.valid_assembly_levels,
                              favor_assembly_level_order=args.favor_assembly_level_order,
                              number_per_taxid=args.genomes_per_taxid,
                              force_download=args.force_download)


if __name__ == "__main__":
    main(my_env, parsed_args)
