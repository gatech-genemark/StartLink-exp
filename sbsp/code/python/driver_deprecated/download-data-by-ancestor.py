# Karl Gemayel
# Georgia Institute of Technology
#
# Created:

import logging
import argparse

from typing import *

# noinspection PyUnresolvedReferences
import pathmagic                        # add path to custom library

# Custom library imports
from sbsp_container.taxonomy_tree import TaxonomyTree
from sbsp_general import Environment
from sbsp_general.data_download import download_data_by_ancestor
from sbsp_io.assembly_summary import read_assembly_summary_into_dataframe

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #

parser = argparse.ArgumentParser("Description of driver.")

parser.add_argument("--ancestor-id", required=True, help="ID of ancestor")
parser.add_argument("--ancestor-id-type", required=True, choices=["taxid", "name_txt"], help="Type of ID for ancestor. See choices.")

parser.add_argument('--pf-taxonomy-tree', required=True, help="Pickle file containing taxonomy tree")
parser.add_argument('--pf-assembly-summary', required=True, help="Assembly summary file")

parser.add_argument('--valid-assembly-levels', choices={"Complete Genome", "Scaffold", "Contig"}, nargs="+")
parser.add_argument('--favor-assembly-level-order', action="store_true", default=False)
parser.add_argument('--number-per-taxid', type=int, default=None)

parser.add_argument('--pd-output', required=True, help="Path to output directory where data will be downloaded")
parser.add_argument('--pf-output-list', required=True, help="Path to list file containing genome names")

parser.add_argument('--dry-run', default=False, action="store_true")

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
ENV = Environment(pd_data=parsed_args.pd_data, pd_work=parsed_args.pd_work,
                  pd_results=parsed_args.pd_results)

# Setup logger
logging.basicConfig(level=parsed_args.loglevel)
logger = logging.getLogger("logger")                    # type: logging.Logger


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None

    taxonomy_tree = TaxonomyTree.load(args.pf_taxonomy_tree)
    df_assembly_summary = read_assembly_summary_into_dataframe(args.pf_assembly_summary)

    download_data_by_ancestor(args.ancestor_id, args.ancestor_id_type, taxonomy_tree, df_assembly_summary,
                              args.pd_output, pf_output_list=args.pf_output_list, dry_run=args.dry_run,
                              valid_assembly_levels=args.valid_assembly_levels,
                              favor_assembly_level_order=args.favor_assembly_level_order,
                              number_per_taxid=args.number_per_taxid,
                              force_download=args.force_download)



if __name__ == "__main__":
    main(ENV, parsed_args)
