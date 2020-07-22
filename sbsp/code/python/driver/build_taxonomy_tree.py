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
from sbsp_container.taxonomy_tree import TaxonomyTree
from sbsp_general import Environment

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #

parser = argparse.ArgumentParser("Build a taxonomy tree object from NCBI Taxonomy Dump files.")

parser.add_argument('--pf-nodes-dmp', required=True, help="NCBI taxonomy nodes dump file")
parser.add_argument('--pf-names-dmp', required=True, help="NCBI taxonomy names dump file")
parser.add_argument('--pf-tree', required=True, help="Path to output tree file")


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
    tax_tree = TaxonomyTree.init_from_file(args.pf_nodes_dmp, args.pf_names_dmp)
    tax_tree.save(args.pf_tree)


if __name__ == "__main__":
    main(my_env, parsed_args)
