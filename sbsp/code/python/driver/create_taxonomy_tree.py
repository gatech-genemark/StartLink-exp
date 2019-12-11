# Karl Gemayel
# Georgia Institute of Technology
#
# Created:

import os
import logging
import argparse
import pandas as pd
from typing import *

# noinspection PyUnresolvedReferences
import pathmagic                        # add path to custom library

# Custom library imports
import sbsp_general
from sbsp_general import Environment
from sbsp_container.taxonomy_tree import TaxonomyTree

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #

parser = argparse.ArgumentParser("Description of driver.")


parser.add_argument('--pf-nodes-dmp', required=True, help="NCBI taxonomy nodes dump file")
parser.add_argument('--pf-save', required=True, help="Path to tree output file")

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

    tax_tree = TaxonomyTree.init_from_file(args.pf_nodes_dmp)
    tax_tree.save(args.pf_save)


if __name__ == "__main__":
    main(ENV, parsed_args)
