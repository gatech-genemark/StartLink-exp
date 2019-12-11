# Karl Gemayel
# Georgia Institute of Technology
#
# Created: April 26, 2019

import logging
import argparse

from typing import *

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import sbsp_log  # runs init in sbsp_log and configures logger

from sbsp_general import Environment
from sbsp_container.taxonomy_tree import TaxonomyTree
from sbsp_general.general import get_value
from sbsp_io.assembly_summary import get_rows_by_key


# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from sbsp_io.general import write_string_to_file
from sbsp_pipeline.pipeline_msa import PipelineMSA

parser = argparse.ArgumentParser("Description of driver.")

parser.add_argument('--pf-taxonomy-tree', required=True, help="Pickle file containing taxonomy tree")
parser.add_argument('--pf-assembly-summary', required=True, help="Assembly summary file")

parser.add_argument('--pf-output', required=True, help="Path to output file")

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


def count_refseq_under_node(curr_node_attributes, children_attributes, attribute_name, **kwargs):
    # type: (Dict[str, Any], List[Dict[str, Any]], str, Dict[str, Any]) -> Any

    refseq_count_per_taxid = get_value(kwargs, "refseq_count_per_taxid", required=True)

    num_refseq = 0
    for c in children_attributes:
        num_refseq += c[attribute_name]

    if curr_node_attributes["taxid"] in refseq_count_per_taxid:
        num_refseq += refseq_count_per_taxid[curr_node_attributes["taxid"]]

    return num_refseq


def print_taxonomy_tree(env, pf_taxonomy_tree, pf_assembly_summary, pf_output):
    tax_tree = TaxonomyTree.load(pf_taxonomy_tree)
    taxid_to_info_list = get_rows_by_key(pf_assembly_summary, key="taxid")

    refseq_count_per_taxid = {
        taxid: len(taxid_to_info_list[taxid]) for taxid in taxid_to_info_list
    }

    out = tax_tree.to_string_tree_with_stats("num_refseq", count_refseq_under_node, {
        "refseq_count_per_taxid": refseq_count_per_taxid
    })

    write_string_to_file(out, pf_output)


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None
    print_taxonomy_tree(env, args.pf_taxonomy_tree, args.pf_assembly_summary, args.pf_output)


if __name__ == "__main__":
    main(my_env, parsed_args)
