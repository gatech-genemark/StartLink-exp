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
from sbsp_container.taxonomy_tree import TaxonomyTree, Attributes, AttributeUpdater
from sbsp_general.general import get_value
from sbsp_io.assembly_summary import get_rows_by_key

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from sbsp_io.general import write_string_to_file, read_rows_to_list
from sbsp_pipeline.pipeline_msa import PipelineMSA

parser = argparse.ArgumentParser("Description of driver.")

parser.add_argument('--pf-taxonomy-tree', required=True, help="Pickle file containing taxonomy tree")
parser.add_argument('--pf-assembly-summary', required=True, help="Assembly summary file")
parser.add_argument('--pf-names-of-interest', required=False, help="Names (one per line) for nodes of interest")
parser.add_argument('--valid-assembly-levels', required=False, default=None, nargs="+",
                    choices=["Complete Genome", "Scaffold", "Contig"])

parser.add_argument('--pf-output', required=True, help="Path to output file")

parser.add_argument('--tag', required=False, help="If set, use tag as name in tree")

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


def count_refseq_under_node(children_attributes, curr_node_attributes, attribute_name, **kwargs):
    # type: (List[Dict[str, Any]], Dict[str, Any], str, Dict[str, Any]) -> Any

    refseq_count_per_taxid = get_value(kwargs, "refseq_count_per_taxid", required=True)
    limit_path_to = get_value(kwargs, "limit_path_to", None)

    num_refseq = 0
    for c in children_attributes:
        num_refseq += c[attribute_name]

    if curr_node_attributes["taxid"] in refseq_count_per_taxid:
        num_refseq += refseq_count_per_taxid[curr_node_attributes["taxid"]]

    if limit_path_to is not None:
        leads_to_node_of_interest = False
        for c in children_attributes:
            if c["leads_to_node_of_interest"] == True:
                leads_to_node_of_interest = True

        if curr_node_attributes["name_txt"] in limit_path_to:
            leads_to_node_of_interest = True

        curr_node_attributes["leads_to_node_of_interest"] = leads_to_node_of_interest

    return num_refseq


def print_taxonomy_tree(env, pf_taxonomy_tree, pf_assembly_summary, pf_output, **kwargs):
    pf_names_of_interest = get_value(kwargs, "pf_names_of_interest", None)
    tax_tree = TaxonomyTree.load(pf_taxonomy_tree)
    taxid_to_info_list = get_rows_by_key(pf_assembly_summary, key="taxid", **kwargs)

    limit_path_to = None
    if pf_names_of_interest:
        limit_path_to = set(read_rows_to_list(pf_names_of_interest))

    refseq_count_per_taxid = {
        taxid: len(taxid_to_info_list[taxid]) for taxid in taxid_to_info_list
    }

    refseq_count_per_taxid = {
        taxid: 1 for taxid in taxid_to_info_list
    }

    def check_if_should_print(attributes):
        # type: (Dict[str, Any]) -> bool
        if "leads_to_node_of_interest" in attributes:
            return attributes["leads_to_node_of_interest"]
        return "num_refseq" in attributes and attributes["num_refseq"] != 0

    out = tax_tree.to_string_tree_with_stats("num_refseq", count_refseq_under_node, {
        "refseq_count_per_taxid": refseq_count_per_taxid,
        "limit_path_to": limit_path_to,
    }, check_if_should_print=check_if_should_print, attribute_format="{:,}", **kwargs)

    write_string_to_file(out, pf_output)


def updater_add_gcfid(attributes, parent_attributes, children_attributes, **kwargs):
    # type: (Attributes, Union[Attributes, None], List[Attributes], Dict[str, Any]) -> None

    assembly_summary_per_taxid = get_value(kwargs, "assembly_summary_per_taxid", required=True)

    tax_id = attributes["tax_id"]
    try:
        info = assembly_summary_per_taxid[tax_id]
        gcf = info["assembly_accession"]
        acc = info["asm_name"].replace(" ", "_")
        gcfid = "{}_{}".format(gcf, acc)
    except KeyError:
        gcfid = "noname"

    attributes["gcfid"] = gcfid


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None
    print_taxonomy_tree(env, args.pf_taxonomy_tree, args.pf_assembly_summary, args.pf_output,
                        valid_assembly_levels=args.valid_assembly_levels,
                        tag_name=args.tag,
                        pf_names_of_interest=args.pf_names_of_interest)


if __name__ == "__main__":
    main(my_env, parsed_args)
