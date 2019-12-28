# Karl Gemayel
# Georgia Institute of Technology
#
# Created:

import os
import logging
import argparse
from datetime import datetime

from typing import *

# noinspection PyUnresolvedReferences
import pathmagic                        # add path to custom library

# Custom library imports
from sbsp_container.genome_list import GenomeInfoList, GenomeInfo
from sbsp_container.taxonomy_tree import TaxonomyTree
from sbsp_general import Environment
from sbsp_general.data_download import set_up_gcfid
from sbsp_general.general import get_value
from sbsp_io.assembly_summary import get_rows_by_key

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from sbsp_options.pbs import PBSOptions

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


def filter_list(list_info, **kwargs):
    # type: (List[Dict[str, Any]], Dict[str, Any]) -> List[Dict[str, Any]]

    if len(list_info) == 0:
        return list()

    possible_assembly_levels = {"Complete Genome", "Scaffold", "Contig"}

    valid_assembly_levels = get_value(kwargs, "valid_assembly_levels", possible_assembly_levels, default_if_none=True)
    favor_assembly_level_order = get_value(kwargs, "favor_assembly_level_order", False)
    number_per_taxid = get_value(kwargs, "number_per_taxid", None)

    list_info_filtered = list()

    def select_from_list(local_list_info, n):
        # type: (List[Dict[str, Any]], Union[None, int]) -> List[Dict[str, Any]]

        if n is None:
            return local_list_info

        if len(local_list_info) <= n:
            return local_list_info

        return local_list_info[0:n]

    list_info = sorted(list_info, reverse=True, key=lambda x: datetime.strptime(x["seq_rel_date"], "%Y/%m/%d"))

    if favor_assembly_level_order:

        for assembly_level in valid_assembly_levels:

            list_info_filtered += select_from_list(
                [x for x in list_info if x["assembly_level"] == assembly_level],
                number_per_taxid - len(list_info_filtered)
            )

            if len(list_info_filtered) == number_per_taxid:
                break
    else:
        list_info_filtered += select_from_list(list_info, number_per_taxid)

    return list_info_filtered


def get_genomes_under_ancestor_with_filters(ancestor_tag, tag_type, pf_taxonomy_tree, pf_assembly_summary, **kwargs):
    # type: (Union[str, int], str, str, str, Dict[str, Any]) -> GenomeInfoList

    tax_tree = TaxonomyTree.load(pf_taxonomy_tree)

    taxid_to_info_list = get_rows_by_key(pf_assembly_summary, key="taxid")

    list_of_genome_infos = list()
    for genome_node in tax_tree.get_possible_genomes_under_ancestor(ancestor_tag, tag_type):

        # find in assembly summary
        tax_id = genome_node["taxid"]
        if tax_id in taxid_to_info_list:
            info_list = taxid_to_info_list[tax_id]

            for gcfid_info in filter_list(info_list, **kwargs):
                list_of_genome_infos.append(gcfid_info)

    gil = GenomeInfoList([
        GenomeInfo("{}_{}".format(d["assembly_accession"], d["asm_name"]), 11) for d in list_of_genome_infos
    ])

    return gil


def download_data_by_ancestor(env, ancestor_tag, tag_type, pf_taxonomy_tree, pf_assembly_summary, pd_output,
                              pf_output_list, **kwargs):
    # type: (Environment, Union[str, int], str, str, str, str, str, Dict[str, Any]) -> None

    pd_output = os.path.abspath(pd_output)

    tax_tree = TaxonomyTree.load(pf_taxonomy_tree)

    taxid_to_info_list = get_rows_by_key(pf_assembly_summary, key="taxid")

    dry_run = get_value(kwargs, "dry_run", False)

    counter = 0

    success_downloads = list()
    for genome_node in tax_tree.get_possible_genomes_under_ancestor(ancestor_tag, tag_type):

        # find in assembly summary
        tax_id = genome_node["taxid"]
        if tax_id in taxid_to_info_list:
            info_list = taxid_to_info_list[tax_id]

            for gcfid_info in filter_list(info_list, **kwargs):
                if dry_run:
                    counter += 1
                    continue
                try:
                    set_up_gcfid(gcfid_info, pd_output)
                    success_downloads.append(gcfid_info)
                except (IOError, OSError):
                    pass

    if dry_run:
        print("Number of genomes: {}".format(counter))
    else:
        gil = GenomeInfoList([
            GenomeInfo("{}_{}".format(d["assembly_accession"], d["asm_name"]), 11) for d in success_downloads
        ])

        gil.to_file(pf_output_list)




# def filter_assembly_summary_by_ancestor( )


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None

    pbs_options = PBSOptions.init_from_dict(env, vars(args))

    # if pbs_options["use-pbs"]:
    #
    #     gil = get_genomes_under_ancestor_with_filters(
    #         args.ancestor_id, args.ancestor_id_type,
    #         args.pf_taxonomy_tree,
    #         args.pf_assembly_summary,
    #         valid_assembly_levels=args.valid_assembly_levels,
    #         favor_assembly_level_order=args.favor_assembly_level_order,
    #         number_per_taxid=args.number_per_taxid
    #     )
    #
    #     pbs = PBS(
    #         env, pbs_options,
    #         splitter=split_dataframe,
    #         merger=merge_genome_info_lists
    #     )
    #
    #     pbs.run(
    #         data={"gil": gil, },
    #         func=download_data,
    #
    #
    #
    #     )

    download_data_by_ancestor(
        env,
        args.ancestor_id,
        args.ancestor_id_type,
        args.pf_taxonomy_tree,
        args.pf_assembly_summary,
        args.pd_output,
        args.pf_output_list,
        dry_run=args.dry_run,
        valid_assembly_levels=args.valid_assembly_levels,
        favor_assembly_level_order=args.favor_assembly_level_order,
        number_per_taxid=args.number_per_taxid
    )



if __name__ == "__main__":
    main(ENV, parsed_args)
