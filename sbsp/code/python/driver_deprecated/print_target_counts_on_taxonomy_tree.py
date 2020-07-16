# Karl Gemayel
# Georgia Institute of Technology
#
# Created: April 26, 2019

import logging
import argparse
import pandas as pd
from typing import *

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import sbsp_log  # runs init in sbsp_log and configures logger

from sbsp_general import Environment
from sbsp_container.taxonomy_tree import TaxonomyTree, Attributes, AttributeUpdater
from sbsp_general.general import get_value
from sbsp_io.assembly_summary import get_rows_by_key, read_assembly_summary, read_assembly_summary_into_dataframe

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from sbsp_io.general import write_string_to_file, read_rows_to_list
from sbsp_io.objects import save_obj, load_obj

parser = argparse.ArgumentParser("Description of driver.")

parser.add_argument('--pf-taxonomy-tree', required=True, help="Pickle file containing taxonomy tree")
parser.add_argument('--pf-assembly-summary', required=True,
                    help="Assembly summary file containing taxa id and GCFID info")
parser.add_argument('--pf-sbsp-output', required=True, help="File containing target information")

parser.add_argument('--pf-names-of-interest', required=False, help="Names (one per line) for nodes of interest")
parser.add_argument('--valid-assembly-levels', required=False, default=None, nargs="+",
                    choices=["Complete Genome", "Scaffold", "Contig"])

parser.add_argument('--max-depth', type=int, default=None)

parser.add_argument('--pf-output', required=True, help="Path to output file")

parser.add_argument('--pf-save-state')
parser.add_argument('--pf-load-state')

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


def get_assembly_info_per_gcfid(df_assembly_summary):
    # type: (pd.DataFrame) -> Dict[str, pd.Series]

    gcfid_to_assembly_info = dict()

    for index, row in df_assembly_summary.iterrows():

        gcf = row["assembly_accession"]
        acc = row["asm_name"].replace(" ", "_")
        gcfid = "{}_{}".format(gcf, acc)

        gcfid_to_assembly_info[gcfid] = row

    return gcfid_to_assembly_info


def count_targets_per_gcfid(pf_sbsp_output):
    # type: (str) -> Dict[str, int]
    df = pd.read_csv(pf_sbsp_output, header=0)
    return dict(df["t-genome"].value_counts())


def set_number_of_targets_per_gcfid(attributes, parent_attributes, children_attributes, **kwargs):
    # type: (Attributes, Attributes, List[Attributes], Dict[str, Any]) -> None

    gcfid_to_number_of_targets = get_value(kwargs, "gcfid_to_number_of_targets", required=True)

    attributes["number_of_targets"] = 0
    if len(children_attributes) > 0:
        attributes["number_of_targets"] = sum(child["number_of_targets"] for child in children_attributes)

    if attributes["gcfid"] in gcfid_to_number_of_targets:
        attributes["number_of_targets"] += gcfid_to_number_of_targets[attributes["gcfid"]]


def set_number_of_targets_per_taxid(attributes, parent_attributes, children_attributes, **kwargs):
    # type: (Attributes, Attributes, List[Attributes], Dict[str, Any]) -> None

    gcfid_to_number_of_targets = get_value(kwargs, "taxid_to_number_of_targets", required=True)

    attributes["number_of_targets"] = 0
    if len(children_attributes) > 0:
        attributes["number_of_targets"] = sum(child["number_of_targets"] for child in children_attributes)

    if attributes["taxid"] in gcfid_to_number_of_targets:
        attributes["number_of_targets"] += gcfid_to_number_of_targets[attributes["taxid"]]


def should_print(attributes):
    # type: (Attributes) -> bool

    if "number_of_targets" in attributes:
        if attributes["number_of_targets"] > 0:
            return True

    return False


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None

    if args.pf_load_state is None:
        gcfid_to_number_of_targets = count_targets_per_gcfid(args.pf_sbsp_output)



        df_assembly_summary = read_assembly_summary_into_dataframe(args.pf_assembly_summary)
        gcfid_to_assembly_info = get_assembly_info_per_gcfid(df_assembly_summary)


        taxid_to_number_of_targets = {
            int(gcfid_to_assembly_info[gcfid]["taxid"]): gcfid_to_number_of_targets[gcfid] for gcfid in gcfid_to_number_of_targets if gcfid in gcfid_to_assembly_info
        }

        tree = TaxonomyTree.load(args.pf_taxonomy_tree)

        tree.update_tree_attributes(
            set_number_of_targets_per_taxid,
            {"taxid_to_number_of_targets": taxid_to_number_of_targets},
            direction="bottom-up"
        )

        if args.pf_save_state is not None:
            save_obj(tree, args.pf_save_state)
    else:
        tree = load_obj(args.pf_load_state)

    tree_string = tree.to_string(check_if_should_print=should_print, attribute_name="number_of_targets", attribute_format="{:,}",tag_name=args.tag, max_depth=args.max_depth)
    write_string_to_file(tree_string, args.pf_output)


if __name__ == "__main__":
    main(my_env, parsed_args)
# GCA_002186515.1_ASM218651v1