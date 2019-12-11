# Karl Gemayel
# Georgia Institute of Technology
#
# Created:

import os
import logging
import argparse

from typing import *

# noinspection PyUnresolvedReferences
import pathmagic                        # add path to custom library

# Custom library imports
from sbsp_container.genome_list import GenomeInfoList, GenomeInfo
from sbsp_container.taxonomy_tree import TaxonomyTree
from sbsp_general import Environment
from sbsp_general.general import run_shell_cmd, get_value
from sbsp_io.assembly_summary import get_rows_by_key
from sbsp_io.general import mkdir_p

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #

parser = argparse.ArgumentParser("Description of driver.")

ancestor_grp = parser.add_mutually_exclusive_group(required=True)
ancestor_grp.add_argument('--ancestor-name', help="Name of ancestor")
ancestor_grp.add_argument('--ancestor-tax-id', type=int, help="Ancestor taxonomy id")

parser.add_argument('--pf-taxonomy-tree', required=True, help="Pickle file containing taxonomy tree")
parser.add_argument('--pf-assembly-summary', required=True, help="Assembly summary file")

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


def set_up_gcfid(gcfid_info, pd_output):
    # type: (Dict[str, Any], str) -> None

    # build name
    gcf = gcfid_info["assembly_accession"]
    acc = gcfid_info["asm_name"].replace(" ", "_")

    gcfid = "{}_{}".format(gcf, acc)

    pd_gcfid = os.path.join(pd_output, gcfid)
    mkdir_p(pd_gcfid)

    ftplink = gcfid_info["ftp_path"]
    fn_sequence = "{}_genomic.fna".format(gcfid)
    fn_labels = "{}_genomic.gff".format(gcfid)

    pf_ftp_sequence = os.path.join(ftplink, "{}.gz".format(fn_sequence))
    pf_ftp_labels = os.path.join(ftplink, "{}.gz".format(fn_labels))

    pf_local_sequence = os.path.join(pd_gcfid, "sequence.fasta")
    pf_local_labels = os.path.join(pd_gcfid, "ncbi.gff")

    run_shell_cmd(
        "pwd; cd {}; wget --quiet {}; wget --quiet {}; gunzip -f {}; gunzip -f {}".format(
            pd_gcfid,
            pf_ftp_sequence,
            pf_ftp_labels,
            "{}.gz".format(fn_sequence),
            "{}.gz".format(fn_labels)
        ),

    )

    run_shell_cmd(
        "cd {}; mv {} {}; mv {} {}".format(
            pd_gcfid,
            fn_sequence, "sequence.fasta",
            fn_labels, "ncbi.gff"
        )
    )








def download_data_by_ancestor(env, ancestor_tag, tag_type, pf_taxonomy_tree, pf_assembly_summary, pd_output,
                              pf_output_list, **kwargs):
    # type: (Environment, Union[str, int], str, str, str, str, str, Dict[str, Any]) -> None

    tax_tree = TaxonomyTree.load(pf_taxonomy_tree)

    taxid_to_info_list = get_rows_by_key(pf_assembly_summary, key="taxid")
    counter = 0
    total = 0

    dry_run = get_value(kwargs, "dry_run", False)

    counter = 0

    success_downloads = list()
    for genome_node in tax_tree.get_possible_genomes_under_ancestor(ancestor_tag, tag_type):

        # find in assembly summary
        tax_id = genome_node["taxid"]
        if tax_id in taxid_to_info_list:
            info_list = taxid_to_info_list[tax_id]

            for gcfid_info in info_list:
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



def main(env, args):
    # type: (Environment, argparse.Namespace) -> None

    pd_output = os.path.abspath(args.pd_output)

    download_data_by_ancestor(env, args.ancestor_tax_id, "taxid", args.pf_taxonomy_tree, args.pf_assembly_summary,
                              pd_output, args.pf_output_list, dry_run=args.dry_run)



if __name__ == "__main__":
    main(ENV, parsed_args)
