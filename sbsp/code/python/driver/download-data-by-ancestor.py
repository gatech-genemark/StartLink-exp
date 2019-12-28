# Karl Gemayel
# Georgia Institute of Technology
#
# Created:

import os
import sys
import shutil
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
from sbsp_general.general import run_shell_cmd, get_value
from sbsp_io.assembly_summary import get_rows_by_key
from sbsp_io.general import mkdir_p

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
    pd_runs = os.path.join(pd_gcfid, "runs")

    try:
        mkdir_p(pd_gcfid)
        mkdir_p(pd_runs)

        ftplink = gcfid_info["ftp_path"]
        fn_sequence = "{}_genomic.fna".format(gcfid)
        fn_labels = "{}_genomic.gff".format(gcfid)

        pf_ftp_sequence = os.path.join(ftplink, "{}.gz".format(fn_sequence))
        pf_ftp_labels = os.path.join(ftplink, "{}.gz".format(fn_labels))

        for not_allowed in {"#", "(", ")", ","}:
            if not_allowed in pf_ftp_sequence or not_allowed in pf_ftp_labels:
                raise ValueError("Invalid character in path")

        for not_allowed in {"#", "(", ")", "/", ":", ","}:
            if not_allowed in fn_sequence or not_allowed in fn_labels:
                raise ValueError("Invalid character in path")

        pf_local_sequence = os.path.join(pd_gcfid, "sequence.fasta")
        pf_local_labels = os.path.join(pd_gcfid, "ncbi.gff")

        # don't re-download. TODO: add option to force re-download
        if os.path.isfile(pf_local_sequence) and os.path.isfile(pf_local_labels):
            return


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
    except (IOError, OSError, ValueError):
        # cleanup failed attempt
        if os.path.exists(pd_gcfid) and os.path.isdir(pd_gcfid):
            shutil.rmtree(pd_gcfid)


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


def download_data_by_ancestor(env, ancestor_tag, tag_type, pf_taxonomy_tree, pf_assembly_summary, pd_output,
                              pf_output_list, **kwargs):
    # type: (Environment, Union[str, int], str, str, str, str, str, Dict[str, Any]) -> None

    pd_output = os.path.abspath(pd_output)

    tax_tree = TaxonomyTree.load(pf_taxonomy_tree)

    taxid_to_info_list = get_rows_by_key(pf_assembly_summary, key="taxid")

    dry_run = get_value(kwargs, "dry_run", False)

    counter = 0

    successful = 0
    failed = 0

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

                    successful += 1
                    sys.stdout.write("Download progress: {} / {} \r".format(successful, successful + failed))
                    sys.stdout.flush()
                except (IOError, OSError):
                    failed += 1
                    sys.stdout.write("Download progress: {} / {} \r".format(successful, successful + failed))
                    sys.stdout.flush()
                    pass

    if dry_run:
        print("Number of genomes: {}".format(counter))
    else:
        gil = GenomeInfoList([
            GenomeInfo("{}_{}".format(d["assembly_accession"], d["asm_name"].replace(" ", "_")), 11) for d in success_downloads
        ])

        gil.to_file(pf_output_list)


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None

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
