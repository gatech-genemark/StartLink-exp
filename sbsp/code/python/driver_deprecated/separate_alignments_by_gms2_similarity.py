# Karl Gemayel
# Georgia Institute of Technology
#
# Created: 1/15/20
import logging
import argparse
import os
import shutil

import numpy
import pandas as pd
from typing import *

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import sbsp_log  # runs init in sbsp_log and configures logger

# Custom imports
from sbsp_container.genome_list import GenomeInfoList, GenomeInfo
from sbsp_general import Environment

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from sbsp_general.general import get_value
from sbsp_general.labels import Labels, Label, Coordinates
from sbsp_general.labels_comparison_detailed import LabelsComparisonDetailed
from sbsp_io.general import read_rows_to_list, mkdir_p
from sbsp_io.labels import read_labels_from_file
from sbsp_viz.general import FigureOptions
from sbsp_viz.labels_venn import venn_diagram_5prime

parser = argparse.ArgumentParser("Description of driver.")

parser.add_argument('--pf-genome-list', required=True, help="List containing genome information")
parser.add_argument('--pf-gcfid-to-pd-sbsp', required=True, help="CSV file containing GCFID to SBSP run directory")
parser.add_argument('--pf-output-summary', required=True, help="Output summary file")

parser.add_argument('--prodigal', default=False, action="store_true")

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


def create_3prime_key_from_fields(accession, left, right, strand):
    # type: (str, int, int, str) -> str
    if strand == "+":
        return "{};{};{};{}".format(
            accession, "", right, strand
        )
    else:
        return "{};{};{};{}".format(
            accession, left, "", strand
        )


def map_key_to_labels(labels):
    # type: (Labels) -> Dict[str, Label]

    result = dict()
    for l in labels:
        key = create_3prime_key_from_fields(
            accession=l.seqname(),
            left=l.left()+1,
            right=l.right()+1,
            strand=l.strand()
        )

        result[key] = l

    return result


def labels_match_5prime_3prime(a, b):
    # type: (Label, Label) -> bool
    return a.left() == b.left() and a.right() == b.right() and a.strand() == b.strand() and a.seqname() == b.seqname()


def copy_files_with_new_indexing(list_pf_source, pd_destination, fn_template="{}.txt"):
    # type: (List[str], str, str) -> None

    file_number = 0
    for pf_source in list_pf_source:

        pf_destination = os.path.join(pd_destination, fn_template.format(file_number))
        shutil.copyfile(pf_source, pf_destination)

        file_number += 1




def compare_gms2_sbsp_ncbi_for_genome_list(env, gil, gcfid_to_pd_sbsp, pf_output_summary, **kwargs):
    # type: (Environment, GenomeInfoList, Dict[str, str], str, Dict[str, Any]) -> None

    prodigal = get_value(kwargs, "prodigal", None)
    list_summary = list()
    list_pf_gms2_sbsp_not_ncbi = list()
    list_pf_gms2_sbsp_ncbi = list()

    for gi in gil:
        logger.info("{}".format(gi.name))
        pd_genome = os.path.join(env["pd-data"], gi.name)
        pf_gms2 = os.path.join(pd_genome, "runs", "gms2", "gms2.gff")
        pf_ncbi = os.path.join(pd_genome, "ncbi.gff")
        pf_sbsp_details = os.path.join(gcfid_to_pd_sbsp[gi.name], "output.csv")

        labels_gms2 = read_labels_from_file(pf_gms2, name="GMS2")
        labels_ncbi = read_labels_from_file(pf_ncbi, name="NCBI")

        key_3prime_to_label_gms2 = map_key_to_labels(labels_gms2)
        key_3prime_to_label_ncbi = map_key_to_labels(labels_ncbi)

        df_sbsp = pd.read_csv(pf_sbsp_details, header=0)

        for index, row in df_sbsp.groupby("q-key", as_index=False).agg("first").iterrows():

            q_key_3prime = create_3prime_key_from_fields(
                accession=row["q-accession"], left=row["q-left-sbsp"], right=row["q-right-sbsp"],
                strand=row["q-strand-sbsp"]
            )


            # make sure key is in both
            if q_key_3prime in key_3prime_to_label_gms2 and q_key_3prime in key_3prime_to_label_ncbi:

                # make sure SBSP 5' matches GMS2
                label_sbsp = Label(
                    Coordinates(row["q-left-sbsp"]-1, row["q-right-sbsp"]-1, row["q-strand-sbsp"]),
                    seqname=row["q-accession"]
                )

                label_gms2 = key_3prime_to_label_gms2[q_key_3prime]

                if labels_match_5prime_3prime(label_sbsp, label_gms2):

                    label_ncbi = key_3prime_to_label_ncbi[q_key_3prime]
                    if labels_match_5prime_3prime(label_sbsp, label_ncbi):
                        list_pf_gms2_sbsp_ncbi.append(row["pf-msa-output"])
                    else:
                        list_pf_gms2_sbsp_not_ncbi.append(row["pf-msa-output"])

    pd_gms2_sbsp_ncbi = os.path.join(env["pd-work"], "sbsp_gms2_ncbi")
    pd_gms2_sbsp_not_ncbi = os.path.join(env["pd-work"], "sbsp_gms2_not_ncbi")

    mkdir_p(pd_gms2_sbsp_ncbi)
    mkdir_p(pd_gms2_sbsp_not_ncbi)

    # copy files
    copy_files_with_new_indexing(list_pf_gms2_sbsp_ncbi, pd_gms2_sbsp_ncbi)
    copy_files_with_new_indexing(list_pf_gms2_sbsp_not_ncbi, pd_gms2_sbsp_not_ncbi)


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None

    gil = GenomeInfoList.init_from_file(args.pf_genome_list)
    df = pd.read_csv(args.pf_gcfid_to_pd_sbsp)

    gcfid_to_pd_sbsp = {
        x["gcfid"]: x["pd-sbsp"] for _, x in df.iterrows()
    }

    compare_gms2_sbsp_ncbi_for_genome_list(
        env,
        gil, gcfid_to_pd_sbsp, args.pf_output_summary,
        prodigal=args.prodigal
    )


if __name__ == "__main__":
    main(my_env, parsed_args)
