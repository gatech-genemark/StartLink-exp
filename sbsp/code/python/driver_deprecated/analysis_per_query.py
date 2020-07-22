# Karl Gemayel
# Georgia Institute of Technology
#
# Created: 1/15/20
import logging
import argparse
import os

import pandas as pd
from typing import *

# noinspection All
from Bio.Seq import Seq

# noinspection PyUnresolvedReferences
import pathmagic
import sbsp_log  # runs init in sbsp_log and configures logger

# Custom imports
from sbsp_container.genome_list import GenomeInfoList, GenomeInfo
from sbsp_general import Environment

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from sbsp_general.general import os_join
from sbsp_general.labels import Label
from sbsp_general.shelf import add_q_key_3p_to_df, map_key_3p_to_label, map_key_3p_to_df_group, labels_match_5p_3p, \
    append_data_frame_to_csv
from sbsp_io.general import remove_p
from sbsp_io.labels import read_labels_from_file

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


def compute_gc_from_file(pf_sequence):
    # type: (str) -> float

    from sbsp_io.sequences import read_fasta_into_hash
    sequences = read_fasta_into_hash(pf_sequence)


    counts = {"A": 0, "C": 0, "G": 0, "T": 0}

    for seq in sequences.values():
        for s in seq:
            if s.upper() in {"A", "C", "G", "T"}:
                counts[s] += 1

    total = sum(counts.values())
    count_gc = counts["G"] + counts["C"]

    if total == 0:
        return 0.0

    return count_gc / float(total)



def distance_to_upstream(df, index, source):
    # type: (pd.DataFrame, pd.Index, str) -> Union[int, None]
    """
    0 means overlap by 1 nt. Positive numbers mean no overlap. Negatives mean overlap
    :param series:
    :param source:
    :return:
    """

    # if no upstream gene
    if df.at[index, "{}-upstream_left".format(source)] == -1:
        return None

    if df.at[index, "{}-strand".format(source)] != df.at[index, "{}-upstream_strand".format(source)]:
        return None

    if df.at[index, "{}-strand".format(source)] == "+":
        d = df.at[index, "{}-left".format(source)] - df.at[index, "{}-upstream_right".format(source)]
    else:
        d = df.at[index, "{}-upstream_left".format(source)] - df.at[index, "{}-right".format(source)]

    return d


def list_of_upstream_distances_from_df(df):
    # type: (pd.DataFrame) -> List[int]
    d = distance_to_upstream(df, df.index[0], "q")
    distances = list()
    if d is not None:
        distances.append(int(d))

    for index in df.index:
        d = distance_to_upstream(df, index, "t")
        if d is not None:
            distances.append(int(d))
    return distances


def analyze_query(key, key_to_label_sbsp, key_to_label_gms2, key_to_label_ncbi, key_to_label_prodigal,
                  key_to_df_sbsp_details):
    # type: (str, Dict[str, Label], Dict[str, Label], Dict[str, Label], Dict[str, Label], Dict[str, pd.DataFrame]) -> Dict[str, Any]

    r = dict()

    l_sbsp = key_to_label_sbsp[key] if key in key_to_label_sbsp else None
    l_gms2 = key_to_label_gms2[key] if key in key_to_label_gms2 else None
    l_ncbi = key_to_label_ncbi[key] if key in key_to_label_ncbi else None
    l_prodigal = key_to_label_prodigal[key] if key in key_to_label_prodigal else None

    r["SBSP"] = l_sbsp is not None
    r["GMS2"] = l_gms2 is not None
    r["NCBI"] = l_ncbi is not None
    r["Prodigal"] = l_prodigal is not None

    r["GMS2=SBSP"] = r["SBSP"] and r["GMS2"] and labels_match_5p_3p(l_sbsp, l_gms2)
    r["GMS2=NCBI"] = r["GMS2"] and r["NCBI"] and labels_match_5p_3p(l_gms2, l_ncbi)
    r["Prodigal=NCBI"] = r["Prodigal"] and r["NCBI"] and labels_match_5p_3p(l_prodigal, l_ncbi)
    r["SBSP=NCBI"] = r["SBSP"] and r["NCBI"] and labels_match_5p_3p(l_sbsp, l_ncbi)


    r["GMS2=SBSP=NCBI"] = r["GMS2=SBSP"] and r["NCBI"] and labels_match_5p_3p(l_sbsp, l_ncbi)
    r["GMS2=SBSP=Prodigal"] = r["GMS2=SBSP"] and r["Prodigal"] and labels_match_5p_3p(l_sbsp, l_prodigal)

    r["(GMS2=SBSP)!=(NCBI=Prodigal)"] = r["GMS2=SBSP"] and r["Prodigal"] and r["NCBI"] and not r["GMS2=SBSP=NCBI"] and labels_match_5p_3p(l_prodigal, l_ncbi)


    # sbsp details
    detail_keys = ["Support", "Predicted-at-step", "Kimura-to-query", "Upstream-distance"]
    for k in detail_keys:
        r[k] = None

    if r["SBSP"]:
        df = key_to_df_sbsp_details[key]
        r["Support"] = len(df)
        r["Predicted-at-step"] = df.at[df.index[0], "predicted-at-step"]
        r["Kimura-to-query"] = list(df["distance"])
        r["Upstream-distance"] = list_of_upstream_distances_from_df(df)

    return r


def analysis_per_query_for_genome(env, gi, pd_sbsp, **kwargs):
    # type: (Environment, GenomeInfo, str, Dict[str, Any]) -> pd.DataFrame

    pd_genome = os_join(env["pd-data"], gi.name)
    pf_gms2 = os_join(pd_genome, "runs", "gms2", "gms2.gff")
    pf_prodigal = os_join(pd_genome, "runs", "prodigal", "prodigal.gff")
    pf_sbsp = os_join(pd_sbsp, "accuracy", "{}.gff".format(gi.name))
    pf_ncbi = os_join(pd_genome, "ncbi.gff")
    pf_sbsp_details = os_join(pd_sbsp, "output.csv")


    # Read all input and sbsp prediction details
    common_options = {"shift": 0}
    labels_sbsp = read_labels_from_file(pf_sbsp, name="SBSP", **common_options)
    labels_gms2 = read_labels_from_file(pf_gms2, name="GMS2", **common_options)
    labels_ncbi = read_labels_from_file(pf_ncbi, name="NCBI", **common_options)
    labels_prodigal = read_labels_from_file(pf_prodigal, name="Prodigal", **common_options)
    df_sbsp_details = pd.read_csv(pf_sbsp_details)
    add_q_key_3p_to_df(df_sbsp_details, "q-key-3p")

    # get keys per label
    key_to_label_sbsp = map_key_3p_to_label(labels_sbsp)
    key_to_label_gms2 = map_key_3p_to_label(labels_gms2)
    key_to_label_ncbi = map_key_3p_to_label(labels_ncbi)
    key_to_label_prodigal = map_key_3p_to_label(labels_prodigal)
    key_to_df_sbsp_details = map_key_3p_to_df_group(df_sbsp_details)

    df_result = pd.DataFrame()


    # Sketch: Dataframe will contain one row per gene (3prime end), for all genes in
    # the union set of SBSP, GMS2, NCBI, and prodigal
    all_key_3p = set(key_to_label_sbsp.keys()).union(
        set(key_to_label_gms2.keys()),
        set(key_to_label_ncbi.keys()),
        set(key_to_label_prodigal)
    )

    list_analysis = list()
    for key in all_key_3p:

        curr_analysis = analyze_query(
            key, key_to_label_sbsp, key_to_label_gms2, key_to_label_ncbi, key_to_label_prodigal, key_to_df_sbsp_details
        )
        list_analysis.append(curr_analysis)

    if len(list_analysis) == 0:
        return pd.DataFrame()

    return pd.DataFrame(list_analysis)


def analysis_per_query(env, gil, gcfid_to_pd_sbsp, pf_output_summary, **kwargs):
    # type: (Environment, GenomeInfoList, Dict[str, str], str, Dict[str, Any]) -> None

    if os.path.isfile(pf_output_summary):
        remove_p(pf_output_summary)

    counter = 0
    for gi in gil:
        logger.info("{} / {}: {}".format(counter, len(gil), gi.name))
        pd_genome = os.path.join(env["pd-data"], gi.name)
        pf_sequence = os.path.join(pd_genome, "sequence.fasta")
        gc = compute_gc_from_file(pf_sequence)

        df = analysis_per_query_for_genome(env, gi, gcfid_to_pd_sbsp[gi.name])
        df["GCFID"] = gi.name
        df["Name"] = gi.attributes["name"] if "name" in gi.attributes else gi.name
        df["Genome GC"] = gc
        df["Ancestor"] = gi.attributes["ancestor"] if "ancestor" in gi.attributes else ""

        append_data_frame_to_csv(df, pf_output_summary)
        counter += 1


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None

    gil = GenomeInfoList.init_from_file(args.pf_genome_list)
    df = pd.read_csv(args.pf_gcfid_to_pd_sbsp)

    gcfid_to_pd_sbsp = {
        x["gcfid"]: x["pd-sbsp"] for _, x in df.iterrows()
    }

    analysis_per_query(
        env,
        gil, gcfid_to_pd_sbsp, args.pf_output_summary,
        prodigal=args.prodigal
    )


if __name__ == "__main__":
    main(my_env, parsed_args)
