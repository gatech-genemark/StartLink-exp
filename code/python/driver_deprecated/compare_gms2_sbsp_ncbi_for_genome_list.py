# Karl Gemayel
# Georgia Institute of Technology
#
# Created: 1/15/20
import collections
import logging
import argparse
import os

import numpy
import pandas as pd
from typing import *

# noinspection All
from Bio.Seq import Seq

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
from sbsp_general.labels import Labels
from sbsp_general.labels_comparison_detailed import LabelsComparisonDetailed
from sbsp_io.general import read_rows_to_list
from sbsp_io.labels import read_labels_from_file
from sbsp_io.sequences import read_fasta_into_hash
from sbsp_viz.general import FigureOptions
from sbsp_viz.labels_venn import venn_diagram_5prime
from stats_start_candidates import count_candidates_on_positive_strand, count_candidates_on_negative_strand

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


def count_candidates_per_gene_and_update_labels(sequences, labels, **kwargs):
    # type: (Dict[str, Seq], Labels, Dict[str, Any]) -> None

    list_stats = list()

    for label in labels:
        try:
            sequence = sequences[label.seqname()]
        except KeyError:
            logger.warning("Could not get sequence for label witn sequence name: {}".format(label.seqname()))
            continue

        if label.strand() == "+":
            stats = count_candidates_on_positive_strand(sequence, label, **kwargs)
        else:
            stats = count_candidates_on_negative_strand(sequence, label, **kwargs)

        if stats is not None:
            total_candidates = stats["upstream"] + stats["downstream"]
            label.set_attribute_value("candidates", total_candidates)


def add_number_of_start_candidates_to_labels(sequences, labels):
    # type: (Dict[str, Seq], Labels) -> None
    count_candidates_per_gene_and_update_labels(sequences, labels)




def compare_gms2_sbsp_ncbi(env, pf_gms2, pf_sbsp, pf_ncbi, **kwargs):
    # type: (Environment, str, str, str, Dict[str, Any]) -> Dict[str, Any]

    venn_title = get_value(kwargs, "venn_title", None)
    pf_venn = get_value(kwargs, "pf_venn", os.path.join(env["pd-work"], "venn.pdf"))
    pf_prodigal = get_value(kwargs, "pf_prodigal", None)

    start_candidate_analysis = get_value(kwargs, "start_candidate_analysis", False)
    gcfid = get_value(kwargs, "gcfid", None)

    predicted_at_step = get_value(kwargs, "predicted_at_step", None)

    labels_gms2 = read_labels_from_file(pf_gms2, name="GMS2")
    labels_sbsp = read_labels_from_file(pf_sbsp, name="SBSP")
    labels_ncbi = read_labels_from_file(pf_ncbi, name="NCBI")

    if start_candidate_analysis:
        # add number of start candidates per gene
        sequences = read_fasta_into_hash(os.path.join(env["pd-data"], gcfid, "sequence.fasta"))
        add_number_of_start_candidates_to_labels(sequences, labels_gms2)

    if predicted_at_step is not None:
        labels_sbsp = Labels(
            [l for l in labels_sbsp if l.get_attribute_value("predicted-at-step") == predicted_at_step],
            name="SBSP"
        )


    lcd = LabelsComparisonDetailed(labels_gms2, labels_sbsp,
                                   name_a="gms2",
                                   name_b="sbsp")

    labels_gms2_sbsp_3p_5p = lcd.intersection("a")

    lcd_2 = LabelsComparisonDetailed(labels_gms2_sbsp_3p_5p, labels_ncbi,
                                     name_a="gms2_sbsp",
                                     name_b="ncbi")

    labels_gms2_sbsp_ncbi_3p_5p = lcd_2.intersection("a")

    # venn_diagram_5prime(labels_gms2, labels_sbsp, labels_ncbi, FigureOptions(
    #     title=venn_title,
    #     save_fig=pf_venn
    # ))

    labels_prodigal = None
    prodigal_info = dict()
    if pf_prodigal is not None:
        labels_prodigal = read_labels_from_file(pf_prodigal, name="Prodigal")
        lcd_prodigal = LabelsComparisonDetailed(labels_gms2_sbsp_3p_5p, labels_prodigal,
                                         name_a="gms2_sbsp",
                                         name_b="prodigal")

        labels_gms2_sbsp_prodigal_3p_5p = lcd_prodigal.intersection("a")

        # Goal: check (GMS2=SBSP) != (Prodigal=NCBI)

        # step1: Prodigal=NCBI
        labels_ncbi_prodigal_3p_5p = LabelsComparisonDetailed(labels_ncbi, labels_prodigal,
                                                              name_a="ncbi", name_b="prodigal").match_3p_5p("a")

        # Get same genes in (GMS2=SBSP) and (Prodigal=NCBI)
        lcd_full = LabelsComparisonDetailed(labels_gms2_sbsp_3p_5p, labels_ncbi_prodigal_3p_5p,
                                            name_a="gms2_sbsp", name_b="ncbi_prodigal")

        labels_match_3p = lcd_full.match_3p("a")
        labels_match_3p_5p = lcd_full.match_3p_5p("a")

        prodigal_info = {
            "(GMS2=SBSP)!=Prodigal": len(labels_gms2_sbsp_3p_5p) - len(labels_gms2_sbsp_prodigal_3p_5p),
            "(GMS2=SBSP)!=(NCBI=Prodigal)": len(labels_match_3p) - len(labels_match_3p_5p),
        }

    return {
        "GMS2": len(labels_gms2),
        "SBSP": len(labels_sbsp),
        "NCBI": len(labels_ncbi),
        "GMS2=SBSP": len(labels_gms2_sbsp_3p_5p),
        "GMS2=SBSP=NCBI": len(labels_gms2_sbsp_ncbi_3p_5p),
        **prodigal_info
    }

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

    if df.at[index, "{}-strand".format(source)] == "+":
        d = df.at[index, "{}-left".format(source)] - df.at[index, "{}-upstream_right".format(source)]
    else:
        d = df.at[index, "{}-upstream_left".format(source)] - df.at[index, "{}-right".format(source)]

    return d

def most_frequent(items):
    # type: (Iterable[Any]) -> Any
    occurence_count = collections.Counter(items)
    return occurence_count.most_common(1)[0][0]


def compute_consistency(distances, pivot, flexibility=0):
    # type: (List[int], int, int) -> float

    if flexibility == 0:
        return len([1 for d in distances if d == pivot]) / float(len(distances))
    else:
        numerator = 0
        for d in distances:
            if d is not None and  d <= pivot + flexibility and d >= pivot - flexibility:
                numerator += 1

        return float(numerator) / len(distances)
    pass


def get_upstream_info(pf_sbsp_details, **kwargs):
    # type: (str, Dict[str, Any]) -> Dict[str, Any]

    df = pd.read_csv(pf_sbsp_details, header=0)
    headers = ["{}-{}".format(x, y) for x in {"q", "t"} for y in {"upstream_left", "upstream_right", "upstream_strand", "left", "right", "strand"}]
    headers += ["q-key"]
    df = df[headers]

    accumulator = 0
    denominator = 0
    list_consistency = list()

    t_to_r_to_rpc = {
        t: {r: list() for r in {0,3}} for t in {"3", "0", "-3", "R"}
    }       # type: Dict[str, Dict[int, List[float]]]

    num_queries_considered = 0
    num_queries_assigned_to_overlap_group = 0


    # for each query
    for q_key, df_group in df.groupby("q-key", as_index=False):     # type: pd.Index, pd.DataFrame

        # skip in small numbers
        if len(df_group) < 10:
            continue

        num_queries_considered += 1

        total_genes = len(df_group) + 1

        number_with_overlap = 0
        number_close = 0
        distances = list()

        d = distance_to_upstream(df_group, df_group.index[0], "q")
        if d is not None:
            distances.append(d)

        if d is not None and d <= 0:
            number_with_overlap += 1
        if d is not None and d <= 3:
            number_close += 1

        for index in df_group.index:
            d = distance_to_upstream(df_group, index, "t")
            if d is not None:
                distances.append(d)

            if d is not None and d <= 0:
                number_with_overlap += 1

            if d is not None and d <= 3:
                number_close += 1

        # if at least one as overlap
        if number_close > 0.2*total_genes:
            num_queries_assigned_to_overlap_group += 1

            accumulator += number_close / float(total_genes)
            list_consistency.append(number_close / float(total_genes))
            denominator += 1

            # choose overlap group
            most_common_distance = most_frequent(distances)

            c_f0 = compute_consistency(distances, most_common_distance, flexibility=0)
            c_f3 = compute_consistency(distances, most_common_distance, flexibility=3)


            if most_common_distance in {3, 0, -3}:
                t_key = str(most_common_distance)
                t_to_r_to_rpc[t_key][0].append(c_f0)
                t_to_r_to_rpc[t_key][3].append(c_f3)
            elif most_common_distance < -3:
                t_key = "R"
                t_to_r_to_rpc[t_key][0].append(c_f0)
                t_to_r_to_rpc[t_key][3].append(c_f3)

    consistency = 0 if denominator == 0 else accumulator / float(denominator)
    pog = 0 if num_queries_considered == 0 else 100 * num_queries_assigned_to_overlap_group / float(num_queries_considered)
    return {
        "Overlap Consistency": consistency,
        "Overlap Consistency List": list_consistency,
        "Percent of Overlap Groups": pog,
        **{
            "D{}, RPC({},{})".format(t, t, r): t_to_r_to_rpc[t][r]
            for t in t_to_r_to_rpc.keys()
            for r in t_to_r_to_rpc[t].keys()
        }
    }


def compare_gms2_sbsp_ncbi_for_genome_list(env, gil, gcfid_to_pd_sbsp, pf_output_summary, **kwargs):
    # type: (Environment, GenomeInfoList, Dict[str, str], str, Dict[str, Any]) -> None

    prodigal = get_value(kwargs, "prodigal", None)
    list_summary = list()
    import timeit
    for gi in gil:
        logger.info("{}".format(gi.name))
        pd_genome = os.path.join(env["pd-data"], gi.name)
        pf_gms2 = os.path.join(pd_genome, "runs", "gms2", "gms2.gff")
        pf_sbsp = os.path.join(gcfid_to_pd_sbsp[gi.name], "accuracy", "{}.gff".format(gi.name))
        pf_ncbi = os.path.join(pd_genome, "ncbi.gff")

        pf_sequence = os.path.join(pd_genome, "sequence.fasta")

        t = timeit.default_timer()
        gc = compute_gc_from_file(pf_sequence)
        logger.info("Compute GC: {}".format((timeit.default_timer() - t)/60.0))

        pf_prodigal = None
        if prodigal:
            pf_prodigal = os.path.join(pd_genome, "runs", "prodigal", "prodigal.gff")

        name = gi.attributes["name"] if "name" in gi.attributes else gi.name
        ancestor = gi.attributes["ancestor"] if "ancestor" in gi.attributes else ""

        t = timeit.default_timer()
        out = compare_gms2_sbsp_ncbi(env, pf_gms2, pf_sbsp, pf_ncbi, pf_prodigal=pf_prodigal,
                               venn_title="{}, {}".format(name, ancestor),
                               pf_venn="venn_{}.pdf".format(gi.name))
        logger.info("Compare GMS2 SBSP NCBI: {}".format((timeit.default_timer() - t)/60.0))

        out["GCFID"] = gi.name
        out["Name"] = name
        out["Ancestor"] = ancestor
        out["GC"] = gc

        pf_sbsp_details = os.path.join(gcfid_to_pd_sbsp[gi.name], "output.csv")
        t = timeit.default_timer()
        if os.path.isfile(pf_sbsp_details):
            out.update(get_upstream_info(pf_sbsp_details))
        logger.info("Upstream Info: {}".format((timeit.default_timer() - t)/60.0))

        # if step information included, do analysis for steps
        valid_steps = ["A", "B", "C", "U"]

        t = timeit.default_timer()
        for v in valid_steps:
            out_step = compare_gms2_sbsp_ncbi(env, pf_gms2, pf_sbsp, pf_ncbi, pf_prodigal=pf_prodigal,
                               venn_title="{}, {}".format(name, ancestor),
                               pf_venn="venn_{}.pdf".format(gi.name),
                                              predicted_at_step=v, gcfid=gi.name)

            out["{}: GMS2=SBSP".format(v)] = out_step["GMS2=SBSP"]
            out["{}: GMS2=SBSP=NCBI".format(v)] = out_step["GMS2=SBSP=NCBI"]
        logger.info("Stepwise Compare GMS2 SBSP NCBI: {}".format((timeit.default_timer() - t)/60.0))

        list_summary.append(out)

    if len(list_summary) > 0:

        ordered_header = ["GCFID", "Name", "Ancestor", "GC", "GMS2", "SBSP", "GMS2=SBSP", "GMS2=SBSP=NCBI"]
        remaining_header = sorted(
            set(list_summary[0].keys()).difference(ordered_header)
        )
        header = ordered_header + remaining_header

        df = pd.DataFrame(list_summary, columns=header)
        df.to_csv(pf_output_summary, index=False)


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
