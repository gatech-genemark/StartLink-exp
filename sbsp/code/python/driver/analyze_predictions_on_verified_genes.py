# Karl Gemayel
# Georgia Institute of Technology
#
# Created: 3/17/20
import copy
import logging
import argparse
import pandas as pd
from typing import *

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import sbsp_log  # runs init in sbsp_log and configures logger

# Custom imports
from sbsp_container.genome_list import GenomeInfoList, GenomeInfo
from sbsp_general import Environment
from sbsp_viz.colormap import ColorMap as CM

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from sbsp_general.general import get_value, os_join
from sbsp_general.labels import Labels
from sbsp_general.labels_comparison_detailed import LabelsComparisonDetailed
from sbsp_general.shelf import map_key_3p_to_df_group, map_key_3p_to_label, next_name, add_q_key_3p_to_df
from sbsp_io.labels import read_labels_from_file
from sbsp_viz import sns
from sbsp_viz.general import FigureOptions

parser = argparse.ArgumentParser("Description of driver.")

parser.add_argument('--pf-genome-list', required=True, help="List containing genome information")
parser.add_argument('--pf-gcfid-to-pd-sbsp', required=True, help="CSV file containing GCFID to SBSP run directory")
parser.add_argument('--fn-prefix', required=False)

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


def get_stats_a_from_b_3p(labels_a, labels_b, output):
    # type: (Labels, Labels, Dict[str, Any]) -> None

    lcd = LabelsComparisonDetailed(labels_a, labels_b)

    a_name = labels_a.name if labels_a.name is not None else "a"
    b_name = labels_b.name if labels_b.name is not None else "b"

    output.update({
        a_name: len(labels_a),
        b_name: len(labels_b),
        "Number 3p match: {} from {}".format(a_name, b_name): len(lcd.match_3p("a")),
        "Percentage 3p match: {} from {}".format(a_name, b_name): round(
            100 * len(lcd.match_3p("a")) / float(len(labels_b)), 2)
    })


def get_stats_a_from_b_3p_by_upstream(labels_a, labels_b, output):
    # type: (Labels, Labels, Dict[str, Any]) -> None

    a_name = labels_a.name if labels_a.name is not None else "a"
    b_name = labels_b.name if labels_b.name is not None else "b"


def get_stats_sn_sp(labels_a, labels_b, output):
    # type: (Labels, Labels, Dict[str, Any]) -> None

    lcd = LabelsComparisonDetailed(labels_a, labels_b)

    a_name = labels_a.name if labels_a.name is not None else "a"
    b_name = labels_b.name if labels_b.name is not None else "b"

    a_total = len(labels_a)
    b_total = len(labels_b)

    match_3p = len(lcd.match_3p("a"))
    match_3p_5p = len(lcd.match_3p_5p("a"))

    sp = 0 if a_total == 0 else round(100 * match_3p / float(a_total), 2)
    sn = 0 if match_3p == 0 else round(100 * match_3p_5p / float(match_3p), 2)

    output.update({
        a_name: a_total,
        b_name: b_total,
        "Number 3p match: {} from {}".format(a_name, b_name): match_3p,
        "Percentage 3p match: {} from {}".format(a_name, b_name): sp,
        "Number 5p-3p match: {} from {}".format(a_name, b_name): match_3p_5p,
        "Percentage 5p-3p match: {} from {}".format(a_name, b_name): sn,
    })


def add_support_to_labels(labels, df):
    # type: (Labels, pd.DataFrame) -> None

    key_to_df = map_key_3p_to_df_group(df)
    key_to_label = map_key_3p_to_label(labels)

    for key, label in key_to_label.items():
        support = 0
        if key in key_to_df:
            df_curr = key_to_df[key]
            support = len(df_curr)

        label.set_attribute_value("support", support)


def get_stats_sn_sp_by_support(labels_ref, labels_pred, stats, tag):
    # type: (Labels, Labels, Dict[str, Any], str) -> None

    min_support = min(lab.get_attribute_value("support") for lab in labels_pred)
    max_support = max(lab.get_attribute_value("support") for lab in labels_pred)

    sorted_labels_pred = [l for l in labels_pred.sort_by_attribute("support", in_place=False)]

    curr_index = 0
    num_labels = len(labels_pred)

    results = list()

    # for each support, compute sn/sp for all labels with support >= s
    for s in range(min_support, max_support+1):

        while curr_index < num_labels and sorted_labels_pred[curr_index].get_attribute_value("support") < s:
            curr_index += 1

        if curr_index == num_labels:
            break

        # compute
        curr_labels = Labels(sorted_labels_pred[curr_index:], name=labels_pred.name)

        results_s = dict()

        get_stats_sn_sp(labels_ref, curr_labels, results_s)
        results_s["Min Support"] = s

        results.append(results_s)

    stats["by_support_{}".format(tag)] = results


def get_stats_sn_sp_by_step_group(labels_ref, labels_pred, stats, tag):
    # type: (Labels, Labels, Dict[str, Any], str) -> None


    step_groups = [("A", "B", "C", "U"), ("A", "B", "U"), ("A", "B"), ("A")]

    step_to_list_labels = {s: list() for s in {"A", "B", "U", "C"}}
    for l in labels_pred:
        step_to_list_labels[l.get_attribute_value("predicted-at-step")].append(l)
    curr_index = 0
    num_labels = len(labels_pred)

    results = list()

    # for each support, compute sn/sp for all labels with support >= s
    for sg in step_groups:

        # get all labels for group
        curr_labels = Labels(
            [l for s in sg for l in step_to_list_labels[s]],
            name=labels_pred.name
        )

        results_s = dict()

        get_stats_sn_sp(labels_ref, curr_labels, results_s)
        results_s["Step Group"] = str(sg)

        results.append(results_s)

    stats["by_step_group_{}".format(tag)] = results


def analyze_predictions_on_verified_genes(env, gi, pd_sbsp, **kwargs):
    # type: (Environment, GenomeInfo, str, Dict[str, Any]) -> Dict[str, Any]
    pd_gcfid = os_join(env["pd-data"], gi.name)

    pf_sbsp = os_join(pd_sbsp, "accuracy", "{}.gff".format(gi.name))
    pf_gms2 = os_join(pd_gcfid, "runs", "gms2", "gms2.gff")
    pf_verified = os_join(pd_gcfid, "verified.gff")
    pf_ncbi = os_join(pd_gcfid, "ncbi.gff")
    pf_sbsp_details = os_join(pd_sbsp, "output.csv")

    kwargs_labels = {"ignore_frameshifted": True, "ignore_partial": True, "shift": 0}

    labels_sbsp = read_labels_from_file(pf_sbsp, name="SBSP", **kwargs_labels)
    labels_verified = read_labels_from_file(pf_verified, name="Verified", **kwargs_labels)
    labels_gms2 = read_labels_from_file(pf_gms2, name="GMS2", **kwargs_labels)
    labels_ncbi = read_labels_from_file(pf_ncbi, name="NCBI", **kwargs_labels)
    df_sbsp_details = pd.read_csv(pf_sbsp_details)
    add_q_key_3p_to_df(df_sbsp_details, "q-key-3p")

    add_support_to_labels(labels_sbsp, df_sbsp_details)

    labels_sbsp_eq_gms2 = LabelsComparisonDetailed(labels_sbsp, labels_gms2).match_3p_5p("a")
    labels_sbsp_eq_gms2.name = "GMS2=SBSP"

    stats = dict()

    # Stats: 3prime match
    get_stats_a_from_b_3p(labels_verified, labels_ncbi, stats)
    get_stats_a_from_b_3p(labels_verified, labels_gms2, stats)
    get_stats_a_from_b_3p(labels_verified, labels_sbsp, stats)
    get_stats_a_from_b_3p_by_upstream(labels_verified, labels_ncbi, stats)

    # SN SP
    get_stats_sn_sp(labels_verified, labels_sbsp, stats)
    get_stats_sn_sp(labels_verified, labels_ncbi, stats)
    get_stats_sn_sp(labels_verified, labels_gms2, stats)

    # Stats: GMS2=SBSP Accuracy on verified
    get_stats_sn_sp(labels_verified, labels_sbsp_eq_gms2, stats)

    # stats by support
    get_stats_sn_sp_by_support(labels_verified, labels_sbsp, stats, "SBSP")

    # stats by support
    get_stats_sn_sp_by_support(labels_verified, labels_sbsp_eq_gms2, stats, "GMS2=SBSP")

    # stats by steps combinations
    get_stats_sn_sp_by_step_group(labels_verified, labels_sbsp, stats, "SBSP")

    # stats by steps combinations
    get_stats_sn_sp_by_step_group(labels_verified, labels_sbsp_eq_gms2, stats, "GMS2=SBSP")

    return stats


def df_to_pf_csv(df, pf_csv):
    # type: (pd.DataFrame, str) -> None
    df.to_csv(pf_csv, index=None)


def mk_pf(name, *prefixes):
    # type: (str, List[Any]) -> str
    out = ""
    if prefixes is not None and len(prefixes) > 0:
        non_empty_prefixes = [str(x) for x in prefixes if len(str(x)) > 0]
        "_".join(non_empty_prefixes)
        out += "_"

    out += name
    return out


def analyze_by_support(df, pd_work, fn_prefix, tag):
    # type: (pd.DataFrame, str, str, str) -> None

    list_df = list()
    for index in df.index:
        curr_df = pd.DataFrame(df.at[index, "by_support_{}".format(tag)])
        curr_df["Genome"] = df.at[index, "Genome"]

        if df.at[index, "Genome"] in {"A. pernix", "Synechocystis"}:
            continue

        list_df.append(curr_df)

    df_acc = pd.concat(list_df)

    sns.lineplot(
        df_acc, "Min Support", "Percentage 3p match: Verified from {}".format(tag), hue="Genome", figure_options=FigureOptions(
            title="Percentage of verified genes predicted\nby {}".format(tag),
            ylabel="Percentage",
            save_fig=next_name(pd_work),
            ylim=[None, 100.5]
        )
    )

    sns.lineplot(
        df_acc, "Min Support", "Percentage 5p-3p match: Verified from {}".format(tag), hue="Genome", figure_options=FigureOptions(
            title="Percentage of predicted {} genes\nwith correct 5' end".format(tag),
            ylabel="Percentage of 5p-3p match",
            save_fig=next_name(pd_work),
            ylim=[90, 100.5]
        )
    )

def analyze_by_step_group(df, pd_work, fn_prefix, tag):
    # type: (pd.DataFrame, str, str, str) -> None

    list_df = list()
    for index in df.index:
        curr_df = pd.DataFrame(df.at[index, "by_step_group_{}".format(tag)])
        curr_df["Genome"] = df.at[index, "Genome"]

        if df.at[index, "Genome"] in {"A. pernix", "Synechocystis"}:
            continue

        list_df.append(curr_df)

    df_acc = pd.concat(list_df)

    sns.catplot(
        df_acc, "Step Group", "Percentage 3p match: Verified from {}".format(tag), hue="Genome", kind="point", figure_options=FigureOptions(
            title="Percentage 3p match versus minimum support",
            ylabel="Percentage of 3p match",
            save_fig=next_name(pd_work),
            ylim=[None, 100.5]
        ),
    )

    sns.catplot(
        df_acc, "Step Group", "Percentage 5p-3p match: Verified from {}".format(tag), kind="point", hue="Genome", figure_options=FigureOptions(
            title="Percentage 5p-3p match versus minimum support",
            ylabel="Percentage of 5p-3p match",
            save_fig=next_name(pd_work),
            ylim=[90, 100.5]
        ),
    )

    print(df_acc.to_string())

def short_name(name):
    # type: (str) -> str
    components = name.split(maxsplit=1)

    if len(components) == 1:
        return name

    return "{}. {}".format(
        components[0][0], components[1]
    )


def print_csvs(env, df, **kwargs):
    # type: (Environment, pd.DataFrame, Dict[str, Any]) -> None

    fn_prefix = get_value(kwargs, "fn_prefix", "", default_if_none=True)
    pd_work = env["pd-work"]

    df["Genome"] = df["Genome"].apply(short_name)

    num = 0

    df_to_pf_csv(
        df[["Genome", "NCBI", "Verified",
            "Number 3p match: Verified from NCBI", "Percentage 3p match: Verified from NCBI",
            "Number 5p-3p match: Verified from NCBI", "Percentage 5p-3p match: Verified from NCBI"
            ]],
        next_name(pd_work, ext="csv")

    )
    num += 1

    df_to_pf_csv(
        df[["Genome", "GMS2", "Verified",
            "Number 3p match: Verified from GMS2", "Percentage 3p match: Verified from GMS2",
            "Number 5p-3p match: Verified from GMS2", "Percentage 5p-3p match: Verified from GMS2"
            ]],
        next_name(pd_work, ext="csv")    )
    num += 1

    df_to_pf_csv(
        df[["Genome", "SBSP", "Verified",
            "Number 3p match: Verified from SBSP", "Percentage 3p match: Verified from SBSP",
            "Number 5p-3p match: Verified from SBSP", "Percentage 5p-3p match: Verified from SBSP"
            ]],
        next_name(pd_work, ext="csv")    )
    num += 1

    df_to_pf_csv(
        df[["Genome", "Verified", "GMS2=SBSP",
            "Number 3p match: Verified from GMS2=SBSP", "Percentage 3p match: Verified from GMS2=SBSP",
            "Number 5p-3p match: Verified from GMS2=SBSP", "Percentage 5p-3p match: Verified from GMS2=SBSP"
            ]],
        next_name(pd_work, ext="csv")    )
    num += 1

    # # by support
    analyze_by_support(df, pd_work, fn_prefix, "SBSP")
    analyze_by_support(df, pd_work, fn_prefix, "GMS2=SBSP")

    analyze_by_step_group(df, pd_work, fn_prefix, "SBSP")
    analyze_by_step_group(df, pd_work, fn_prefix, "GMS2=SBSP")


def analyze_predictions_on_verified_genes_for_genome_list(env, gil, gcfid_to_pd_sbsp, **kwargs):
    # type: (Environment, GenomeInfoList, Dict[str, str], Dict[str, Any]) -> None

    fn_prefix = get_value(kwargs, "fn_prefix", "", default_if_none=True)

    info_per_gcfid = dict()

    for gi in gil:
        gcfid = gi.name
        # if gcfid != "Escherichia_coli_K_12_substr__MG1655_uid57779":
        #     continue
        try:
            pd_sbsp = gcfid_to_pd_sbsp[gcfid]
            info_per_gcfid[gcfid] = analyze_predictions_on_verified_genes(env, gi, pd_sbsp, **kwargs)
            info_per_gcfid[gcfid]["Genome"] = gi.attributes["name"]
            # info_per_gcfid["gi"] = gi
        except KeyError:
            logger.warning("Couldn't get SBSP directory for: {}".format(gcfid))

    list_stats = [info_per_gcfid[x] for x in info_per_gcfid.keys()]

    df = pd.DataFrame(list_stats)

    print_csvs(env, df, **kwargs)

    # df.to_csv(os_join(env["pd-work"], "{}summary.csv".format(fn_prefix)))


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None

    gil = GenomeInfoList.init_from_file(args.pf_genome_list)
    df = pd.read_csv(args.pf_gcfid_to_pd_sbsp)

    gcfid_to_pd_sbsp = {
        x["gcfid"]: x["pd-sbsp"] for _, x in df.iterrows()
    }

    analyze_predictions_on_verified_genes_for_genome_list(env, gil, gcfid_to_pd_sbsp,
                                                          fn_prefix=args.fn_prefix)


if __name__ == "__main__":
    main(my_env, parsed_args)
