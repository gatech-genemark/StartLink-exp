# Karl Gemayel
# Georgia Institute of Technology
#
# Created: 3/17/20
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

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from sbsp_general.general import os_join
from sbsp_general.labels import Labels, Label
from sbsp_general.shelf import compute_gc_from_file, map_key_3p_to_label
from sbsp_io.labels import read_labels_from_file

parser = argparse.ArgumentParser("Description of driver.")

parser.add_argument('--pf-genome-lists', required=True, nargs="+")
parser.add_argument('--list-names', required=True, nargs="+")
parser.add_argument('--dn-tools', required=True, nargs="+")
parser.add_argument('--tool-names', required=False, nargs="+")
parser.add_argument('--pf-output', required=True)

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

def read_labels_for_multiple_tools(env, gi, list_dn_tools, list_tool_names):
    # type: (Environment, GenomeInfo, List[str], List[str]) -> Dict[str, Labels]

    common_options = {"shift": 0, "ignore_frameshifted": True, "ignore_partial": True, "ignore_hypothetical": True}

    labels = dict()
    for name, dn_tool in zip(list_dn_tools, list_tool_names):
        pf_labels = os_join(env["pd-runs"], gi.name, dn_tool, f"{dn_tool}.gff")
        labels[name] =  read_labels_from_file(pf_labels, name="SBSP", **common_options)

    return labels


def _single_analysis(indexed_labels):
    # type: (Dict[str, Dict[str, Label]]) -> Dict[str, Any]
    return {
        name: len(lab) for name, lab in indexed_labels.items()
    }


def _pairwise_analysis(indexed_labels):
    # type: (Dict[str, Dict[str, Label]]) -> Dict[str, Any]

    result = dict()
    unique_names = sorted(set(indexed_labels.keys()))

    for i in range(len(unique_names)-1):
        n1 = unique_names[i]
        for j in range(i+1, len(unique_names)):
            n2 = unique_names[j]

            # find keys by 3p
            common_3p = set(indexed_labels[n1].keys()).intersection(indexed_labels[n2].keys())
            num_common_5p = sum(
                1 for c in common_3p if indexed_labels[n1][c].get_5prime() == indexed_labels[n2][c].get_5prime()
            )

            result.update({
                f"3p:{n1}={n2}": len(common_3p),
                f"5p:{n1}={n2}": num_common_5p
            })
    
    return result


def _all_together_analysis(indexed_labels):
    # type: (Dict[str, Dict[str, Label]]) -> Dict[str, Any]

    result = dict()
    sorted_names = sorted(indexed_labels.keys())

    common_3p = None
    tag = None

    for name in sorted_names:
        index = indexed_labels[name]

        if common_3p is None:
            common_3p = set(index.keys())
            tag = f"{name}"
        else:
            common_3p.intersection(index.keys())
            tag += f"={name}"

    def all_equal(items):
        # type: (List[Any]) -> bool
        for l_i in range(len(items)-1):
            for l_j in range(l_i + 1, len(items)):
                if items[l_i] != items[l_j]:
                    return False
        return True

    num_common_5p = sum(
        1 for c in common_3p if all_equal([indexed_labels[n][c].get_5prime() for n in indexed_labels.keys()])
    )

    return {
        f"3p:{tag}": len(common_3p),
        f"5p:{tag}": num_common_5p
    }


def stats_for_gi(env, gi, list_dn_tools, list_tool_names):
    # type: (Environment, GenomeInfo, List[str], List[str]) -> Dict[str, Any]
    result = {
        "Genome": gi.name
    }
    # compute GC
    result["GC"] = compute_gc_from_file(os_join(env["pd-data"], gi.name, "sequence.fasta"))

    # read all labels
    labels = read_labels_for_multiple_tools(env, gi, list_dn_tools, list_tool_names)

    # indexed labels
    indexed_labels = {
        name: map_key_3p_to_label(lab) for name, lab in labels.items()
    }

    # single analysis
    result.update(_single_analysis(indexed_labels))

    # pairwise analysis
    result.update(_pairwise_analysis(indexed_labels))
    
    # all together
    result.update(_all_together_analysis(indexed_labels))


    return result


def stats_for_gil(env, gil, list_dn_tools, list_tool_names):
    # type: (Environment, GenomeInfoList, List[str], List[str]) -> pd.DataFrame

    list_entries = list()
    for gi in gil:
        entry = stats_for_gi(env, gi, list_dn_tools, list_tool_names)
        list_entries.append(entry)

    return pd.DataFrame(list_entries)


def stats_tools_5prime(env, list_gil, list_names, list_dn_tools, list_tool_names, pf_output):
    # type: (Environment, List[GenomeInfoList], List[str], List[str], List[str], str) -> None

    # for each gil
    list_df = list()
    for name, gil in zip(list_names, list_gil):
        df_tmp = stats_for_gil(env, gil, list_dn_tools, list_tool_names)

        df_tmp["Genome"] = name
        list_df.append(df_tmp)

    df = pd.concat(list_df, sort=False)

    df.to_csv(pf_output, index=False)


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None
    list_gil = [
        GenomeInfoList.init_from_file(pf) for pf in args.pf_genome_lists
    ]

    list_names = args.list_names

    if len(list_names) != len(list_gil):
        raise Warning(f"Names and genome lists must have the same length {len(list_names)} != {len(list_gil)}")

    list_dn_tools = args.dn_tools
    list_tool_names = args.tool_names
    if len(list_dn_tools) != len(list_tool_names):
        raise Warning(f"Tools and dirs must have the same length {len(list_dn_tools)} != {len(list_tool_names)}")

    stats_tools_5prime(env, list_gil, list_names, list_dn_tools, list_tool_names, args.pf_output)



if __name__ == "__main__":
    main(my_env, parsed_args)
