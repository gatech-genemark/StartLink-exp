# Karl Gemayel
# Georgia Institute of Technology
#
# Created: 3/17/20
import re
import logging
import argparse
from subprocess import CalledProcessError
import matplotlib.pyplot as plt
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
from sbsp_general.general import os_join, get_value
from sbsp_general.labels_comparison_detailed import LabelsComparisonDetailed
from sbsp_general.shelf import run_gms2_prediction_with_model, fix_names, next_name
from sbsp_io.general import mkdir_p
from sbsp_io.labels import read_labels_from_file
from sbsp_viz import sns
from sbsp_viz.general import FigureOptions

parser = argparse.ArgumentParser("Description of driver.")

parser.add_argument('--pf-genome-list', required=True, help="Genome information list")

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

def key_to_gms2_tags(key):
    # type: (str) -> Set[str]

    return {
        "Start Codons": {"ATG", "GTG", "TTG"},
        "Stop Codons ": {"TAA", "TGA", "TAG"},
        "Start Context": {"SC_RBS", "SC_PROMOTER"},
        "RBS": {"RBS"},
        "Promoter": {"PROMOTER"},
        "Native": {"TO_NATIVE", "TO_MGM"}
    }[key]


def clean_up_start_context(mod_string, t, delete_sc):
    # type: (str, str, bool) -> str
    mod_string = re.sub(r"\$" + t + r" [^\n]+", "", mod_string)
    mod_string = re.sub(r"\$" + t + r"_ORDER [^\n]+", "", mod_string)
    mod_string = re.sub(r"\$" + t + r"_WIDTH [^\n]+", "", mod_string)
    mod_string = re.sub(r"\$" + t + r"_MARGIN [^\n]+", "", mod_string)
    mod_string = re.sub(r"\$" + t + r"_MAT[^\$]+", "", mod_string)

    if not delete_sc:
        mod_string += f"\n${t} 1\n${t}_ORDER 0\n${t}_WIDTH 1\n${t}_MARGIN -1\n${t}_MAT\nA 0.25\nC 0.25\nG 0.25\nT 0.25\n"
    return mod_string


def turn_off_components(pf_mod_original, pf_new_mod, components_off, native_coding_off=True):
    # type: (str, str, Iterable[str], bool) -> None

    with open (pf_mod_original, "r") as f:
        mod_string = f.read()

        turn_off_sc_for_rbs = "RBS" in components_off
        turn_off_sc_for_promoter = "Promoter" in components_off

        # change model based on components
        for coff in components_off:
            tags = key_to_gms2_tags(coff)

            if coff in {"Start Codons", "Stop Codons"}:
                for t in tags:
                    if re.findall(r"\$" + t + r"[\s\n]", mod_string):
                        mod_string = re.sub(r"\$" + t + r"\s+\d+\.\d+", f"${t} 0.333", mod_string)
            elif coff in {"RBS", "Promoter"}:
                for t in tags:
                    if re.findall(r"\$" + t + r"[\s\n]", mod_string):
                        mod_string = re.sub(r"\$" + t + r"\s+1", f"${t} 0", mod_string)
            elif coff in {"Start Context"}:
                for t in tags:
                    if re.findall(r"\$" + t + r"[\s\n]", mod_string):
                        delete_sc = turn_off_sc_for_rbs if t == "SC_RBS" else turn_off_sc_for_promoter
                        mod_string = clean_up_start_context(mod_string, t, delete_sc)

        if native_coding_off:
            mod_string = re.sub(r"\$TO_NATIVE" + r"\s+\d+\.\d+", f"$TO_NATIVE 0.0", mod_string)
            mod_string = re.sub(r"\$TO_MGM" + r"\s+\d+\.\d+", f"$TO_MGM 1.0", mod_string)

        with open(pf_new_mod, "w") as f_out:
            f_out.write(mod_string)



def run_gms2_with_component_toggles_and_get_accuracy(env, gi, components_off, **kwargs):
    # type: (Environment, GenomeInfo, Set[str], Dict[str, Any]) -> Dict[str, Any]

    pf_mod_original = os_join(env["pd-runs"], gi.name, "gms2", "GMS2.mod")
    pf_reference = os_join(env["pd-data"], gi.name, "verified.gff")
    pf_sequence = os_join(env["pd-data"], gi.name, "sequence.fasta")
    pf_prediction = os_join(env["pd-work"], "prediction.gff")

    native_coding_off = get_value(kwargs, "native_coding_off", True)

    pf_new_mod = os_join(env["pd-work"], "model.mod")
    turn_off_components(pf_mod_original, pf_new_mod, components_off, native_coding_off=native_coding_off)

    done = False
    while not done:
        try:
            run_gms2_prediction_with_model(pf_sequence, pf_new_mod, pf_prediction)
            done = True
        except CalledProcessError:
            pass

    # compare with verified
    lcd = LabelsComparisonDetailed(read_labels_from_file(pf_reference), read_labels_from_file(pf_prediction))

    return {
        "Error": 100 - 100 * len(lcd.match_3p_5p('a')) / len(lcd.match_3p('a'))
    }


def component_in_model_file(env, gi, component):
    # type: (Environment, GenomeInfo, str) -> bool
    pf_mod = os_join(env["pd-runs"], gi.name, "gms2", "GMS2.mod")
    with open (pf_mod, "r") as f:
        mod_string = f.read()

        for t in key_to_gms2_tags(component):
            if re.findall(r"\$" + t + r"[\s\n]", mod_string):
            # if (r"$" + t + r"") in mod_string:
                return True
        return False


def analyze_gms2_components_on_verified_set_for_gi(env, gi):
    # type: (Environment, GenomeInfo) -> pd.DataFrame

    list_entries = list()

    start_components = {
        "Start Codons", "Start Context", "RBS", "Promoter",
    }

    pd_gi = os_join(env["pd-work"], gi.name)
    mkdir_p(pd_gi)

    # for each component to keep on
    for component_on in sorted(start_components) + ["MGM2*", "MGM", "GMS2"]:
        components_off = start_components.difference({component_on})

        if component_on == "MGM2*" or component_on == "GMS2":
            components_off = set()
        elif component_on == "MGM":
            pass
        elif not component_in_model_file(env, gi, component_on) and component_on not in {"MGM2*", "MGM", "GMS2"}:
            continue

        native_coding_off = False if component_on == "GMS2" else True

        pd_gi_component = os_join(pd_gi, component_on).replace(" ", "")
        mkdir_p(pd_gi_component)

        env_dup = env.duplicate({"pd-work": pd_gi_component})

        if component_on == "Start Context":
            component_on = {component_on}  # "rbs", "promoter"}
            components_off.remove("RBS")
            components_off.remove("Promoter")
        else:
            component_on = {component_on}


        results = run_gms2_with_component_toggles_and_get_accuracy(env_dup, gi, components_off,
                                                                   native_coding_off=native_coding_off)

        list_entries.append({
            "Genome": gi.name,
            "Component": next(iter(component_on)).replace("_", "-"),
            # **{con: True for con in component_on},                             # current component is on
            # **{coff: False for coff in components_off},     # all others are off
            **results
        })



    return pd.DataFrame(list_entries)


def analyze_gms2_components_on_verified_set(env, gil):
    # type: (Environment, GenomeInfoList) -> None

    # run different components
    list_df = list()
    for gi in gil:
        list_df.append(
            analyze_gms2_components_on_verified_set_for_gi(env, gi)
        )

    df = pd.concat(list_df, ignore_index=True, sort=False)
    df["Genome"] = df.apply(fix_names, axis=1)
    print(df.to_csv())


    fig, ax = plt.subplots(figsize=(12,4))
    sns.barplot(df, "Genome", "Error", hue="Component",
                ax=ax,
                figure_options=FigureOptions(
                    save_fig=next_name(env["pd-work"])
                ),
                sns_kwargs={
                    "hue_order": reversed(["GMS2", "MGM2*", "Start Context", "RBS", "Start Codons", "Promoter", "MGM"]),
                    "palette": CM.get_map("gms2_components")

                })




def main(env, args):
    # type: (Environment, argparse.Namespace) -> None

    gil = GenomeInfoList.init_from_file(args.pf_genome_list)

    analyze_gms2_components_on_verified_set(env, gil)


if __name__ == "__main__":
    main(my_env, parsed_args)
