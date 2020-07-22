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
from sbsp_alg.sbsp_steps import compute_conservation_in_region, count_number_of_5prime_candidates_at_position
from sbsp_argparse.sbsp import add_sbsp_options
from sbsp_general import Environment

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from sbsp_container.msa import MSAType
from sbsp_ml.msa_features import ScoringMatrix
from sbsp_options.sbsp import SBSPOptions

parser = argparse.ArgumentParser("Description of driver.")

parser.add_argument('--pf-data', required=True, help="SBSP details output")
parser.add_argument('--pf-output', required=True, help="Output file")
add_sbsp_options(parser)


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

def get_block_positions_with_region_label(msa_t, **kwargs):
    # type: (MSAType, Dict[str, Any]) -> Dict[str, int]

    pos_annotated = msa_t.get_mark_position("ref")

    output = dict()     # type: Dict[str, int]

    def get_first_nongap(begin, end, direction, skip_nongap=0):
        # type: (int, int, str, int) -> Union[int, None]
        curr = begin
        non_gap_skipped = 0

        while curr != end:
            if msa_t[0][curr] != "-":
                if non_gap_skipped >= skip_nongap:
                    return curr
                else:
                    non_gap_skipped += 1
            if direction == "downstream":
                curr += 1
            else:
                curr -= 1
        return None


    # get close and far downstreams
    output["Dc"] = get_first_nongap(pos_annotated + 1, msa_t.alignment_length(), "downstream")
    output["Df"] = get_first_nongap(pos_annotated + 1, msa_t.alignment_length(), "downstream", skip_nongap=20)

    output["Uc"] = get_first_nongap(pos_annotated - 1, -1, "upstream")
    output["Uf"] = get_first_nongap(pos_annotated - 1, -1, "upstream", skip_nongap=20)

    return output


def generate_block_features_for_starts_for_msa(msa_t, sbsp_options, **kwargs):
    # type: (MSAType, SBSPOptions, Dict[str, Any]) -> List[Dict[str, Any]]

    name_to_pos = get_block_positions_with_region_label(msa_t, **kwargs)
    scorer = ScoringMatrix("identity")
    common = {"skip_gaps_in_query": True, "score_on_all_pairs": True}

    output = list()
    for name, pos in name_to_pos.items():
        if pos is None:
            continue
        if name in {"Dc", "Df"}:
            try:
                cons = compute_conservation_in_region(msa_t, pos, pos+10, scorer, direction="downstream", **common)
            except ValueError:
                continue
        else:
            try:
                cons = compute_conservation_in_region(msa_t, pos-11, pos, scorer, direction="upstream", **common)
            except ValueError:
                continue

        output.append({"region": name, "score": cons, "type": "block"})

    return output


def get_5prime_positions_with_region_label(msa_t, **kwargs):
    # type: (MSAType, Dict[str, Any]) -> List[Dict[str, int]]
    pos_annotated = msa_t.get_mark_position("ref")

    output = list()  # type: List[Dict[str, int]]

    def get_first_x_5primes(begin, end, direction, skip_nongap=0):
        # type: (int, int, str, int) -> List[int]
        curr = begin
        non_gap_skipped = 0
        candidates = list()
        while curr != end:
            if msa_t[0][curr].upper():
                if non_gap_skipped >= skip_nongap:
                    candidates.append(curr)
                else:
                    non_gap_skipped += 1
            if direction == "downstream":
                curr += 1
            else:
                curr -= 1

            if len(candidates) >= 3:
                break
        return candidates

    # get close and far downstreams
    l_downstream = get_first_x_5primes(pos_annotated +1, msa_t.alignment_length(), "downstream", skip_nongap=5)
    l_upstream = get_first_x_5primes(pos_annotated -1 , -1, "upstream", skip_nongap=5)

    output += [{"region": "Downstream", "pos": x} for x in l_downstream]
    output += [{"region": "Upstream", "pos": x} for x in l_upstream]
    output += [{"region": "Annotated", "pos": pos_annotated}]

    return output


def generate_5prime_features_for_starts_for_msa(msa_t, sbsp_options, **kwargs):
    # type: (MSAType, SBSPOptions, Dict[str, Any]) -> List[Dict[str, Any]]
    list_region_pos = get_5prime_positions_with_region_label(msa_t, **kwargs)
    output = list()
    for region_pos in list_region_pos:
        number_5prime = count_number_of_5prime_candidates_at_position(msa_t, region_pos["pos"], sbsp_options)
        score = number_5prime / float(msa_t.number_of_sequences())

        output.append({
            "region": region_pos["region"], "score": score, "type": "5prime"
        })

    return output


def generate_features_for_starts_for_single_pf_msa(series, sbsp_options, **kwargs):
    # type: (pd.Series, SBSPOptions, Dict[str, Any]) -> List[Dict[str, Any]]

    msa_t = MSAType.init_from_file(series["pf-msa-output"])
    if msa_t.get_mark_position("ref") is None:
        return None

    # block features
    list_blk = generate_block_features_for_starts_for_msa(msa_t, sbsp_options, **kwargs)


    # five prime features
    list_5prime = generate_5prime_features_for_starts_for_msa(msa_t, sbsp_options, **kwargs)

    return list_blk + list_5prime


def generate_features_for_starts_from_msa(df, sbsp_options, pf_output, **kwargs):
    # type: (pd.DataFrame, SBSPOptions, str, Dict[str, Any]) -> None

    list_df = list()           # type: List[pd.DataFrame]

    for _, df_group in df.groupby("pf-msa-output", as_index=False):
        series = df_group.iloc[0]
        list_entries = generate_features_for_starts_for_single_pf_msa(series, sbsp_options, **kwargs)
        if list_entries is None:
            continue

        df_tmp = pd.DataFrame(list_entries)
        df_tmp["Genome"] = series["q-genome"]
        list_df.append(df_tmp)

    df_features = pd.concat(list_df)

    df_features.to_csv(pf_output, index=False)


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None
    sbsp_options = SBSPOptions.init_from_dict(env, vars(args))

    df_data = pd.read_csv(args.pf_data, header=0)
    generate_features_for_starts_from_msa(df_data, sbsp_options, args.pf_output)


if __name__ == "__main__":
    main(my_env, parsed_args)
