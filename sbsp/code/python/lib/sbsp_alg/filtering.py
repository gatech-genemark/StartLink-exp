from typing import *
import pandas as pd

from sbsp_general import Environment
from sbsp_general.general import get_value
from sbsp_io.general import write_string_to_file
from sbsp_options.sbsp import SBSPOptions


def filter_orthologs(env, pf_data, pf_output, **kwargs):
    # type: (Environment, str, str, Dict[str, Any]) -> str
    msa_options = get_value(kwargs, "msa_options", SBSPOptions(env))  # type: SBSPOptions
    pf_filter_stats = get_value(kwargs, "pf_filter_stats", None)
    filter_non_group_only = get_value(kwargs, "filter_non_group_only", True)

    from sbsp_alg.msa import filter_df, print_filter_stats_to_file, print_filter_stats

    column_distance = "k2p-distance"
    if msa_options.safe_get("column-distance"):
        column_distance = msa_options.safe_get("column-distance")

    df = pd.read_csv(pf_data, header=0)

    filter_stats = {}
    df = filter_df(df, msa_options, filter_stats=filter_stats, filter_non_group_only=filter_non_group_only,
                   column_distance=column_distance)

    # df = df_add_labeled_sequences(env, df,
    #                               source="both",
    #                               suffix_coordinates=suffix_coordinates,
    #                               suffix_gene_sequence=tag_msa,
    #                               upstream_length_nt=upstream_length_nt,
    #                               downstream_length_nt=downstream_length_nt)

    if pf_filter_stats:
        print_filter_stats_to_file(filter_stats, pf_filter_stats)
    else:
        print_filter_stats(filter_stats)

    df.to_csv(pf_output, index=False)

    write_string_to_file("\n".join(set(df["q-3prime"].unique())), pf_output + "_map")
    return pf_output
