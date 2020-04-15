import logging
import pandas as pd
from typing import *
from sbsp_general.msa_2 import MSAType, MultipleSeqAlignment, MSASinglePointMarker
from sbsp_general.general import get_value
from sbsp_options.sbsp import SBSPOptions
from sbsp_alg.msa_2 import get_positions_of_query_candidate_starts_in_msa_with_class
from sbsp_ml import msa_features_2 as msa_features

logger = logging.getLogger(__name__)


def try_and_except(func, *args, **kwargs):
    # type: (Callable, List[Any], Dict[str, Any]) -> Any

    return func(*args, **kwargs)

def generate_features_for_position(msa_t, position, msa_options, **kwargs):
    # type: (MSAType, int, SBSPOptions, Dict[str, Any]) -> Dict[str, Any]

    primary_defaults_for_features = get_value(kwargs, "primary_defaults_for_features", None)
    secondary_default_for_all = get_value(kwargs, "secondary_default_for_all", None)

    try:
        f_upstream_identity = msa_features.compute_upstream_score(
            msa_t, position, msa_options, scoring_function=msa_features.ScoringMatrix("identity")
        )
    except ValueError:
        raise ValueError("Cannot compute upstream score")

    try:
        f_upstream_identity_all_pairs = msa_features.compute_upstream_score(
            msa_t, position, msa_options, scoring_function=msa_features.ScoringMatrix("identity"),
            score_on_all_pairs=True
        )
    except ValueError:
        raise ValueError("Cannot compute upstream score with all pairs")

    try:
        f_upstream_blosum = msa_features.compute_upstream_score(
            msa_t, position, msa_options, scoring_function=msa_features.ScoringMatrix("blosum45")
        )
    except ValueError:
        raise ValueError("Cannot compute downstream score")

    try:
        f_upstream_blosum_all_pairs = msa_features.compute_upstream_score(
            msa_t, position, msa_options, scoring_function=msa_features.ScoringMatrix("blosum45"),
            score_on_all_pairs=True
        )
    except ValueError:
        raise ValueError("Cannot compute downstream score")

    try:
        f_downstream_identity = msa_features.compute_downstream_score(
            msa_t, position, msa_options, scoring_function=msa_features.ScoringMatrix("identity")
        )
    except ValueError:
        raise ValueError("Cannot compute downstream identity")

    try:
        f_downstream_blosum = msa_features.compute_downstream_score(
            msa_t, position, msa_options, scoring_function=msa_features.ScoringMatrix("blosum45")
        )
    except ValueError:
        raise ValueError("Cannot compute downstream blosum")

    try:
        f_5prime_score = msa_features.compute_5prime_score(
            msa_t, position, msa_options
        )
    except ValueError:
        raise ValueError("Cannot compute 5-prime score")

    try:
        f_5prime_penalized_score = msa_features.compute_simple_saas(msa_t, position)
    except ValueError:
        raise ValueError("Cannot compute 5-prime penalized score")

    return {
        "upstream-identity": f_upstream_identity,
        "upstream-identity-all-pairs": f_upstream_identity_all_pairs,
        "upstream-blosum": f_upstream_blosum,
        "upstream-blosum-all-pairs": f_upstream_blosum_all_pairs,
        "downstream-identity": f_downstream_identity,
        "downstream-blosum": f_downstream_blosum,
        "5prime-score": f_5prime_score,
        "5prime-penalized-score": f_5prime_penalized_score
    }




def generate_features_for_msa(msa_t, msa_options, **kwargs):
    # type: (MSAType, SBSPOptions, Dict[str, Any]) -> pd.DataFrame

    # TODO: Clean up and use MSAType only


    # Sketch:
    # Get candidate starts, with label relative to reference
    # Compute features for each candidate

    candidate_positions_in_msa, candidate_class = get_positions_of_query_candidate_starts_in_msa_with_class(
        msa_t, **kwargs
    )

    list_features = list()

    for curr_pos, curr_class in zip(candidate_positions_in_msa, candidate_class):

        try:
            curr_features = generate_features_for_position(msa_t, curr_pos, msa_options)
            curr_features["class"] = curr_class
            list_features.append(curr_features)

        except ValueError:
            pass

    return pd.DataFrame(list_features)




def generate_features_for_msa_from_file(pf_msa, msa_options, **kwargs):
    # type: (str, SBSPOptions, Dict[str, Any]) -> pd.DataFrame

    # read msa from file
    msa_t = MSAType.init_from_file(pf_msa)

    # select start
    df_features = generate_features_for_msa(msa_t, msa_options, **kwargs)
    df_features["pf-msa-output"] = pf_msa
    return df_features


def generate_features_for_msa_from_df(df_data, msa_options, **kwargs):
    # type: (pd.DataFrame, SBSPOptions, Dict[str, Any]) -> pd.DataFrame

    df_features = pd.DataFrame()

    for index, row in df_data.groupby("pf-msa-output", as_index=False).agg("first").iterrows():
        pf_msa = row["pf-msa-output"]
        curr_df_features = generate_features_for_msa_from_file(pf_msa, msa_options, **kwargs)

        df_features = df_features.append(curr_df_features, ignore_index=True)

    return df_features

