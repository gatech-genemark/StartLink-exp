# Karl Gemayel
# Georgia Institute of Technology
#
# Created: 3/16/20
import datetime
import logging
import argparse
from math import ceil

import pandas as pd
from copy import copy
from typing import *

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import sbsp_log  # runs init in sbsp_log and configures logger

# Custom imports
from sbsp_container.genome_list import GenomeInfoList
from sbsp_general import Environment

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from sbsp_general.general import get_value

parser = argparse.ArgumentParser("Choose a sample of genomes from a list.")

parser.add_argument('--pf-genome-list', required=True, help="Genome list")
parser.add_argument('--pf-output', required=True, help="Output file")
parser.add_argument('-n', required=True, type=int, help="Number of genomes")

parser.add_argument('--pf-gcfid-to-ignore', required=False)

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


def sample_from_genome_list(gil, n, **kwargs):
    # type: (GenomeInfoList, int, Dict[str, Any]) -> GenomeInfoList

    gcfid_to_ignore = get_value(kwargs, "gcfid_to_ignore", None)
    if gcfid_to_ignore is None:
        gcfid_to_ignore = set()

    if len(gil) < n:
        return copy(gil)

    for gi in gil:
        if len(gi.attributes["annotation_date"].strip()) == 0:
            gi.attributes["annotation_date"] = "01/01/0001"
        # datetime.datetime.strptime(gi.attributes["annotation_date"], "%m/%d/%Y")

    df = pd.DataFrame({
        "parent_id": gi.attributes["parent_id"],
        "gc": float(gi.attributes["gc"]),
        "gc_int": int(float(gi.attributes["gc"])),
        "num_genes": int(gi.attributes["num_genes"]),
        "annotation_date": datetime.datetime.strptime(gi.attributes["annotation_date"], "%m/%d/%Y"),
        "gcfid": gi.name,
        "gi": gi
    } for gi in gil)

    if len(gcfid_to_ignore) > 0:
        df = df[~df["gcfid"].isin(gcfid_to_ignore)]

    # remove anything before 2020
    df = df[df["annotation_date"] >= datetime.datetime.strptime("01/01/2020", "%m/%d/%Y")]
    df = df[df["num_genes"] < 6000]


    num_unique_parents = len(df["parent_id"].unique())
    num_per_parent = ceil(n / float(num_unique_parents))

    # select roughly "num_per_parent" per parent (to normalize oversampling of some genomes)

    df_normalized = pd.DataFrame(columns=df.columns)

    # get all groups in ascending order (by group size)
    list_df_group = sorted(
        [df_group for key, df_group in df.groupby("parent_id", as_index=False)],
        key=lambda x: len(x)
    )       # type: List[pd.DataFrame]

    list_normalized = list()
    missed = 0          # in case a group doesn't have enough data, allow gathering from other groups
    for df_group in list_df_group:

        if len(df_group) <= num_per_parent:
            list_normalized += [r for _, r in df_group.iterrows()]

            missed += num_per_parent - len(df_group)

        # otherwise, sample as much as possible
        else:

            recovered_from_missed = 0

            # cannot recover all
            if len(df_group) < num_per_parent + missed:
                amount_to_sample = len(df_group)
                recovered_from_missed = amount_to_sample - num_per_parent
            else:
                amount_to_sample = num_per_parent + missed
                recovered_from_missed = missed

            missed -= recovered_from_missed

            df_sampled = df_group.sample(amount_to_sample)

            list_normalized += [r for _, r in df_sampled.iterrows()]

    # loop over parents by size

    df_normalized = pd.DataFrame(list_normalized)
    df_sampled = df_normalized.sample(min(n, len(df_normalized)))

    return GenomeInfoList([r["gi"] for _, r in df_sampled.iterrows()])


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None

    gil = GenomeInfoList.init_from_file(args.pf_genome_list)

    gcfid_to_ignore = None
    if args.pf_gcfid_to_ignore is not None:
        df_gti = pd.read_csv(args.pf_gcfid_to_ignore)
        gcfid_to_ignore = set(df_gti["gcfid"])

    gil_n = sample_from_genome_list(gil, n=args.n, gcfid_to_ignore=gcfid_to_ignore)
    gil_n.to_file(args.pf_output)


if __name__ == "__main__":
    main(my_env, parsed_args)
