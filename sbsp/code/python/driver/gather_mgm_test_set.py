# Karl Gemayel
# Georgia Institute of Technology
#
# Created: 3/17/20
import logging
import argparse
import os

import pandas as pd
from typing import *

# noinspection All
from Bio.Seq import Seq
from tqdm import tqdm

import pathmagic

# noinspection PyUnresolvedReferences
import sbsp_log  # runs init in sbsp_log and configures logger

# Custom imports
from sbsp_container.genome_list import GenomeInfoList, GenomeInfo
from sbsp_container.gms2_mod import GMS2Mod
from sbsp_general import Environment

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from sbsp_general.MotifModel import MotifModel
from sbsp_general.general import os_join, get_value
from sbsp_general.labels import Labels, Label
from sbsp_io.general import remove_p
from sbsp_io.labels import read_labels_from_file
from sbsp_io.sequences import read_fasta_into_hash

parser = argparse.ArgumentParser("Gather sequences, and motif scores.")

parser.add_argument('--pf-genome-list', required=True, help="List of genomes")
parser.add_argument('--pf-output', required=True, help="Output file")


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


def append_data_frame_to_csv(df, pf_output):
    # type: (pd.DataFrame, str) -> None
    if df is not None and len(df) > 0:
        if not os.path.isfile(pf_output):
            df.to_csv(pf_output, index=False)
        else:
            df.to_csv(pf_output, mode="a", index=False, header=False)


def extract_upstream_sequences(labels, sequences, **kwargs):
    # type: (Labels, Dict[str, Seq], Dict[str, Any]) -> List[Tuple[Label, Seq]]

    upstream_length_nt = get_value(kwargs, "upstream_length", 20, default_if_none=True)
    sequence_type = get_value(kwargs, "sequence_type", choices=["aa", "nt"], default="nt",
                              default_if_none=True)
    reverse_complement = get_value(kwargs, "reverse_complement", default=True)

    result = list()

    for l in labels:
        if l.seqname() in sequences.keys():
            # positive strand
            frag = None
            if l.strand() == "+":
                frag_left = l.left()-upstream_length_nt
                frag_right = l.left()-1
                if frag_left >= 0 and frag_right < len(sequences[l.seqname()]):
                    frag = sequences[l.seqname()][frag_left:frag_right+1]
            # negative strand
            else:
                frag_left = l.right() + 1
                frag_right = l.right() + upstream_length_nt
                if frag_left >= 0 and frag_right < len(sequences[l.seqname()]):
                    frag = sequences[l.seqname()][frag_left:frag_right + 1]
                    if reverse_complement:
                        frag = frag.reverse_complement()

            if frag is not None:
                result.append((l, frag))
    return result


def gather_upstream_sequences_for_genome(env, gi, **kwargs):
    # type: (Environment, GenomeInfo, Dict[str, Any]) -> pd.DataFrame

    list_entries = list()       # type: List[Dict[str, Any]]

    # read sequences
    pf_sequences = os_join(env["pd-data"], gi.name, "sequence.fasta")
    pf_labels = os_join(env["pd-runs"], gi.name, "gms2", "gms2.gff")

    sequences = read_fasta_into_hash(pf_sequences)
    labels = read_labels_from_file(pf_labels)

    upstream_info = extract_upstream_sequences(labels, sequences)

    for info in upstream_info:
        label = info[0]      # type: Label
        frag = info[1]      # type: Seq

        list_entries.append({
            "left": label.left() + 1,
            "right": label.right() + 1,
            "strand": label.strand() + 1,
            "GCFID": gi.name,
            "Accession": label.seqname(),
            "upstream_nt": str(frag)
        })

    return pd.DataFrame(list_entries)

def create_motif_model_from_gms2_model(mod, key):
    # type: (GMS2Mod, str) -> Union[MotifModel, None]
    if f"{key}_MAT" not in mod.items:
        return None

    value_motif = mod.items[f"{key}_MAT"]
    value_spacer = mod.items[f"{key}_POS_DISTR"] if f"{key}_POS_DISTR" in mod.items else None

    return MotifModel(value_motif, value_spacer)


def gather_mgm_test_set_for_genome(env, gi, **kwargs):
    # type: (Environment, GenomeInfo, Dict[str, Any]) -> pd.DataFrame

    # get upstream sequences
    df = gather_upstream_sequences_for_genome(env, gi)

    pf_mod = os_join(env["pd-runs"], gi.name, "gms2", "GMS2.mod")
    mod = GMS2Mod.init_from_file(pf_mod)

    m_rbs = create_motif_model_from_gms2_model(mod, "RBS")
    m_promoter = create_motif_model_from_gms2_model(mod, "PROMOTER")

    names = ["RBS", "PROMOTER"]
    models = [m_rbs, m_promoter]

    # add score columns to dataframe
    score_column_names = [x + y for x in ["RBS", "Promoter"] for y in ["motif", "spacer", "both"]]
    df.reindex(columns=[*(df.columns.tolist()+score_column_names)], fill_value=None)

    for idx in df.index:
        frag = df.at[idx, "upstream_nt"]

        for name, model in zip(names, models):
            if model is not None:
                for c in ["motif", "spacer", "both"]:
                    result = model.find_best_position_and_score(frag, component=c)
                    pos = result[0]
                    score = result[1]
                    df.at[idx, f"{name}_{c}_score"] = score
                    df.at[idx, f"{name}_{c}_position"] = len(frag) - pos - model.motif_width()

    return df


def gather_mgm_test_set(env, gil, pf_output, **kwargs):
    # type: (Environment, GenomeInfoList, str, Dict[str, Any]) -> None

    remove_p(pf_output)     # start clean

    for gi in tqdm(gil, total=len(gil)):

        df = gather_mgm_test_set_for_genome(env, gi, **kwargs)
        append_data_frame_to_csv(df, pf_output)


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None

    gil = GenomeInfoList.init_from_file(args.pf_genome_list)
    gather_mgm_test_set(env, gil, args.pf_output)


if __name__ == "__main__":
    main(my_env, parsed_args)
