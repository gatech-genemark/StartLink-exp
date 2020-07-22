import logging
import os
from typing import *

from Bio.Align.Applications import ClustalOmegaCommandline

from sbsp_container.msa import MSAType
from sbsp_general.general import get_value
from sbsp_io.general import write_string_to_file, remove_p

logger = logging.getLogger(__name__)


def write_sequence_list_to_fasta_file(sequences, pf_sequences):
    # type: (List[Seq], str) -> None

    data = ""
    for i in range(len(sequences)):
        data += ">{}\n{}\n".format(i, sequences[i])

    write_string_to_file(data, pf_sequences)


def run_msa_on_sequence_file(pf_fasta, sbsp_options, pf_msa, **kwargs):
    # type: (str, SBSPOptions, str, Dict[str, Any]) -> None

    gapopen = sbsp_options.safe_get("msa-gapopen")
    gapext = sbsp_options.safe_get("msa-gapext")

    num_processors = get_value(kwargs, "num_processors", None)
    output_order = get_value(kwargs, "outputorder", "input-order")
    gapopen = get_value(kwargs, "gapopen", None)

    # clustalw_cline = ClustalwCommandline(
    #     "clustalw2", infile=pf_fasta, outfile=pf_msa,
    #     gapopen=gapopen,
    #     gapext=gapext,
    #     outorder="input"
    # )

    logger.debug("Number of processors for MSA: {}".format(num_processors))
    other_options = dict()
    if num_processors is not None:
        other_options["threads"] = num_processors
    # if gapopen is not None:
    #     other_options["gapopen"] = gapopen

    clustalw_cline = ClustalOmegaCommandline(
        "clustalo", infile=pf_fasta, outfile=pf_msa,
        # gapopen=gapopen,
        # gapext=gapext,
        outputorder=output_order,
        force=True,
        outfmt="clustal",
        **other_options
    )

    clustalw_cline()


def run_msa_on_sequences(env, sequences, sbsp_options, **kwargs):
    # type: (Environment, List[Seq], SBSPOptions, Dict[str, Any]) -> MSAType

    pd_work = env["pd-work"]
    fn_tmp_prefix = get_value(kwargs, "fn_tmp_prefix", "", default_if_none=True)

    # write sequences to file
    pf_fasta = os.path.join(pd_work, "{}tmp_sequences.fasta".format(fn_tmp_prefix))
    remove_p(pf_fasta)
    write_sequence_list_to_fasta_file(sequences, pf_fasta)

    # run msa
    pf_msa = os.path.join(pd_work, "{}tmp_msa.txt".format(fn_tmp_prefix))
    run_msa_on_sequence_file(pf_fasta, sbsp_options, pf_msa, **kwargs)

    msa_t = MSAType.init_from_file(pf_msa)

    remove_p(pf_msa, pf_fasta)

    return msa_t