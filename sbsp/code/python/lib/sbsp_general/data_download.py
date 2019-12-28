import os
import logging
import pandas as pd
from typing import *

from sbsp_container.genome_list import GenomeInfoList, GenomeInfo
from sbsp_general import Environment
from sbsp_general.general import get_value, run_shell_cmd
from sbsp_io.general import mkdir_p

logger = logging.getLogger(__name__)


def set_up_gcfid(gcfid_info, pd_output):
    # type: (Dict[str, Any], str) -> None

    # build name
    gcf = gcfid_info["assembly_accession"]
    acc = gcfid_info["asm_name"].replace(" ", "_")

    gcfid = "{}_{}".format(gcf, acc)

    pd_gcfid = os.path.join(pd_output, gcfid)
    pd_runs = os.path.join(pd_gcfid, "runs")
    mkdir_p(pd_gcfid)
    mkdir_p(pd_runs)

    ftplink = gcfid_info["ftp_path"]
    fn_sequence = "{}_genomic.fna".format(gcfid)
    fn_labels = "{}_genomic.gff".format(gcfid)

    pf_ftp_sequence = os.path.join(ftplink, "{}.gz".format(fn_sequence))
    pf_ftp_labels = os.path.join(ftplink, "{}.gz".format(fn_labels))

    pf_local_sequence = os.path.join(pd_gcfid, "sequence.fasta")
    pf_local_labels = os.path.join(pd_gcfid, "ncbi.gff")

    # don't re-download. TODO: add option to force re-download
    if os.path.isfile(pf_local_sequence) and os.path.isfile(pf_local_labels):
        return

    run_shell_cmd(
        "pwd; cd {}; wget --quiet {}; wget --quiet {}; gunzip -f {}; gunzip -f {}".format(
            pd_gcfid,
            pf_ftp_sequence,
            pf_ftp_labels,
            "{}.gz".format(fn_sequence),
            "{}.gz".format(fn_labels)
        ),

    )

    run_shell_cmd(
        "cd {}; mv {} {}; mv {} {}".format(
            pd_gcfid,
            fn_sequence, "sequence.fasta",
            fn_labels, "ncbi.gff"
        )
    )


def download_data_from_assembly_summary(df_assembly_summary, pd_output, **kwargs):
    # type: (pd.DataFrame, str, Dict[str, Any]) -> GenomeInfoList

    pf_output = get_value(kwargs, "pf_output", None)

    pd_output = os.path.abspath(pd_output)
    success_downloads = list()
    for _, gcfid_info in df_assembly_summary.iterrows():
        try:
            set_up_gcfid(gcfid_info, pd_output)
            success_downloads.append(gcfid_info)
        except (IOError, OSError):
            pass

    gil = GenomeInfoList([
        GenomeInfo("{}_{}".format(d["assembly_accession"], d["asm_name"]), 11) for d in success_downloads
    ])

    if pf_output is not None:
        gil.to_file(pf_output)

    return gil
