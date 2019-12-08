import os
import logging
from typing import *

from sbsp_general import Environment
from sbsp_general.blast import run_blast, convert_blast_output_to_csv
from sbsp_general.general import get_value
from sbsp_io.general import mkdir_p

log = logging.getLogger(__name__)


def get_orthologs_from_files(env, pf_q_list, pf_t_list, pf_output, **kwargs):
    # type: (Environment, str, str, str, Dict[str, Any]) -> str

    clean = get_value(kwargs, "clean", False)

    # pf_q_list = data["pf-q-list"]
    # pf_t_list = data["pf-t-list"]

    pd_work = env["pd-work"]

    mkdir_p(pd_work)

    # run blast
    fn_blast_out = "blast.xml"
    pf_blast_out = os.path.join(pd_work, fn_blast_out)

    run_blast(env, pf_q_list, pf_t_list, pf_blast_out, **kwargs)

    # convert blast output to csv
    convert_blast_output_to_csv(pf_blast_out, pf_output, select_best_alignment_per_qt_pair=True)

    if clean:
        try:
            os.remove(pf_blast_out)
        except OSError:
            pass

    return pf_output



