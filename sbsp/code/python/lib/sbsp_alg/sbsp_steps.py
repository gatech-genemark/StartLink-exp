import copy
import os
import logging
from typing import *

from sbsp_alg.ortholog_finder import get_orthologs_from_files
from sbsp_general import Environment
from sbsp_options.pipeline_sbsp import PipelineSBSPOptions
from sbsp_parallelization.pbs import PBS
from sbsp_pbs_data.mergers import merge_identity
from sbsp_pbs_data.splitters import split_query_genomes_target_genomes_one_vs_group

log = logging.getLogger(__name__)


def sbsp_step_get_orthologs(env, pipeline_options):
    # type: (Environment, PipelineSBSPOptions) -> Dict[str, Any]
    """
    Given a list of query and target genomes, find the set of related genes
    for each query
    """

    log.debug("Running: sbsp_step_get_orthologs")

    print(env["pd-work"])
    output = {
        "pf-list-output": os.path.join(env["pd-work"], "pbs-summary.txt")
    }

    if pipeline_options.perform_step("get-orthologs"):

        if pipeline_options.use_pbs():
            pbs_options = copy.deepcopy(pipeline_options["pbs-options"])
            pbs_options["pd-head"] = os.path.abspath(env["pd-work"])
            if pbs_options["pd-root-compute"] is None:
                pbs_options["pd-root-compute"] = os.path.abspath(env["pd-work"])

            pbs = PBS(env, pbs_options,
                      splitter=split_query_genomes_target_genomes_one_vs_group,
                      merger=merge_identity
            )

            output = pbs.run(
                data={"pf_q_list": pipeline_options["pf-q-list"], "pf_t_list": pipeline_options["pf-t-list"],
                      "pf_output_template": os.path.join(pbs_options["pd-head"], pipeline_options["fn-orthologs"] + "_{}")},
                func=get_orthologs_from_files,
                func_kwargs={
                    "env": env,
                }
            )

    return output


def sbsp_step_compute_features(env, pipeline_options, state):
    # type: (Environment, PipelineSBSPOptions, List[str]) -> Dict[str, Any]
    """
    Given a list of query and target genomes, find the set of related genes
    for each query
    """

    log.debug("Running: sbsp_step_compute_features")

    output = {
        "pf-list-output": os.path.join(env["pd-work"], "output_summary.txt")
    }

    print(state)
    # TODO

    return output

def sbsp_step_filter(env, pipeline_options):
    # type: (Environment, PipelineSBSPOptions) -> Dict[str, Any]
    """
    Given a list of query and target genomes, find the set of related genes
    for each query
    """

    log.debug("Running: sbsp_step_filter")

    output = {
        "pf-list-output": os.path.join(env["pd-work"], "output_summary.txt")
    }

    # TODO

    return output

def sbsp_step_msa(env, pipeline_options):
    # type: (Environment, PipelineSBSPOptions) -> Dict[str, Any]
    """
    Given a list of query and target genomes, find the set of related genes
    for each query
    """

    log.debug("Running: sbsp_step_msa")

    output = {
        "pf-list-output": os.path.join(env["pd-work"], "output_summary.txt")
    }

    # TODO

    return output

def sbsp_step_accuracy(env, pipeline_options):
    # type: (Environment, PipelineSBSPOptions) -> Dict[str, Any]
    """
    Given a list of query and target genomes, find the set of related genes
    for each query
    """

    log.debug("Running: sbsp_step_accuracy")

    output = {
        "pf-list-output": os.path.join(env["pd-work"], "output_summary.txt")
    }

    # TODO

    return output




