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

    output = {
        "pf-list-output": os.path.join(env["pd-work"], "output_summary.txt")
    }

    if pipeline_options.perform_step("get-orthologs"):

        if pipeline_options.use_pbs():
            pbs = PBS(env, pipeline_options["pbs-options"],
                      splitter=split_query_genomes_target_genomes_one_vs_group,
                      merger=merge_identity
            )

            pbs.run(
                data={"pf-q-list": pipeline_options["pf-q-list"], "pf-t-list": pipeline_options["pf-t-list"]},
                func=get_orthologs_from_files,
                func_kwargs=dict()
            )

    return output


def sbsp_step_compute_features(env, pipeline_options):
    # type: (Environment, PipelineSBSPOptions) -> Dict[str, Any]
    """
    Given a list of query and target genomes, find the set of related genes
    for each query
    """

    log.debug("Running: sbsp_step_compute_features")

    output = {
        "pf-list-output": os.path.join(env["pd-work"], "output_summary.txt")
    }

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




