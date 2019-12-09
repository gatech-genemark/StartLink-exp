import copy
import os
import logging
from typing import *

from sbsp_alg.feature_computation import compute_features
from sbsp_alg.filtering import filter_orthologs
from sbsp_alg.ortholog_finder import get_orthologs_from_files
from sbsp_general import Environment
from sbsp_io.general import read_rows_to_list
from sbsp_options.pbs import PBSOptions
from sbsp_options.pipeline_sbsp import PipelineSBSPOptions
from sbsp_parallelization.pbs import PBS
from sbsp_pbs_data.mergers import merge_identity
from sbsp_pbs_data.splitters import *

log = logging.getLogger(__name__)


def duplicate_pbs_options_with_updated_paths(env, pbs_options):
    # type: (Environment, PBSOptions) -> PBSOptions
    pbs_options = copy.deepcopy(pbs_options)
    pbs_options["pd-head"] = os.path.abspath(env["pd-work"])
    if pbs_options["pd-root-compute"] is None:
        pbs_options["pd-root-compute"] = os.path.abspath(env["pd-work"])

    return pbs_options

def sbsp_step_get_orthologs(env, pipeline_options):
    # type: (Environment, PipelineSBSPOptions) -> Dict[str, Any]
    """
    Given a list of query and target genomes, find the set of related genes
    for each query
    """

    log.debug("Running: sbsp_step_get_orthologs")

    output = {
        "pf-list-output": os.path.join(env["pd-work"], "pbs-summary.txt")
    }



    if pipeline_options.use_pbs():

        pbs_options = duplicate_pbs_options_with_updated_paths(env, pipeline_options["pbs-options"])

        pbs = PBS(env, pbs_options,
                  splitter=split_query_genomes_target_genomes_one_vs_group,
                  merger=merge_identity
        )

        if pipeline_options.perform_step("get-orthologs"):
            output = pbs.run(
                data={"pf_q_list": pipeline_options["pf-q-list"], "pf_t_list": pipeline_options["pf-t-list"],
                      "pf_output_template": os.path.join(pbs_options["pd-head"], pipeline_options["fn-orthologs"] + "_{}")},
                func=get_orthologs_from_files,
                func_kwargs={
                    "env": env,
                }
            )
        else:
            # read data from file
            list_pf_output_packages = read_rows_to_list(os.path.join(env["pd-work"], "pbs-summary.txt"))
            output = pbs.merge_output_package_files(list_pf_output_packages)


    return output


def sbsp_step_compute_features(env, pipeline_options, list_pf_previous):
    # type: (Environment, PipelineSBSPOptions, List[str]) -> Dict[str, Any]
    """
    Given a list of query and target genomes, find the set of related genes
    for each query
    """

    log.debug("Running: sbsp_step_compute_features")

    output = {
        "pf-list-output": os.path.join(env["pd-work"], "pbs-summary.txt")
    }

    if pipeline_options.use_pbs():

        pbs_options = duplicate_pbs_options_with_updated_paths(env, pipeline_options["pbs-options"])

        pbs = PBS(env, pbs_options,
                  splitter=split_list,
                  merger=merge_identity
                  )

        if pipeline_options.perform_step("compute-features"):

            output = pbs.run(
                data={"list_pf_data": list_pf_previous,
                      "pf_output_template": os.path.join(pbs_options["pd-head"],
                                                         pipeline_options["fn-compute-features"] + "_{}")},
                func=compute_features,
                func_kwargs={
                    "env": env,
                }
            )
        else:
            # read data from file
            list_pf_output_packages = read_rows_to_list(os.path.join(env["pd-work"], "pbs-summary.txt"))
            output = pbs.merge_output_package_files(list_pf_output_packages)

    return output

def sbsp_step_filter(env, pipeline_options, list_pf_previous):
    # type: (Environment, PipelineSBSPOptions, List[str]) -> Dict[str, Any]
    """
    Given a list of query and target genomes, find the set of related genes
    for each query
    """

    log.debug("Running: sbsp_step_filter")

    output = {
        "pf-list-output": os.path.join(env["pd-work"], "pbs-summary.txt")
    }

    if pipeline_options.use_pbs():
        pbs_options = duplicate_pbs_options_with_updated_paths(env, pipeline_options["pbs-options"])

        pbs = PBS(env, pbs_options,
                  splitter=split_list,
                  merger=merge_identity
                  )


        if pipeline_options.perform_step("filter"):

            output = pbs.run(
                data={"list_pf_data": list_pf_previous,
                      "pf_output_template": os.path.join(pbs_options["pd-head"],
                                                         pipeline_options["fn-filter"] + "_{}")},
                func=filter_orthologs,
                func_kwargs={
                    "env": env,
                    "msa_options": pipeline_options["msa-options"]
                }
            )
        else:
            # read data from file
            list_pf_output_packages = read_rows_to_list(os.path.join(env["pd-work"], "pbs-summary.txt"))
            output = pbs.merge_output_package_files(list_pf_output_packages)

    return output


def sbsp_step_msa(env, pipeline_options, list_pf_previous):
    # type: (Environment, PipelineSBSPOptions, List[str]) -> Dict[str, Any]
    """
    Given a list of query and target genomes, find the set of related genes
    for each query
    """

    log.debug("Running: sbsp_step_msa")

    output = {
        "pf-list-output": os.path.join(env["pd-work"], "pbs-summary.txt")
    }

    if pipeline_options.use_pbs():
        pbs_options = duplicate_pbs_options_with_updated_paths(env, pipeline_options["pbs-options"])

        pbs = PBS(env, pbs_options,
                  splitter=split_list_and_remerge_by_key,
                  merger=merge_identity
                  )

        if pipeline_options.perform_step("msa"):

            output = pbs.run(
                data={"list_pf_data": list_pf_previous, "group_key": "q-3prime",
                      "pf_output_template": os.path.join(pbs_options["pd-head"],
                                                         pipeline_options["fn-msa"] + "_{}")},
                func=filter_orthologs,
                func_kwargs={
                    "env": env,
                    "msa_options": pipeline_options["msa_options"]
                }
            )
        else:
            # read data from file
            list_pf_output_packages = read_rows_to_list(os.path.join(env["pd-work"], "pbs-summary.txt"))
            output = pbs.merge_output_package_files(list_pf_output_packages)

    return output


def sbsp_step_accuracy(env, pipeline_options, list_pf_previous):
    # type: (Environment, PipelineSBSPOptions, List[str]) -> Dict[str, Any]
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




