import copy
import os
import logging
from typing import *

import sbsp_io
from sbsp_alg.feature_computation import compute_features
from sbsp_alg.filtering import filter_orthologs
from sbsp_alg.msa import run_sbsp_msa, get_files_per_key, run_sbsp_msa_from_multiple, \
    run_sbsp_msa_from_multiple_for_multiple_queries
from sbsp_io.general import mkdir_p
from sbsp_general.general import get_value
from sbsp_alg.ortholog_finder import get_orthologs_from_files
from sbsp_alg.sbsp_compute_accuracy import pipeline_step_compute_accuracy, separate_msa_outputs_by_stats
from sbsp_general import Environment
from sbsp_io.general import read_rows_to_list
from sbsp_io.msa_2 import add_true_starts_to_msa_output
from sbsp_options.pbs import PBSOptions
from sbsp_options.pipeline_sbsp import PipelineSBSPOptions
from sbsp_parallelization.pbs import PBS
from sbsp_pbs_data.mergers import merge_identity
from sbsp_pbs_data.splitters import *

log = logging.getLogger(__name__)


def duplicate_pbs_options_with_updated_paths(env, pbs_options, **kwargs):
    # type: (Environment, PBSOptions, Dict[str, Any]) -> PBSOptions
    keep_on_head = get_value(kwargs, "keep_on_head", False, default_if_none=True)

    pbs_options = copy.deepcopy(pbs_options)
    pbs_options["pd-head"] = os.path.abspath(env["pd-work"])

    if keep_on_head:
        pbs_options["pd-root-compute"] = os.path.abspath(env["pd-work"])
    elif pbs_options["pd-root-compute"] is None:
        pbs_options["pd-root-compute"] = os.path.abspath(env["pd-work"])

    return pbs_options


def run_step_generic(env, pipeline_options, step_name, splitter, merger, data, func, func_kwargs, **kwargs):
    # type: (Environment, PipelineSBSPOptions, str, Callable, Callable, Dict[str, Any], Callable, Dict[str, Any], Dict[str, Any]) -> Dict[str, Any]

    output = {
        "pf-list-output": os.path.join(env["pd-work"], "pbs-summary.txt")
    }

    if pipeline_options.use_pbs():
        pbs_options = duplicate_pbs_options_with_updated_paths(env, pipeline_options["pbs-options"])

        if pbs_options.safe_get("pd-data-compute"):
            env = env.duplicate({"pd-data": pbs_options["pd-data-compute"]})

        pbs = PBS(env, pbs_options,
                  splitter=splitter,
                  merger=merger
                  )

        if pipeline_options.perform_step(step_name):

            output = pbs.run(
                data=data,
                func=func,
                func_kwargs=func_kwargs
            )
        else:
            # read data from file
            list_pf_output_packages = read_rows_to_list(os.path.join(env["pd-work"], "pbs-summary.txt"))
            output = pbs.merge_output_package_files(list_pf_output_packages)

    return output




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

        if pbs_options.safe_get("pd-data-compute"):
            env = env.duplicate({"pd-data": pbs_options["pd-data-compute"]})

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
                    "fn_q_labels": pipeline_options["fn-q-labels"], "fn_t_labels": pipeline_options["fn-t-labels"],
                    "max_evalue": pipeline_options.safe_get("max-evalue"),
                    "clean": True
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

        if pbs_options.safe_get("pd-data-compute"):
            env = env.duplicate({"pd-data": pbs_options["pd-data-compute"]})

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
                    "clean": True
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

        if pbs_options.safe_get("pd-data-compute"):
            env = env.duplicate({"pd-data": pbs_options["pd-data-compute"]})

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
                    "msa_options": pipeline_options["msa-options"],
                    "clean": True
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
        pbs_options = duplicate_pbs_options_with_updated_paths(env, pipeline_options["pbs-options"], keep_on_head=False)

        if pbs_options.safe_get("pd-data-compute"):
            env = env.duplicate({"pd-data": pbs_options["pd-data-compute"]})

        # pbs = PBS(env, pbs_options,
        #           splitter=split_list_and_remerge_by_key,
        #           merger=merge_identity
        #           )

        pbs = PBS(env, pbs_options,
                  splitter=split_q3prime_to_list_of_data_files,
                  merger=merge_identity
                  )

        if pipeline_options.perform_step("build-msa"):

            # get files per 3prime key
            q3prime_to_list_pf = get_files_per_key(list_pf_previous)
            pd_msa = os.path.join(pbs_options["pd-head"], "msa")
            mkdir_p(pd_msa)


            output = pbs.run(
                data={"q3prime_to_list_pf": q3prime_to_list_pf,
                      "pf_output_template": os.path.join(pbs_options["pd-head"],
                                                         pipeline_options["fn-msa"] + "_{}")},
                func=run_sbsp_msa_from_multiple_for_multiple_queries,
                func_kwargs={
                    "env": env,
                    "msa_options": pipeline_options["msa-options"],
                    "clean": True,
                    "pd_msa_final": pd_msa
                }
            )


            # output = pbs.run(
            #     data={"list_pf_data": list_pf_previous, "group_key": "q-3prime",
            #           "pf_output_template": os.path.join(pbs_options["pd-head"],
            #                                              pipeline_options["fn-msa"] + "_{}")},
            #     func=run_sbsp_msa,
            #     func_kwargs={
            #         "env": env,
            #         "msa_options": pipeline_options["msa-options"]
            #     }
            # )
        else:
            # read data from file
            list_pf_output_packages = read_rows_to_list(os.path.join(env["pd-work"], "pbs-summary.txt"))
            output = pbs.merge_output_package_files(list_pf_output_packages)

    return output


def sbsp_step_accuracy(env, pipeline_options, list_pf_previous):
    # type: (Environment, PipelineSBSPOptions, List[str]) -> List[str]
    """
    Given a list of query and target genomes, find the set of related genes
    for each query
    """

    log.debug("Running: sbsp_step_accuracy")

    mkdir_p(env["pd-work"])

    df = pd.concat([pd.read_csv(f, header=0) for f in list_pf_previous], ignore_index=True)

    df = pipeline_step_compute_accuracy(env, df, pipeline_options)

    df.to_csv(pipeline_options["pf-output"])


    # copy labels
    add_true_starts_to_msa_output(env, df, fn_q_labels_true=pipeline_options["fn-q-labels-true"])
    add_true_starts_to_msa_output(env, df, msa_nt=True, fn_q_labels_true=pipeline_options["fn-q-labels-true"])
    separate_msa_outputs_by_stats(env, df, pipeline_options["dn-msa-output"])

    return list()




