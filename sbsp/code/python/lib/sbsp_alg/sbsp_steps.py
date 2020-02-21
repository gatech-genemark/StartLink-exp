import copy
import os
import logging
from typing import *

from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SubsMat import MatrixInfo as matlist


import sbsp_io
from sbsp_alg.feature_computation import compute_features
from sbsp_alg.filtering import filter_orthologs
from sbsp_alg.msa import run_sbsp_msa, get_files_per_key, run_sbsp_msa_from_multiple, \
    run_sbsp_msa_from_multiple_for_multiple_queries, perform_msa_on_df, move_files_using_scp
from sbsp_general.blast import run_blast
from sbsp_io.general import mkdir_p
from sbsp_general.general import get_value
from sbsp_alg.ortholog_finder import get_orthologs_from_files, extract_labeled_sequences_for_genomes, \
    unpack_fasta_header, select_representative_hsp, create_info_for_query_target_pair, \
    compute_distance_based_on_local_alignment, compute_distance_based_on_global_alignment_from_sequences, \
    run_blast_on_sequences
from sbsp_alg.sbsp_compute_accuracy import pipeline_step_compute_accuracy, separate_msa_outputs_by_stats
from sbsp_general import Environment
from sbsp_io.general import read_rows_to_list
from sbsp_io.msa_2 import add_true_starts_to_msa_output
from sbsp_io.sequences import read_fasta_into_hash, write_fasta_hash_to_file
from sbsp_options.msa import MSAOptions
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
                data={"pf_q_list": pipeline_options["pf-q-list"],
                      "pf_output_template": os.path.join(pbs_options["pd-head"], pipeline_options["fn-orthologs"] + "_{}")},
                func=get_orthologs_from_files,
                func_kwargs={
                    "env": env,
                    "pf_t_db": pipeline_options["pf-t-db"],
                    "fn_q_labels": pipeline_options["fn-q-labels"], "fn_t_labels": pipeline_options["fn-t-labels"],
                    "max_evalue": pipeline_options.safe_get("max-evalue"),
                    "clean": True,
                    "sbsp_options": pipeline_options.safe_get("msa-options")
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


def run_sbsp_steps(env, data, pf_t_db, pf_output, msa_options, **kwargs):
    # type: (Environment, Dict[str, Seq], str, str, MSAOptions, Dict[str, Any]) -> str

    sequences = data
    tag = get_value(kwargs, "tag", None)
    hsp_criteria = get_value(kwargs, "hsp_criteria", None)
    msa_output_start = get_value(kwargs, "msa_output_start", 0)
    pd_msa_final = get_value(kwargs, "pd_msa_final", env["pd-work"])

    distance_min = 0.001
    distance_max = 0.5


    # write sequences to a file

    pf_q_sequences = os.path.join(env["pd-work"], "query_sequences_{}.fasta".format(msa_output_start))
    write_fasta_hash_to_file(sequences, pf_q_sequences)

    # run blast
    pf_blast_output = os.path.join(env["pd-work"], "blast_output.xml")
    run_blast_on_sequences(env, pf_q_sequences, pf_t_db, pf_blast_output, sbsp_options=msa_options)

    if not os.path.isfile(pf_blast_output):
        return pf_output

    # open blast stream
    try:
        f_blast_results = open(pf_blast_output, "r")
    except OSError as e:
        log.warning("Could not open blast results file: {}".format(pf_blast_output))
        raise e

    # start clean
    if os.path.isfile(pf_output):
        os.remove(pf_output)

    matrix = matlist.blosum62
    import sbsp_alg.phylogeny
    sbsp_alg.phylogeny.add_stop_codon_to_blosum(matrix)

    records = NCBIXML.parse(f_blast_results)
    header_written = False

    # for each blast query
    for r in records:

        query_info = unpack_fasta_header(r.query)

        list_entries = list()

        # for each alignment to a target protein for the current query
        for alignment in r.alignments:

            hsp = select_representative_hsp(alignment, hsp_criteria)

            target_info = unpack_fasta_header(alignment.title)

            # get nucleotide sequence that corresponds to protein
            original_q_nt = Seq(query_info["lorf_nt"][query_info["offset"]:])
            original_t_nt = Seq(target_info["lorf_nt"][target_info["offset"]:])

            distance = compute_distance_based_on_local_alignment(query_info, target_info, hsp,
                                                                 original_q_nt=original_q_nt,
                                                                 original_t_nt=original_t_nt,
                                                                 **kwargs)

            original_q_aa = original_q_nt.translate()
            original_t_aa = original_t_nt.translate()

            #global_distance, global_length, global_length_without_gaps = compute_distance_based_on_global_alignment_from_sequences(
            #    original_q_aa, original_t_aa, original_q_nt, original_t_nt, matrix
            #)
            global_distance = global_length = global_length_without_gaps = 0

            # FIXME: thresholds should be from input configuration files
            if distance > distance_min and distance < distance_max:
            #if True:

                output_info = create_info_for_query_target_pair(
                    query_info, target_info, hsp,
                    distance_blast=distance,
                    distance=distance,
                    global_distance=global_distance,
                    global_length=global_length,
                    global_length_without_gaps=global_length_without_gaps,
                    local_distance=distance,
                    local_length=hsp.align_length,
                    local_length_without_gaps=sum([
                        1 for i in range(len(hsp.query)) if hsp.query[i] != "-" and hsp.sbjct[i] != "-"
                    ])
                )

                # clean up
                output_info["q-prot-pos-5prime-in-frag-msa"] = query_info["offset"] / 3
                output_info["q-nucl-pos-5prime-in-frag-msa"] = query_info["offset"]
                output_info["q-prot-position-of-5prime-in-msa-fragment-no-gaps"] = query_info["offset"] / 3
                output_info["q-nucl-position-of-5prime-in-msa-fragment-no-gaps"] = query_info["offset"]
                output_info["t-prot-pos-5prime-in-frag-msa"] = target_info["offset"] / 3
                output_info["t-nucl-pos-5prime-in-frag-msa"] = target_info["offset"]
                output_info["t-prot-position-of-5prime-in-msa-fragment-no-gaps"] = target_info["offset"] / 3
                output_info["t-nucl-position-of-5prime-in-msa-fragment-no-gaps"] = target_info["offset"]

                output_info["q-prot-msa"] = Seq(query_info["lorf_nt"]).translate()._data
                output_info["t-prot-msa"] = Seq(target_info["lorf_nt"]).translate()._data

                output_info["q-nucl-msa"] = Seq(query_info["lorf_nt"])._data
                output_info["t-nucl-msa"] = Seq(target_info["lorf_nt"])._data

                list_entries.append(output_info)

        # run MSA on remaining targets

        if len(list_entries) == 0:
            continue

        print("{};{};{};{}".format(query_info["accession"], query_info["left"], query_info["right"], query_info["strand"]))

        df_entries = pd.DataFrame(list_entries)
        df_results = perform_msa_on_df(env, df_entries, msa_options=msa_options, msa_output_start=msa_output_start)

        # for each query in blast
        if pd_msa_final is not None:
            try:
                move_files_using_scp(df_results, pd_msa_final)
            except Exception:
                pass

        # write/append result to file
        if df_results is not None and len(df_results) > 0:
            if not os.path.isfile(pf_output):
                df_results.to_csv(pf_output, index=False)
            else:
                df_results.to_csv(pf_output, mode="a", index=False, header=False)

    return pf_output


def sbsp_steps(env, pipeline_options):
    # type: (Environment, PipelineSBSPOptions, List[str]) -> Dict[str, Any]
    """
    Given a list of query and target genomes, find the set of related genes
    for each query
    """

    log.debug("Running: sbsp steps")

    output = {
        "pf-list-output": os.path.join(env["pd-work"], "pbs-summary.txt")
    }

    # read input sequences
    q_gil = GenomeInfoList.init_from_file(pipeline_options["pf-q-list"])

    pf_aa = os.path.join(env["pd-work"], "query.faa")
    extract_labeled_sequences_for_genomes(env, q_gil, pf_aa, ignore_frameshifted=True, reverse_complement=True, ignore_partial=True )
    q_sequences = read_fasta_into_hash(pf_aa, stop_at_first_space=False)

    if pipeline_options.use_pbs():
        pbs_options = duplicate_pbs_options_with_updated_paths(env, pipeline_options["pbs-options"], keep_on_head=False)

        if pbs_options.safe_get("pd-data-compute"):
            env = env.duplicate({"pd-data": pbs_options["pd-data-compute"]})

        pbs = PBS(env, pbs_options,
                  splitter=split_dict,
                  merger=merge_identity
                  )

        if pipeline_options.perform_step("build-msa"):

            pd_msa = os.path.join(pbs_options["pd-head"], "msa")
            mkdir_p(pd_msa)


            output = pbs.run(
                data={"dict": q_sequences,
                      "pf_output_template": os.path.join(pbs_options["pd-head"],
                                                         pipeline_options["fn-msa"] + "_{}")},
                func=run_sbsp_steps,
                func_kwargs={
                    "env": env,
                    "pf_t_db": pipeline_options["pf-t-db"],
                    "msa_options": pipeline_options["msa-options"],
                    "clean": True,
                    "pd_msa_final": pd_msa
                }
            )

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




