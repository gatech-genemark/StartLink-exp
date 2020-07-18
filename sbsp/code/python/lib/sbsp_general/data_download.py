import os
import logging
import shutil
import subprocess
from datetime import datetime

import pandas as pd
from typing import *

from tqdm import tqdm

from sbsp_container.genome_list import GenomeInfoList, GenomeInfo
from sbsp_container.taxonomy_tree import TaxonomyTree
from sbsp_general.general import get_value, run_shell_cmd
from sbsp_io.assembly_summary import get_rows_by_key, get_rows_by_key_from_dataframe, filter_entries_with_equal_taxid
from sbsp_io.general import mkdir_p, print_progress

logger = logging.getLogger(__name__)

def create_ftplink_from_gcf_acc(gcf, acc):
    # type: (str, str) -> str
    link = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all"
    gc, remaining = gcf.split("_")
    if len(remaining) < 9:
        raise ValueError("Wrong info {}".format(gcf))

    return os.path.join(link, gc, remaining[0:3], remaining[3:6], remaining[6:9], "{}_{}".format(gcf, acc))

def files_are_different(pf_1, pf_2):
    # type: (str, str) -> bool

    try:
        output = run_shell_cmd(
            "diff {} {}".format(pf_1, pf_2)
        )

        return len(output.strip()) != 0
    except Exception:
        return True

def download_assembly_summary_entry(entry, pd_output, **kwargs):
    # type: (Dict[str, Any], str, Dict[str, Any]) -> Dict[str, Any]

    force_download = get_value(kwargs, "force_download", None, valid={"all", "annotation_changed"})

    # build name
    gcf = entry["assembly_accession"]
    acc = entry["asm_name"].replace(" ", "_")

    output = {
            "assembly_accession": gcf,
            "asm_name": acc,
            "name": entry["name"],
            "parent_id": entry["parent_id"] if "parent_id" in entry else "",
            "genetic_code": entry["genetic_code"]
            }

    ftplink = entry["ftp_path"]

    # if genbank and has refseq, prefer refseq
    if "GCA" in gcf and entry["gbrs_paired_asm"] != "na" and len(entry["gbrs_paired_asm"]) > 0:
        gcf = entry["gbrs_paired_asm"]
        output["assembly_accession"] = gcf
        ftplink = create_ftplink_from_gcf_acc(gcf, acc)

    gcfid = "{}_{}".format(gcf, acc)
    pd_gcfid = os.path.join(pd_output, gcfid)
    pd_runs = os.path.join(pd_gcfid, "runs")

    try:

        mkdir_p(pd_gcfid)
        mkdir_p(pd_runs)


        fn_sequence = "{}_genomic.fna".format(gcfid)
        fn_labels = "{}_genomic.gff".format(gcfid)

        pf_ftp_sequence = os.path.join(ftplink, "{}.gz".format(fn_sequence))
        pf_ftp_labels = os.path.join(ftplink, "{}.gz".format(fn_labels))

        for not_allowed in {"#", "(", ")", ","}:
            if not_allowed in pf_ftp_sequence or not_allowed in pf_ftp_labels:
                raise ValueError("Invalid character in path")

        for not_allowed in {"#", "(", ")", "/", ":", ","}:
            if not_allowed in fn_sequence or not_allowed in fn_labels:
                raise ValueError("Invalid character in path")

        pf_local_sequence = os.path.join(pd_gcfid, "sequence.fasta")
        pf_local_labels = os.path.join(pd_gcfid, "ncbi.gff")

        # don't re-download. TODO: add option to force re-download
        if force_download != "any" and os.path.isfile(pf_local_sequence) and os.path.isfile(pf_local_labels):
            if force_download is None:
                return output

            if force_download == "annotation_changed":
                run_shell_cmd(
                    "cd {}; mkdir temporary; cd temporary; wget --quiet {}; gunzip -f {};".format(
                        pd_gcfid,
                        pf_ftp_labels,
                        "{}.gz".format(fn_labels)
                    )
                )

                update = files_are_different(
                    pf_1=os.path.join(pd_gcfid, "temporary", fn_labels),
                    pf_2=os.path.join(pd_gcfid, "ncbi.gff")
                )

                if update:
                    run_shell_cmd("cd {}; mv {} ../ncbi.gff".format(
                        os.path.join(pd_gcfid, "temporary"),
                        fn_labels
                    ))

                    # download sequence file again
                    run_shell_cmd(
                        "pwd; cd {}; wget --quiet {}; gunzip -f {};".format(
                            pd_gcfid,
                            pf_ftp_sequence,
                            "{}.gz".format(fn_sequence),
                        ),
                    )

                    run_shell_cmd(
                        "cd {}; mv {} {};".format(
                            pd_gcfid,
                            fn_sequence, "sequence.fasta",
                        )
                    )

                # cleanup
                run_shell_cmd(
                    "cd {}; rm -r temporary".format(
                        pd_gcfid
                    )
                )
            elif force_download == "no_download":
                return output
            else:       # FIXME: it's getting out of control. Create different lists: updated, all valid, etc...
                raise ValueError("nope")
        else:
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
    except (IOError, OSError, ValueError, subprocess.CalledProcessError):
        # cleanup failed attempt
        if os.path.exists(pd_gcfid) and os.path.isdir(pd_gcfid):
            shutil.rmtree(pd_gcfid)
        raise ValueError("Could not download data for genome: {}".format(gcfid)) from None

    return output

def compute_gc_from_file(pf_sequence):
    # type: (str) -> float
    return 0.0

    from sbsp_io.sequences import read_fasta_into_hash
    sequences = read_fasta_into_hash(pf_sequence)


    counts = {"A": 0, "C": 0, "G": 0, "T": 0}

    for seq in sequences.values():
        for s in seq:
            if s.upper() in {"A", "C", "G", "T"}:
                counts[s] += 1

    total = sum(counts.values())
    count_gc = counts["G"] + counts["C"]

    if total == 0:
        return 0.0

    return round(100 * count_gc / float(total), 2)

def count_cds(pf_labels):
    # type: (str) -> int

    return int(run_shell_cmd("grep -c CDS {}".format(pf_labels), do_not_log=True))

def get_annotation_date(pf_labels):
    return run_shell_cmd("grep -m 1 \"annotation-date\" {}".format(pf_labels) +
                         r" | awk '{print $2}'", do_not_log=True).strip()


def get_genome_specific_attributes(pd_data, info):
    # type: (str, Dict[str, Any]) -> Dict[str, Any]
    gcfid = "{}_{}".format(info["assembly_accession"], info["asm_name"])
    logger.debug("Genome specific: {}".format(gcfid))

    pd_gcfid = os.path.join(pd_data, gcfid)
    pf_sequences = os.path.join(pd_gcfid, "sequence.fasta")
    pf_labels = os.path.join(pd_gcfid, "ncbi.gff")

    gc = compute_gc_from_file(pf_sequences)
    try:
        num_genes = count_cds(pf_labels)
    except Exception:
        num_genes = 0
    annotation_date = get_annotation_date(pf_labels)

    return {"gc": gc, "num_genes": num_genes, "annotation_date": annotation_date}




def download_data_from_assembly_summary(df_assembly_summary, pd_output, **kwargs):
    # type: (pd.DataFrame, str, Dict[str, Any]) -> GenomeInfoList
    """
    Attempt to download all genomes from assembly summary.
    :param df_assembly_summary: Data frame containing assembly summary entries
    :param pd_output: Path to download directory
    :param kwargs:
        - pf_output: path to output file which will contain list of downloaded genomes
    :return: Genome information list of successfully downloaded entries
    """

    pf_output_list = get_value(kwargs, "pf_output_list", None)
    attributes = get_value(kwargs, "attributes", dict(), default_if_none=True)


    df_assembly_summary = filter_entries_with_equal_taxid(
        df_assembly_summary, **kwargs
    )

    pd_output = os.path.abspath(pd_output)
    success_downloads = list()
    total = 0
    for _, gcfid_info in tqdm(df_assembly_summary.iterrows(), "Downloading", total=len(df_assembly_summary)):
        total += 1
        logger.debug("Trying {}".format(gcfid_info["assembly_accession"]))

        try:
            gcfid_info = download_assembly_summary_entry(gcfid_info, pd_output, **kwargs)
            success_downloads.append(gcfid_info)

            # print_progress("Download", len(success_downloads), total)
        except (IOError, OSError, ValueError):
            # print_progress("Download", len(success_downloads), total)
            pass

    gil = GenomeInfoList([
        GenomeInfo(
            "{}_{}".format(d["assembly_accession"], d["asm_name"]),
            d["genetic_code"],
            attributes={
                "name": d["name"],
                "parent_id": d["parent_id"],
                **get_genome_specific_attributes(pd_output, d),
                **attributes
            }
        ) for d in success_downloads
    ])

    if pf_output_list is not None:
        gil.to_file(pf_output_list)

    return gil


def filter_assembly_summary_by_ancestor(ancestor_tag, tag_type, taxonomy_tree, df_assembly_summary):
    # type: (str, str, TaxonomyTree, pd.DataFrame) -> pd.DataFrame

    taxid_to_list_of_rows = get_rows_by_key_from_dataframe(df_assembly_summary, key="taxid")

    list_rows = list()
    df_filtered = pd.DataFrame(columns=df_assembly_summary.columns)

    for genome_node in tqdm(taxonomy_tree.get_possible_genomes_under_ancestor(ancestor_tag, tag_type), "Searching for nodes under ancestor"):

        # find rows for taxid in assembly summary
        tax_id = genome_node["taxid"]

        # append rows to dataframe
        if tax_id in taxid_to_list_of_rows:
            info_list = taxid_to_list_of_rows[tax_id]

            # add name to all genomes
            for i in range(len(info_list)):
                info_list[i]["name"] = genome_node["name_txt"].replace(",", " ")
                info_list[i]["parent_id"] = genome_node["parent_id"]
                info_list[i]["genetic_code"] = genome_node["genetic_code"]
            list_rows += info_list

    if len(list_rows) > 0:
        df_filtered = df_filtered.append(list_rows, sort=False)

    return df_filtered


def download_data_by_ancestor(ancestor_tag, tag_type, taxonomy_tree, df_assembly_summary, pd_output, **kwargs):
    # type: (str, str, TaxonomyTree, pd.DataFrame, str, Dict[str, Any]) -> GenomeInfoList

    # get assembly summary entries for genomes under ancestor
    df_assembly_summary_filtered = filter_assembly_summary_by_ancestor(
        ancestor_tag, tag_type, taxonomy_tree, df_assembly_summary
    )

    return download_data_from_assembly_summary(df_assembly_summary_filtered, pd_output,
                                               attributes={"ancestor": ancestor_tag}, **kwargs)
