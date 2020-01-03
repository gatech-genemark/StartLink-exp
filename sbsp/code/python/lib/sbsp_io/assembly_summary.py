import ftplib
import os
from ftplib import FTP
import logging

import pandas as pd
from typing import *

from sbsp_general.general import get_value
logger = logging.getLogger(__name__)


def read_assembly_summary(fname):
    # type: (str) -> dict

    # Output structure:
    # {
    #    "column_names": [column1, column2, ..., columnN],
    #    "data": [
    #        {"column1": value1, "column2": value2, ..., "columnN": valueN},
    #        {"column1": value1, "column2": value2, ..., "columnN": valueN},
    #        ...
    #    ]
    # }

    with open(fname, "r") as f:

        columns = list()
        data = list()

        for line in f:

            line = line.rstrip("\n")
            line = line.lstrip("\t")

            if len(line) == 0:
                continue

            if line.strip()[0] == "#":
                # check for header
                if "assembly_accession" in line:
                    columns = line[1:].strip().split("\t")      # skip hash and then split into column names
            else:
                # data line
                data_arr = line.split("\t")

                name_data_pairs = dict()

                for pair in zip(columns, data_arr):
                    name_data_pairs[pair[0]] = pair[1]

                    if pair[0] == "taxid":
                        name_data_pairs[pair[0]] = int(pair[1])

                data += [name_data_pairs]

        return {"column_names": columns, "data": data}


def read_assembly_summary_into_dataframe(pf_assembly_summary):
    # type: (str) -> pd.DataFrame
    info = read_assembly_summary(pf_assembly_summary)
    return pd.DataFrame(data=info["data"], columns=info["column_names"])


def write_assembly_summary(summary, fname):
    # type: (dict, str) -> None

    # Output structure:
    # {
    #    "column_names": [column1, column2, ..., columnN],
    #    "data": [
    #        {"column1": value1, "column2": value2, ..., "columnN": valueN},
    #        {"column1": value1, "column2": value2, ..., "columnN": valueN},
    #        ...
    #    ]
    # }

    with open (fname, "w") as f:

        column_names = summary["column_names"]
        data = summary["data"]

        # create header line
        out = "#\t" + "\t".join(column_names)

        f.write(out + "\n")

        for d in data:

            out = ""
            prefix = ""
            for name in column_names:

                out += prefix + str(d[name])
                prefix = "\t"

            f.write(out + "\n")


def get_rows_by_key_from_dataframe(df_assembly_summary, key="taxid", **kwargs):
    # type: (pd.DataFrame, str, Dict[str, Any]) -> Dict[int, List[Dict[str, Any]]]

    valid_assembly_levels = get_value(kwargs, "valid_assembly_levels", {"Complete Genome", "Scaffold", "Contig"},
                                      default_if_none=True)

    result = dict()

    for _, d in df_assembly_summary.iterrows():
        taxid = int(d[key])

        if d["assembly_level"] not in valid_assembly_levels:
            continue

        if taxid not in result.keys():
            result[taxid] = list()
        result[taxid].append(d)

    return result


def get_rows_by_key(pf_assembly_summary, key="taxid", **kwargs):
    # type: (str, str, Dict[str, Any]) -> Dict[int, List[Dict[str, Any]]]

    valid_assembly_levels = get_value(kwargs, "valid_assembly_levels", {"Complete Genome", "Scaffold", "Contig"},
                                      default_if_none=True)

    df_assembly_summary = read_assembly_summary_into_dataframe(pf_assembly_summary)

    return get_rows_by_key_from_dataframe(df_assembly_summary, key, **kwargs)


def ncbi_ftp_file_exists(ftp_gff, ftp_socket):
    # type: (str, FTP) -> bool

    path_on_remote = ftp_gff.split("ftp.ncbi.nlm.nih.gov")[1]
    try:
        if len(ftp_socket.nlst(path_on_remote)) > 0:
            return True
    except ftplib.error_temp:
        return False
    except BrokenPipeError:
        return False
    except EOFError:
        return False

    return False


def filter_entries_with_equal_taxid(df_assembly_summary, **kwargs):
    # type: (pd.DataFrame, Dict[str, Any]) -> pd.DataFrame

    possible_assembly_levels = {"Complete Genome", "Scaffold", "Contig"}

    valid_assembly_levels = get_value(kwargs, "valid_assembly_levels", possible_assembly_levels, default_if_none=True)
    favor_refseq_version = get_value(kwargs, "favor_refseq_version", False)
    favor_assembly_level_order = get_value(kwargs, "favor_assembly_level_order", False)
    number_per_taxid = get_value(kwargs, "number_per_taxid", None)

    ftp = FTP("ftp.ncbi.nlm.nih.gov")
    ftp.login()

    if len(df_assembly_summary) == 0:
        return pd.DataFrame()

    def select_from_list(local_list_info, n):
        # type: (pd.DataFrame, Union[None, int]) -> pd.DataFrame
        candidate_list = local_list_info
        # if n is not None and len(candidate_list) > n:
        #
        #     if len(local_list_info) <= n:
        #         candidate_list local_list_info
        #
        # # return local_list_info.iloc[0:n]

        if n is None:
            n = len(local_list_info)
        elif n > len(local_list_info):
            n = len(local_list_info)

        num_valid = 0
        final_df = pd.DataFrame(columns=local_list_info.columns)

        for index, row in local_list_info.iterrows():
            if num_valid == n:
                break
            should_add = False

            # if is refseq (guaranteed download)
            if "GCF" in row["assembly_accession"]:
                logger.debug("Refseq")
                should_add = True

            # if is genbank and has refseq
            elif row["gbrs_paired_asm"] != "na" and len(row["gbrs_paired_asm"]) > 0:
                should_add = True
                logger.debug("Genbank but has Refseq")
            # if only genbank but has GFF file
            else:

                gcf = row["assembly_accession"]
                acc = row["asm_name"].replace(" ", "_")

                gcfid = "{}_{}".format(gcf, acc)
                ftp_gff = os.path.join(row["ftp_path"],  "{}_genomic.gff.gz".format(gcfid))

                logger.debug("Genbank only, {}".format(ftp_gff))

                if ncbi_ftp_file_exists(ftp_gff, ftp):
                    should_add = True
                    logger.debug("Genbank only, with GFF")
                else:
                    logger.debug("Genbank only, No GFF - ignoring")

            if should_add:
                final_df = final_df.append(row, ignore_index=True)

        return final_df

    df_filtered = pd.DataFrame(columns=df_assembly_summary.columns)

    # for each group of entries of equal taxid
    for taxid, df_group in df_assembly_summary.groupby("taxid", as_index=False):

        # sort by latest first
        df_group["sort-by"] = pd.to_datetime(df_group["seq_rel_date"], format="%Y/%m/%d")
        df_group = df_group.sort_values("sort-by", ascending=False)
        df_group.drop("sort-by", inplace=True, axis=1)

        if favor_assembly_level_order:

            # go in order of assembly levels
            for assembly_level in valid_assembly_levels:

                df_filtered = df_filtered.append(select_from_list(
                    df_group[df_assembly_summary["assembly_level"] == assembly_level],
                    number_per_taxid - len(df_filtered)
                ))

                if len(df_filtered) == number_per_taxid:
                    break
        else:
            df_filtered = df_filtered.append(select_from_list(df_group, number_per_taxid))

    return df_filtered
