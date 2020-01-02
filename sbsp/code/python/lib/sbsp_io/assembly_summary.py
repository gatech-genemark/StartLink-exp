import pandas as pd
from typing import *

from sbsp_general.general import get_value


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


def filter_entries_with_equal_taxid(df_assembly_summary, **kwargs):
    # type: (pd.DataFrame, Dict[str, Any]) -> pd.DataFrame

    possible_assembly_levels = {"Complete Genome", "Scaffold", "Contig"}

    valid_assembly_levels = get_value(kwargs, "valid_assembly_levels", possible_assembly_levels, default_if_none=True)
    favor_assembly_level_order = get_value(kwargs, "favor_assembly_level_order", False)
    number_per_taxid = get_value(kwargs, "number_per_taxid", None)

    if len(df_assembly_summary) == 0:
        return pd.DataFrame()

    def select_from_list(local_list_info, n):
        # type: (pd.DataFrame, Union[None, int]) -> pd.DataFrame

        if n is None:
            return local_list_info

        if len(local_list_info) <= n:
            return local_list_info

        return local_list_info.iloc[0:n]

    df_filtered = pd.DataFrame(columns=df_assembly_summary.columns)

    # for each group of entries of equal taxid
    for taxid, df_group in df_assembly_summary.groupby("taxid", as_index=False):

        # sort by latest first
        df_group["sort-by"] = pd.to_datetime(df_group["seq_rel_date"], format="%Y/%m/%d")
        df_group = df_group.sort_values("sort-by", ascending=False)
        df_group.drop("sort-by", inplace=True)

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
