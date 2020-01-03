from __future__ import print_function
import os
import errno
import random
import string
import json
import sys
from typing import *


def generate_random_non_existing_filename(pd_work):
    # type: (str) -> str
    while True:
        fn_random = ''.join(random.choice(string.ascii_lowercase) for x in range(10))
        pf_random = os.path.join(pd_work, fn_random)

        if not os.path.exists(pf_random):
            return pf_random


def write_string_to_file(a_string, pf_out):

    with open(pf_out, "w") as f:
        f.write(a_string)


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def split_file_with_header(pf_in, split_tag, num_splits, delimiter=",", pd_work=None):
    # type: (str, str, int, str) -> list
    """
    Split a file (that has a header line) into multiple files, each with the same header

    :param pf_in: Input file
    :param split_tag: Tag associated with output files. Output file names are the input file name, followed by an
    underscore, then the split_tag, then an underscore, then the file number
    :param num_splits: Number of split files
    :param delimiter: CSV file delimiter
    :return:
    """

    import os
    import csv

    import math
    import sbsp_general.general

    total_num_lines = int(sbsp_general.general.run_shell_cmd("wc -l {}".format(pf_in)).strip().split()[0])

    # compute number of lines per split file
    if num_splits == 1:
        num_lines_per_file = total_num_lines
    else:
        num_lines_per_file = math.ceil(total_num_lines / num_splits)

    row_limit = num_lines_per_file

    fn_base = os.path.basename(pf_in)
    if pd_work is None:
        pd_work = os.path.dirname(pf_in)
    extension = os.path.splitext(pf_in)[1]

    list_pf = list()

    with open(pf_in, "r") as filehandler:
        reader = csv.reader(filehandler, delimiter=delimiter)
        curr_file_num = 1
        pf_curr = os.path.join(
            pd_work,
            "{}_{}_{}{}".format(fn_base, split_tag, curr_file_num, extension)
        )

        current_out_writer = csv.writer(open(pf_curr, 'w'), delimiter=delimiter)
        list_pf.append(pf_curr)

        current_limit = row_limit

        headers = next(reader)
        current_out_writer.writerow(headers)

        for i, row in enumerate(reader):
            if i + 1 > current_limit:
                curr_file_num += 1
                current_limit = row_limit * curr_file_num

                pf_curr = os.path.join(
                    pd_work,
                    "{}_{}_{}{}".format(fn_base, split_tag, curr_file_num, extension)
                )

                current_out_writer = csv.writer(open(pf_curr, 'w'), delimiter=delimiter)
                current_out_writer.writerow(headers)
                list_pf.append(pf_curr)

            current_out_writer.writerow(row)

        return list_pf


def split_file_with_header_and_group(pf_in, split_tag, num_splits, group_by, delimiter=",", pd_work=None):
    # type: (str, str, int, str, str, str) -> list[str]

    import pandas as pd

    df = pd.read_csv(pf_in, header=0, delimiter=delimiter)

    if group_by not in df.columns.values:
        raise ValueError("Unknown column to group by ({})".format(group_by))
    import math

    block_start = 0
    max_block_size = max(int(math.ceil(len(df) / float(num_splits))), 1)

    df['original-order'] = range(1, len(df) + 1)

    prev_group = None

    curr_group_size = 0

    df.sort_values(by=[group_by], inplace=True)

    file_num = 1
    fn_base = os.path.basename(pf_in)
    if pd_work is None:
        pd_work = os.path.dirname(pf_in)
    extension = os.path.splitext(pf_in)[1]

    list_pf = list()

    while True:

        block_end = min(block_start + max_block_size, len(df))

        # make sure all group is in, otherwise move down

        while df.iloc[block_end-1].loc[group_by] == df.iloc[block_end-2].loc[group_by] and block_end < len(df):
            block_end += 1

        df_curr = df.iloc[block_start:block_end-1]

        # write to file
        pf_curr = os.path.join(
            pd_work,
            "{}_{}_{}{}".format(fn_base, split_tag, file_num, extension)
        )

        df_curr.to_csv(pf_curr, index=False)
        list_pf.append(pf_curr)

        file_num += 1
        block_start = block_end

        if block_start >= len(df):
            break

    return list_pf


def read_rows_to_list(pf_data):
    # type: (str) -> List[str]

    with open(pf_data, "r") as f:

        mylist = list()

        for line in f:
            mylist.append(line.strip())

        return mylist



def merge_files_with_headers(list_pf, pf_out, warn_missing=False):
    # type: (list, str) -> None

    fout = open(pf_out, "w")

    if len(list_pf) == 0:
        fout.close()
        return

    # first file:
    idx_first_file = 0
    while idx_first_file < len(list_pf):

        try:
            f = open(list_pf[idx_first_file])

            for line in f:
                fout.write(line)

            idx_first_file += 1
            break
        except IOError:
            idx_first_file += 1

    # now the rest:
    for num in range(idx_first_file, len(list_pf)):

        try:
            f = open(list_pf[num])

            next(f)  # skip the header

            for line in f:
                fout.write(line)

            f.close()  # not really needed
        except Exception as e:
            if warn_missing:
                print("Missing file:", list_pf[num])
            else:
                raise e

    fout.close()


def read_genome_names(pf_list):
    # type: (str) -> list[str]
    names = list()

    with open(pf_list, "r") as f:

        for line in f:
            name = line.strip().split()[0]
            names.append(name)

    return names


def read_json(pf_json):
    # type: (str) -> Dict[str, Any]
    with open(pf_json, "r") as file:
        my_string = file.read()
        return json.loads(my_string)


def write_json(my_dict, pf_json):
    # type: (Dict[str, Any], str) -> None
    json.dump(my_dict, pf_json)


def print_progress(name, numer, denom=None):
    # type: (str, int, Union[int, None]) -> None

    sys.stdout.write("{} progress: {} / {} \r".format(name, numer, denom))

    sys.stdout.write("{} progress: {}".format(name, numer))
    if denom is not None:
        sys.stdout.write(" / {}".format(denom))
    sys.stdout.write("\r")

    sys.stdout.flush()