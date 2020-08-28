from __future__ import print_function
import os
import errno
import random
import re
import string
import json
import sys
from typing import *

from Bio import SeqIO

from sbsp_general.general import os_join


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


def remove_p(*path):
    # type: (List[str]) -> None
    for p in path:
        if os.path.isfile(p):
            os.remove(p)



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


def read_fasta_into_hash(fnsequences, stop_at_first_space=True):
    # type: (str) -> dict[str, Bio.Seq.Seq]
    # read all sequences
    sequences = {}  # will contain all protein sequences from file
    try:
        for record in SeqIO.parse(fnsequences, "fasta"):
            if stop_at_first_space:
                sequences[record.name.strip().split()[0]] = record.seq
            else:
                sequences[record.description.strip()] = record.seq

    except Exception:
        pass

    return sequences


def create_attribute_dict(attribute_string, delimiter=";"):
    # type: (str, str) -> Dict[str, Any]

    attributes = dict()

    for current in attribute_string.strip().split(sep=delimiter):
        if len(current.strip().split(maxsplit=1)) == 2:
            k, v = current.strip().split(maxsplit=1)
            attributes[k] = v

    return attributes


def read_gff(pf_labels, shift=-1):
    # type: (str, int) -> Dict[str, List[Dict[str, Any]]]


    labels = dict() # type: Dict[str, List[Dict[str, Any]]]

    pattern = re.compile(r"([^\t]+)\t([^\t]+)\t(CDS)\t(\d+)\t(\d+)\t([^\t]+)\t([+-])\t([^\t]+)\t([^\t]+)")

    with open(pf_labels, "r") as f:

        for line in f:

            line = line.strip()

            m = pattern.match(line)
            if m:

                attributes = create_attribute_dict(m.group(9))
                seqname = m.group(1).strip()

                label = {
                    "left" : int(m.group(4)) + shift,
                    "right" : int(m.group(5)) + shift,
                    "strand" : m.group(7),
                    "seqname" : m.group(1),
                    "attributes": attributes
                }

                if seqname not in labels.keys():
                    labels[seqname] = list()

                labels[seqname].append(label)
    return labels

def read_lst(pf_labels, shift=-1):
    # type: (str, int) -> Dict[str, List[Dict[str, Any]]]


    labels = dict() # type: Dict[str, List[Dict[str, Any]]]

    # pattern = re.compile(r"([^\t]+)\t([^\t]+)\t(CDS)\t(\d+)\t(\d+)\t([^\t]+)\t([+-])\t([^\t]+)\t([^\t]+)")
    pattern = re.compile(r"([^\s]+)\s+([+-])\s+(\d+)\s+(\d+)\s+(\d+)\s+(.+)$")

    # out = str(counter)
    # out += " " + str(l["strand"])
    # out += " " + str(l["left"] + shift)
    # out += " " + str(l["right"] + shift)
    # out += " " + str(l["right"] - l["left"] + 1)
    # out += " " "nativebac" + " AGGAGG 6 1"
    # out += " " + l["attributes"]["gene_type"]

    seqname = None

    with open(pf_labels, "r") as f:

        for line in f:

            line = line.strip()

            if line.startswith("SequenceID:"):
                seqname = line.split(":", maxsplit=1)[1].strip()
                continue
            elif len(line.strip()) == 0 or seqname is None:
                continue

            m = pattern.match(line)
            if m:

                attributes = m.group(6)

                label = {
                    "left" : int(m.group(3)) + shift,
                    "right" : int(m.group(4)) + shift,
                    "strand" : m.group(2),
                    "seqname" : seqname,
                    "attributes": attributes
                }

                if seqname not in labels.keys():
                    labels[seqname] = list()

                labels[seqname].append(label)
    return labels

def convert_multi_fasta_into_single_fasta(sequences, labels, new_seqname):
    # type: (Dict[str, str], Dict[str, List[Dict[str, Any]]], str) -> [Dict[str, str], Dict[str, List[Dict[str, Any]]]]

    sorted_seqnames = sorted(sequences.keys())

    joined_seq = {new_seqname: ""}
    joined_labels = {new_seqname: list()}

    curr_offset = 0

    for seqname in sorted_seqnames:

        if seqname not in labels.keys():
            continue

        seqname_sequence = sequences[seqname]
        seqname_labels = labels[seqname]

        for l in seqname_labels:
            new_l = l.copy()
            new_l["left"] += curr_offset
            new_l["right"] += curr_offset
            new_l["seqname"] = new_seqname

            joined_labels[new_seqname].append(new_l)


        joined_seq[new_seqname] += seqname_sequence
        curr_offset += len(seqname_sequence)

    return [joined_seq, joined_labels]


def write_lst(labels, pf_output, shift=+1):
    # type: (Dict[str, List[Dict[str, Any]]], str, int) -> None

    # seqname, source, feature, left, right, score, strand, frame, attributes
    with open(pf_output, "w") as f:
        f.write("# GeneMark.hmm-2 LST format\n"
"# GeneMark.hmm-2 prokaryotic version: 1.14\n"
"# File with sequence: tmpseq.fna\n"
"# File with native parameters: itr_1.mod\n"
"# Native species name and build: gms2-training\n"
"# File with MetaGeneMark parameters: /storage4/karl/sbsp/biogem/sbsp/bin_external/gms2/mgm_11.mod\n"
"# translation table: 11\n"
"# output date start: Mon Jun  8 09:26:44 2020\n\n")
        for seqname, seqname_labels in labels.items():
            f.write(f"SequenceID: {seqname}\n")
            counter = 1

            for counter, l in enumerate(seqname_labels):

                out = str(counter)
                out += " " + str(l["strand"])
                out += " " + str(l["left"] + shift)
                out += " " + str(l["right"] + shift)
                out += " " + str(l["right"] - l["left"] + 1)
                out += " " "nativebac" + " AGGAGG 6 1"
                out += " " + l["attributes"]["gene_type"] if "gene-type" in l["attributes"] else " ."


                f.write(out + "\n")


def write_fasta_hash_to_file(fasta, pf_output):
    # type: (Dict[str, Any], str) -> None

    output = ""
    for header in fasta.keys():
        output += ">{}\n{}\n".format(
            header, fasta[header]
        )

    write_string_to_file(output, pf_output)



def convert_multi_fasta_to_single(env, pf_sequences, pf_labels):
    # pf_sequence = sys.argv[1]
    # pf_labels = sys.argv[2]
    # pd_work = sys.argv[3]

    org_seq = read_fasta_into_hash(pf_sequences)
    org_labels = read_gff(pf_labels, shift=0)

    new_seq, new_labels = convert_multi_fasta_into_single_fasta(org_seq, org_labels, "anydef")
    pd_work = env["pd-work"]

    pf_new_seq = os_join(pd_work, "sequence_joined")
    pf_new_labels = os_join(pd_work, "labels_joined_lst")
    import os
    # write_gff(new_labels, os.path.join(pd_work, "labels_joined"), shift=0)
    write_lst(new_labels, pf_new_labels, shift=0)
    write_fasta_hash_to_file(new_seq, pf_new_seq)

    return pf_new_seq, pf_new_labels