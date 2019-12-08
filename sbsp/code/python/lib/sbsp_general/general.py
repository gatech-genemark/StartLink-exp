import os
import re
import copy
import subprocess
import logging
from typing import *
from collections import namedtuple
import copy
import pandas as pd


ENV = []


def get_file_dir(file):
    return os.path.dirname(os.path.realpath(file))


def run_shell_cmd(cmd, do_not_log=False):
    # type: (str, bool) -> str
    if not do_not_log:
        logging.debug(cmd)

    return subprocess.check_output(cmd, shell=True).decode("utf-8")


def is_tool(name):
    # check if bash tool exists
    try:
        devnull = open(os.devnull)
        subprocess.Popen([name], stdout=devnull, stderr=devnull).communicate()
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            return False
    return True

def except_if_not_in_set(element, options):
    # type: (object, iter) -> None

    if options is None:
        raise ValueError("Options cannot be None")

    if element not in options:
        raise ValueError("Element ({}) not in options ({})", element, options)


def get_key_value_from_definition_line(key, defline):
    m = re.match(".*:" + str(key) + "=([^:]*)", defline)

    if m:
        return m.group(1)

    raise ValueError("Key " + str(key) + " not in definition line")

# FASTA header parsing
def get_genome_name_from_defition_line(defline):
    p = re.compile(r"^.*tag=([^;]+);(.*);(.*).*$")

    x = p.match(str(defline))
    if x:
        return x.group(1)

    raise ValueError("Definition line doesn't contain 'tag' key")

def get_fields_from_definition_line(defline):

    fields = {}

    # RE Groups
    # Accession
    #
    # >NC_000913:tag=Escherichia_coli_K_12_substr__MG1655_uid57779;11;:gc=52.41:pos=4179996;4180493;+:cds=4179996;4180493;+:type=prot:key=NC_000913;4180492;+

    m = re.search('>?([^:]+)', defline)
    if m:
        fields["accession"] = m.group(1)

    # get all keys in string
    keys = re.findall(r":([^=]+)=", defline)

    # for each key, get value
    for k in keys:
        fields[k] = get_key_value_from_definition_line(k, defline)

    return fields


def GetGenomeNameFromDefLine(defline):
    p = re.compile(r"^.*tag=([^;]+);(.+);(.*).*$")

    x = p.match(str(defline))
    if x:
        return x.group(1)

    raise ValueError("Definition line doesn't contain 'tag' key")

def expand_definition_line(definition_line, delimiter=":", key_prefix=""):
    KT = namedtuple("KeyTypePair", ["name", "type", "default"])


    def expand_field(field, component_setup, delimiter, key_prefix=""):
        # type: (str, dict, str) -> dict

        elements = field.split(delimiter)

        if len(elements) != len(component_setup):
            raise ValueError("Error in parsing field")

        output = dict()

        for e, c in zip(elements, component_setup):
            if len(e) > 0:
                output["{}{}".format(key_prefix, c.name)] = c.type(e)
            else:
                output["{}{}".format(key_prefix, c.name)] = c.default


        return output



    fields = get_fields_from_definition_line(definition_line)

    components = {
        "cds": [KT("left", int, ""), KT("right", int, ""), KT("strand", str, "")],
        "gc": [KT("gc", float, "")],
        "tag": [KT("genome", str, ""), KT("gcode", int, ""), KT("gc", float, "")],
        "type": [KT("type", str, "")]
    }

    output = dict()

    output["{}accession".format(key_prefix)] = fields["accession"]

    for key in fields.keys():
        if key in components.keys():
            output.update(expand_field(fields[key], components[key], ";", key_prefix))

    return output





def construct_key_to_index_dict(keys):

    key_to_index = {}
    index = 0

    for k in keys:
        key_to_index[k] = index
        index += 1

    return key_to_index


def merge_2d_dictionaries(dict1, dict2, func_merge_elements, in_place=False):
    # type: (dict, dict, callable, bool) -> dict

    # if in_place set to true, dict1 will be updated (and still returned)

    new_dict = dict1
    if not in_place:
        new_dict = copy.deepcopy(dict1)

    for a in dict2:

        if a not in new_dict:
            new_dict[a] = dict()

        for b in dict2[a]:

            element1 = None
            if b in new_dict[a]:
                element1 = new_dict[a][b]

            dict1[a][b] = func_merge_elements(element1, dict2[a][b])

    return new_dict


def get_value_DEPRECATED(key_value_pairs, key, default_value):

    return key_value_pairs[key] if key in key_value_pairs else default_value

def get_value(key_value_pairs, key, default_value=None, **kwargs):
    # type: (Dict[str, Any], str, Any, Dict[str, Any]) -> Any
    # TODO: implement

    invalid = get_value_DEPRECATED(kwargs, "invalid", None)
    perform_copy = get_value_DEPRECATED(kwargs, "copy", False)
    default_if_none = get_value_DEPRECATED(kwargs, "default_if_none", False)
    value_type = get_value_DEPRECATED(kwargs, "type", None)

    value = key_value_pairs[key] if key in key_value_pairs else default_value

    if default_if_none and value is None:
        return default_value

    if perform_copy:
        value = copy.deepcopy(value)

    if value_type:
        value = value_type(value)

    return value


def create_gene_key(genome=None, accession=None, left=None, right=None, strand=None, delimiter=";"):
    # type: (object, object, object, object, object, str) -> str

    return "{}{}{}{}{}{}{}{}{}".format(
        genome, delimiter,
        accession, delimiter,
        left, delimiter,
        right, delimiter,
        strand
    )

def create_3prime_key(genome=None, accession=None, left=None, right=None, strand=None, delimiter=";"):
    if strand == "+":
        return "{};{};{};{};{}".format(genome, accession, "",
                                       right, strand)
    else:
        return "{};{};{};{};{}".format(genome, accession, left,
                                       "", strand)


def df_add_gene_key(df, source, key_suffix, coordinates_suffix=None):
    # type: (pd.DataFrame, str, str) -> None

    except_if_not_in_set(source, ["q-", "t-"])

    gene_key = "{}{}".format(source, key_suffix)
    left_key = "{}left{}".format(source, coordinates_suffix)
    right_key = "{}right{}".format(source, coordinates_suffix)
    strand_key = "{}strand{}".format(source, coordinates_suffix)
    genome_key = "{}genome".format(source)
    accession_key = "{}accession".format(source)

    # create key for gene
    df[gene_key] = df.apply(
        lambda r: "{};{};{};{}".format(r[genome_key], r[accession_key], r[right_key], r[strand_key])
        if r[strand_key] == "+"
        else "{};{};{};{}".format(r[genome_key], r[accession_key], r[left_key], r[strand_key]), axis=1)

def df_add_5prime(df, source, key_suffix="5prime", coordinates_suffix=None):
    # type: (pd.DataFrame, str, str) -> None

    except_if_not_in_set(source, ["q-", "t-"])
    if coordinates_suffix is None:
        coordinates_suffix = ""

    result_key = "{}{}".format(source, key_suffix)
    left_key = "{}left{}".format(source, coordinates_suffix)
    right_key = "{}right{}".format(source, coordinates_suffix)
    strand_key = "{}strand{}".format(source, coordinates_suffix)

    # create key for gene
    df[result_key] = df.apply(
        lambda r: "{}".format(r[left_key]) if r[strand_key] == "+"
        else "{}".format(r[right_key]), axis=1)

def df_add_5prime_3prime_key(df, source, key_suffix="5prime-3prime", coordinates_suffix=None):
    # type: (pd.DataFrame, str, str) -> None

    except_if_not_in_set(source, ["q-", "t-"])
    if coordinates_suffix is None:
        coordinates_suffix = ""

    result_key = "{}{}".format(source, key_suffix)
    left_key = "{}left{}".format(source, coordinates_suffix)
    right_key = "{}right{}".format(source, coordinates_suffix)
    strand_key = "{}strand{}".format(source, coordinates_suffix)
    genome_key = "{}genome".format(source)
    accession_key = "{}accession".format(source)

    df[result_key] = df.apply(
        lambda r: "{};{};{};{};{}".format(r[genome_key], r[accession_key], r[left_key],
                                          r[right_key], r[strand_key]), axis=1)

def df_add_3prime_key(df, source, key_suffix="3prime", coordinates_suffix=None):
    # type: (pd.DataFrame, str, str) -> None

    except_if_not_in_set(source, ["q-", "t-"])
    if coordinates_suffix is None:
        coordinates_suffix = ""

    gene_key = "{}{}".format(source, key_suffix)
    left_key = "{}left{}".format(source, coordinates_suffix)
    right_key = "{}right{}".format(source, coordinates_suffix)
    strand_key = "{}strand{}".format(source, coordinates_suffix)
    genome_key = "{}genome".format(source)
    accession_key = "{}accession".format(source)

    # create key for gene
    df[gene_key] = df.apply(
        lambda r: "{};{};{};{}".format(r[genome_key], r[accession_key], r[right_key], r[strand_key])
        if r[strand_key] == "+"
        else "{};{};{};{}".format(r[genome_key], r[accession_key], r[left_key], r[strand_key]), axis=1)

def df_add_3prime_key_new(df, source, key_suffix="3prime", coordinates_suffix=None, use_genome=False):
    # type: (pd.DataFrame, str, str) -> None

    except_if_not_in_set(source, ["q-", "t-"])
    if coordinates_suffix is None:
        coordinates_suffix = ""

    gene_key = "{}{}".format(source, key_suffix)
    left_key = "{}left{}".format(source, coordinates_suffix)
    right_key = "{}right{}".format(source, coordinates_suffix)
    strand_key = "{}strand{}".format(source, coordinates_suffix)
    accession_key = "{}accession".format(source)

    genome_key = None
    if use_genome:
        genome_key = "{}genome".format(source)
    # create key for gene
    df[gene_key] = df.apply (
        lambda r: "{};{};{};{};{}".format(r[genome_key] if genome_key is not None else None, r[accession_key], "", r[right_key], r[strand_key])
        if r[strand_key] == "+"
        else "{};{};{};{};{}".format(r[genome_key] if genome_key is not None else None,
                                      r[accession_key], r[left_key], "", r[strand_key]), axis=1)




def add_gaps_to_nt_based_on_aa(seq_nt, seq_aa_with_gaps, preserve_case=True):
    # type: (str, str) -> str

    # make sure number of nt is 3 times number of amino acids (exclude gaps)
    num_aa = len(seq_aa_with_gaps) - seq_aa_with_gaps.count('-')
    # if len(seq_nt) != 3 * num_aa:
    #     raise ValueError("Number of nucleotides ({}) should be 3 times the number of amino acids ({})".format(
    #         len(seq_nt), num_aa))

    seq_nt_with_gaps = ""

    pos_in_nt = 0
    pos_in_aa_with_gaps = 0

    while pos_in_aa_with_gaps < len(seq_aa_with_gaps):

        curr_aa = seq_aa_with_gaps[pos_in_aa_with_gaps]

        # if gap
        if curr_aa == "-":
            seq_nt_with_gaps += "---"       # 3 nt gaps = 1 aa gap
        else:
            if preserve_case:
                if curr_aa.isupper():
                    seq_nt_with_gaps += seq_nt[pos_in_nt:pos_in_nt + 3].upper()
                else:
                    seq_nt_with_gaps += seq_nt[pos_in_nt:pos_in_nt + 3].lower()
            else:
                seq_nt_with_gaps += seq_nt[pos_in_nt:pos_in_nt+3]       # add next 3 nucleotides
            pos_in_nt += 3

        pos_in_aa_with_gaps += 1

    return seq_nt_with_gaps


def list_find_first(a_list, a_filter):
    # type: (List[Any], Callable) -> Any
    for x in a_list:
        if a_filter(x):
            return x
    return None
