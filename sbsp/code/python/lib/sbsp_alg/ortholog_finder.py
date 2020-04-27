import os
import logging
import re

import numpy as np
from typing import *

from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast.Record import Alignment, HSP
from Bio.Seq import Seq
from Bio.SubsMat import MatrixInfo as matlist

from sbsp_alg.feature_computation import add_gaps_to_nt_based_on_aa
from sbsp_alg.phylogeny import k2p_distance, global_alignment_aa_with_gap
from sbsp_container.genome_list import GenomeInfoList, GenomeInfo
from sbsp_general import Environment
from sbsp_general.blast import run_blast, convert_blast_output_to_csv, create_blast_database, run_blast_alignment
from sbsp_general.general import get_value
from sbsp_general.labels import Labels, Label, create_gene_key_from_label
from sbsp_io.general import mkdir_p
from sbsp_io.labels import read_labels_from_file
from sbsp_io.sequences import read_fasta_into_hash
from sbsp_options.sbsp import SBSPOptions

log = logging.getLogger(__name__)

valid_starts_pos = ["ATG", "GTG", "TTG"]
valid_starts_neg = ["CAT", "CAC", "CAA"]

valid_stops_pos = ["TAA", "TGA", "TAG"]
valid_stops_neg = ["TTA", "TCA", "CTA"]


def is_valid_start(codon, strand):
    if strand == "+":
        return codon in valid_starts_pos
    else:
        return codon in valid_starts_neg


def is_valid_stop(codon, strand):
    if strand == "+":
        return codon in valid_stops_pos
    else:
        return codon in valid_stops_neg


def get_orthologs_from_files_deprecated(env, pf_q_list, pf_t_list, pf_output, **kwargs):
    # type: (Environment, str, str, str, Dict[str, Any]) -> str

    clean = get_value(kwargs, "clean", False)

    # pf_q_list = data["pf-q-list"]
    # pf_t_list = data["pf-t-list"]

    pd_work = env["pd-work"]

    mkdir_p(pd_work)

    # run blast
    fn_blast_out = "blast.xml"
    pf_blast_out = os.path.join(pd_work, fn_blast_out)

    run_blast(env, pf_q_list, pf_t_list, pf_blast_out, **kwargs)

    # convert blast output to csv
    convert_blast_output_to_csv(pf_blast_out, pf_output, select_best_alignment_per_qt_pair=True)

    if clean:
        try:
            os.remove(pf_blast_out)
        except OSError:
            pass

    return pf_output


def pack_fasta_header(label, gi, **kwargs):
    # type: (Label, GenomeInfo, Dict[str, Any]) -> str

    def create_tags_from_dict(**local_kwargs):
        # type: (Dict[str, Any]) -> str
        out = ""
        for k, v in local_kwargs.items():
            out += ";{}={}".format(k, v)
        return out

    return "{} accession={};genome={};left={};right={};strand={}".format(
        label.seqname(),
        label.seqname(),
        gi.name,
        label.left() + 1,
        label.right() + 1,
        label.strand(),
    ) + create_tags_from_dict(**kwargs)


def unpack_fasta_header(header):
    # type: (str) -> Dict[str, Any]

    fields = {}

    def get_key_value_from_definition_line(key, l_defline):
        m = re.match(".*(?:^|;)" + str(key) + "=([^;]*)", l_defline)

        if m:
            return m.group(1)

        raise ValueError("Key " + str(key) + " not in definition line")

    # remove first accession
    header = header.strip().split(maxsplit=1)[1]

    # get all keys in string
    keys = re.findall(r"([^;=]+)=", header)

    type_mapper = {
        "left": int,
        "right": int,
        "gc": float,
        "offset": int,
        "upstream_left": int,
        "upstream_right": int
    }

    # for each key, get value
    for k in keys:
        fields[k] = get_key_value_from_definition_line(k, header)

        if k in type_mapper.keys():
            fields[k] = type_mapper[k](fields[k])

    return fields


def select_representative_hsp(alignment, hsp_criteria):
    # type: (Alignment, str) -> Union[HSP, None]

    selected_hsp = None
    max_length = None
    for hsp in alignment.hsps:
        if max_length is None:
            selected_hsp = hsp
            max_length = hsp.align_length
        else:
            if hsp.align_length > max_length:
                selected_hsp = hsp
                max_length = hsp.align_length

    return selected_hsp


def map_aligned_aa_to_aligned_nt(q_aligned_seq_aa, original_q_nt, q_start_aa, q_end_aa, offset_nt=0):
    # type: (Seq, Seq, int, int) -> Seq
    if len(original_q_nt) < (q_end_aa - q_start_aa + 1) * 3:
        raise ValueError("Nucleotide sequence length less than aligned amino acid fragment: {} < {}".format(
            len(original_q_nt), (q_end_aa - q_start_aa + 1) * 3
        ))

    output = Seq("")
    pos_nt_no_gaps = offset_nt + q_start_aa * 3
    max_pos_nt_no_gaps = offset_nt + (q_end_aa + 1) * 3

    for pos_aa in range(min(200, len(q_aligned_seq_aa))):
        curr_aa = q_aligned_seq_aa[pos_aa]

        if curr_aa == "-":
            output += "---"
        else:
            output += original_q_nt[pos_nt_no_gaps:pos_nt_no_gaps + 3]
            pos_nt_no_gaps += 3

        if pos_nt_no_gaps >= max_pos_nt_no_gaps:
            break

    return output


def compute_distance_based_on_local_alignment(query_info, target_info, hsp, **kwargs):
    # type: (Dict[str, Any], Dict[str, Any], HSP, Dict[str, Any]) -> float

    original_q_nt = get_value(kwargs, "original_q_nt", required=True)
    original_t_nt = get_value(kwargs, "original_t_nt", required=True)
    original_q_nt_offset = get_value(kwargs, "original_q_nt_offset", default=0)
    original_t_nt_offset = get_value(kwargs, "original_t_nt_offset", default=0)

    # aligned fragments (aa)
    q_aligned_seq_aa = hsp.query
    t_aligned_seq_aa = hsp.sbjct

    # indices of where alignment starts in original sequences
    q_start, q_end = hsp.query_start - 1, hsp.query_end - 2  # -2 to make inclusive
    t_start, t_end = hsp.sbjct_start - 1, hsp.sbjct_end - 1  # -2 to make inclusive

    # aligned fragments (nt)
    try:
        q_aligned_seq_nt = map_aligned_aa_to_aligned_nt(q_aligned_seq_aa, original_q_nt, q_start, q_end, offset_nt=original_q_nt_offset)
        t_aligned_seq_nt = map_aligned_aa_to_aligned_nt(t_aligned_seq_aa, original_t_nt, t_start, t_end, offset_nt=original_t_nt_offset)
    except ValueError:
        return 100  # FIXME: the hell is going on

    # compute distance metric
    try:
        distance = k2p_distance(q_aligned_seq_nt, t_aligned_seq_nt)
    except ValueError:
        distance = 100

    return distance


def create_info_for_query_target_pair(query_info, target_info, hsp, **kwargs):
    # type: (Dict[str, Any], Dict[str, Any], HSP, Dict[str, Any]) -> Dict[str, Any]
    output = {
        "evalue": hsp.expect,
    }

    for k, v in kwargs.items():
        output[k] = v

    source_to_info = {"q": query_info, "t": target_info}
    for source in ["q", "t"]:

        for key in ["genome", "accession", "left", "right", "strand", "upstream_left", "upstream_right", "upstream_strand", "offset"]:
            output["{}-{}".format(source, key)] = source_to_info[source][key]

    output["q-prot-pos-5prime-in-frag-msa"] = query_info["offset"] / 3
    output["q-nucl-pos-5prime-in-frag-msa"] = query_info["offset"]
    output["q-prot-position-of-5prime-in-msa-fragment-no-gaps"] = query_info["offset"] / 3
    output["q-nucl-position-of-5prime-in-msa-fragment-no-gaps"] = query_info["offset"]
    output["t-prot-pos-5prime-in-frag-msa"] = target_info["offset"] / 3
    output["t-nucl-pos-5prime-in-frag-msa"] = target_info["offset"]
    output["t-prot-position-of-5prime-in-msa-fragment-no-gaps"] = target_info["offset"] / 3
    output["t-nucl-position-of-5prime-in-msa-fragment-no-gaps"] = target_info["offset"]

    output["q-key"] = "{};{};{};{}".format(
        query_info["accession"], query_info["left"], query_info["right"], query_info["strand"]
    )

    #output["q-prot-msa"] = Seq(query_info["lorf_nt"]).translate()._data
    #output["t-prot-msa"] = Seq(target_info["lorf_nt"]).translate()._data

    output["q-lorf_nt"] = Seq(query_info["lorf_nt"])._data
    output["t-lorf_nt"] = Seq(target_info["lorf_nt"])._data

    return output


def compute_distance_based_on_global_alignment_from_sequences(q_sequence, t_sequence, q_sequence_nt, t_sequence_nt,
                                                              matrix):
    [q_align, t_align, _, _, _] = \
        global_alignment_aa_with_gap(q_sequence, t_sequence, matrix)

    q_align_nt = add_gaps_to_nt_based_on_aa(q_sequence_nt, q_align)
    t_align_nt = add_gaps_to_nt_based_on_aa(t_sequence_nt, t_align)

    # count number of positions without gaps
    len_without_gaps = sum([1 for i in range(len(q_align)) if q_align[i] != "-" and t_align[i] != "-"])

    try:
        distance = k2p_distance(q_align_nt, t_align_nt)
    except ValueError:
        distance = 100

    return (distance, len(q_align), len_without_gaps)


def parse_filter_and_convert_to_csv(pf_blast_results, pf_output, **kwargs):
    # type: (str, str, Dict[str, Any]) -> None

    hsp_criteria = get_value(kwargs, "hsp_criteria", None)
    pf_q_original_nt = get_value(kwargs, "pf_q_original_nt", required=True)
    pf_t_original_nt = get_value(kwargs, "pf_t_original_nt", required=True)
    pf_q_original_aa = get_value(kwargs, "pf_q_original_aa", required=True)
    pf_t_original_aa = get_value(kwargs, "pf_t_original_aa", required=True)
    distance_min = get_value(kwargs, "distance_min", 0.001, default_if_none=True)
    distance_max = get_value(kwargs, "distance_max", 0.4, default_if_none=True)

    # open csv file for writing
    try:
        f_output = open(pf_output, "w")
    except OSError as e:
        log.warning("Could not open csv file for writing converted blast output: {}".format(pf_output))
        raise e

    try:
        f_blast_results = open(pf_blast_results, "r")
    except OSError as e:
        log.warning("Could not open blast results file: {}".format(pf_blast_results))
        raise e

    # read original nucleotide sequences (for computing distances)
    q_original_sequences_nt = read_fasta_into_hash(pf_q_original_nt)
    t_original_sequences_nt = read_fasta_into_hash(pf_t_original_nt)

    # read original sequences for trying out pairwise alignment ;)
    q_original_sequences_aa = read_fasta_into_hash(pf_q_original_aa)
    t_original_sequences_aa = read_fasta_into_hash(pf_t_original_aa)

    matrix = matlist.blosum62
    import sbsp_alg.phylogeny
    sbsp_alg.phylogeny.add_stop_codon_to_blosum(matrix)

    # open blast stream
    records = NCBIXML.parse(f_blast_results)
    header_written = False

    # for each blast query
    for r in records:

        query_info = unpack_fasta_header(r.query)
        num_selected_targets_for_query = 0

        # for each alignment to a target protein for the current query
        for alignment in r.alignments:

            if num_selected_targets_for_query >= 100:
                logger.debug("Stopping at 100 targets (from {}) for query".format(len(r.alignments)))
                break

            hsp = select_representative_hsp(alignment, hsp_criteria)

            target_info = unpack_fasta_header(alignment.hit_id)

            original_q_nt = q_original_sequences_nt[r.query]
            original_t_nt = t_original_sequences_nt[alignment.hit_id]

            distance = compute_distance_based_on_local_alignment(query_info, target_info, hsp,
                                                                 original_q_nt=original_q_nt,
                                                                 original_t_nt=original_t_nt,
                                                                 **kwargs)

            original_q_aa = q_original_sequences_aa[r.query]
            original_t_aa = t_original_sequences_aa[alignment.hit_id]

            #global_distance, global_length, global_length_without_gaps = compute_distance_based_on_global_alignment_from_sequences(
            #    original_q_aa, original_t_aa, original_q_nt, original_t_nt, matrix
            #)
            global_distance = global_length = global_length_without_gaps = 0

            # FIXME: thresholds should be from input configuration files
            if distance > distance_min and distance < distance_max:
            #if True:
                num_selected_targets_for_query  += 1

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

                sorted_header = sorted(output_info.keys())

                # if header not yet written, write it
                if not header_written:
                    f_output.write("{}\n".format(",".join(sorted_header)))
                    header_written = True

                # write values in sorted order
                f_output.write("{}\n".format(
                    ",".join([str(output_info[x]) for x in sorted_header])
                ))

    f_output.close()


# def convert_blast_output_to_csv(pf_blast_output, pf_csv, select_best_alignment_per_qt_pair=True, delimiter=",", **kwargs):
#     # type: (str, str, bool) -> None
#
#     def get_lowest_evalue_of_hsps(alignment):
#         return np.min([hsp.expect for hsp in alignment.hsps])
#
#     def set_alignment_if_lowest_evalue(data, genome_name, alignment):
#         # type: (dict, str, dict) -> None
#
#         if genome_name not in data:
#             data[genome_name] = alignment
#         # otherwise, check if it has a lower E-value than what's already there
#         else:
#             evalue_existing = get_lowest_evalue_of_hsps(data[genome_name])
#             evalue_new = get_lowest_evalue_of_hsps(alignment)
#
#             if evalue_new < evalue_existing:
#                 data[genome_name] = alignment
#
#     # open output file for writing
#     with open(pf_csv, "w") as f_csv:
#
#         header = None
#
#         import sbsp_io.blast
#         records = sbsp_io.blast.read_hits(pf_blast_output)
#
#         for r in records:
#
#             q_def = r.query  # query definition line
#             q_genome = sbsp_general.general.get_genome_name_from_defition_line(q_def)
#             if select_best_alignment_per_qt_pair:
#                 # Structure:
#                 #   Key     : target genome name
#                 #   Value   : alignment with the lowest evalue for that target
#                 best_alignment_per_genome = {}
#
#                 # for each target genome, only select the best (lowest e-value) alignment
#                 for alignment in r.alignments:
#                     t_def = alignment.hit_id
#                     t_genome = sbsp_general.general.get_genome_name_from_defition_line(t_def)
#
#                     # keep track of the best alignment for this genome
#                     set_alignment_if_lowest_evalue(best_alignment_per_genome, t_genome, alignment)
#
#                 for t_genome in best_alignment_per_genome.keys():
#                     t_def = best_alignment_per_genome[t_genome].hit_id
#                     curr_alignment = best_alignment_per_genome[t_genome]
#                     evalue = get_lowest_evalue_of_hsps(curr_alignment)
#
#
#                     if not blast_on_query_database:
#                         # now print
#                         info = {
#                             "q-def": q_def, "t-def" : t_def
#                         }
#
#                         q_info = sbsp_general.general.expand_definition_line(q_def, key_prefix="q-")
#                         t_info = sbsp_general.general.expand_definition_line(t_def, key_prefix="t-")
#                     else:
#                         info = {
#                             "q-def": t_def, "t-def": q_def
#                         }
#
#                         q_info = sbsp_general.general.expand_definition_line(q_def, key_prefix="t-")
#                         t_info = sbsp_general.general.expand_definition_line(t_def, key_prefix="q-")
#
#                     def is_frame_shifted(tmp_info, key_prefix):
#                         # type: (Dict[str, Any], str) -> bool
#                         return (int(tmp_info[key_prefix+"right"]) - int(tmp_info[key_prefix+"left"]) + 1) % 3 != 0
#
#                     if is_frame_shifted(q_info, "q-") or is_frame_shifted(t_info, "t-"):
#                         continue
#
#
#                     info.update(q_info)
#                     info.update(t_info)
#
#                     info["evalue"] = evalue
#
#                     # if first line, write header
#                     if header is None:
#                         header = sorted(info.keys())
#                         line = "{}".format(delimiter).join(header) + "\n"
#                         f_csv.write(line)
#                     # if not first line, make sure header is consistent
#                     else:
#                         curr_header = sorted(info.keys())
#                         if curr_header != header:
#                             raise ValueError("CSV rows don't have the same columns")
#
#                     # write data
#                     line = "{}".format(delimiter).join([str(info[h]) for h in header]) + "\n"
#                     f_csv.write(line)


def append_sequences_to_file(sequences, f):
    # type: (Dict[str, Seq], IO[AnyStr]) -> None
    for header, sequence in sequences.items():
        f.write(">{}\n{}\n".format(
            header, sequence
        ))


def get_pf_sequences_for_genome(env, gi, **kwargs):
    # type: (Environment, GenomeInfo, Dict[str, Any]) -> str
    fn_sequences = get_value(kwargs, "fn_sequences", "sequence.fasta")
    return os.path.join(env['pd-data'], gi.name, fn_sequences)


def get_pf_labels_for_genome(env, gi, **kwargs):
    # type: (Environment, GenomeInfo, Dict[str, Any]) -> str
    fn_labels = get_value(kwargs, "fn_labels", "ncbi.gff")
    return os.path.join(env['pd-data'], gi.name, fn_labels)


def get_lorf(label, sequences):
    # type: (Label, Dict[str, Seq]) -> Seq

    if label.strand() == "+":

        curr_pos = label.left()
        pos_lorf = curr_pos

        while curr_pos >= 0:
            codon = sequences[label.seqname()][curr_pos:curr_pos + 3]._data
            if is_valid_start(codon, label.strand()):
                pos_lorf = curr_pos
            if is_valid_stop(codon, label.strand()):
                break

            curr_pos -= 3

        lorf_seq = sequences[label.seqname()][pos_lorf:label.right() + 1]

    else:

        curr_pos = label.right()
        pos_lorf = curr_pos
        seq_len = len(sequences[label.seqname()])

        while curr_pos < seq_len:
            codon = sequences[label.seqname()][curr_pos - 2:curr_pos + 1]._data
            if is_valid_start(codon, label.strand()):
                pos_lorf = curr_pos
            if is_valid_stop(codon, label.strand()):
                break

            curr_pos += 3

        lorf_seq = sequences[label.seqname()][label.left():pos_lorf + 1]

    return lorf_seq


def extract_labeled_sequence(label, sequences, **kwargs):
    # type: (Label, Dict[str, Seq], Dict[str, Any]) -> Seq
    reverse_complement = get_value(kwargs, "reverse_complement", False)
    lorf = get_value(kwargs, "lorf", False)

    if lorf:
        frag = get_lorf(label, sequences)
    else:
        frag = sequences[label.seqname()][label.left():label.right() + 1]

    if label.strand() == "-" and reverse_complement:
        frag = frag.reverse_complement()

    return frag

    # return "{}:tag={};11;:gc={}:pos={};{};{}:cds={};{};{}:type={}:key={};{};{}".format(
    #     label.seqname(),
    #     gi.name,
    #     gc,
    #     label.left() + 1,
    #     label.right() + 1,
    #     label.strand(),
    #     label.left() + 1,
    #     label.right() + 1,
    #     label.strand(),
    #     seq_type,
    #     label.seqname(),
    #     label.right() if label.strand() == "+" else label.left(),
    #     label.strand()
    # )


def get_upstream_label_per_label(labels):
    # type: (Labels) -> Dict[str, Label]

    sorted_labels = [l for l in labels.sort_by("left")]
    num_labels = len(sorted_labels)

    gene_key_to_upstream_label = dict()  # type: Dict[str, Label]

    for i, label in enumerate(sorted_labels):

        prev_label = None

        if label.strand() == "+":
            if i > 0:
                prev_i = i - 1
                prev_label = sorted_labels[prev_i]
        else:
            if i < num_labels - 1:
                prev_i = i + 1
                prev_label = sorted_labels[prev_i]

        gene_key_to_upstream_label[create_gene_key_from_label(label)] = prev_label

    return gene_key_to_upstream_label


def extract_labeled_sequences(sequences, labels, **kwargs):
    # type: (Dict[str, Seq], Labels, Dict[str, Any]) -> Dict[str, Seq]

    func_fasta_header_creator = get_value(kwargs, "func_fhc", None)
    kwargs_fasta_header_creator = get_value(kwargs, "kwargs_fhc", None)

    dict_labeled_sequences = dict()  # type: Dict[str, Seq]

    gene_key_to_upstream_label = get_upstream_label_per_label(labels)

    for i, label in enumerate(labels):
        labeled_sequence = extract_labeled_sequence(label, sequences, **kwargs)
        lorf_nt = extract_labeled_sequence(label, sequences, lorf=True, **kwargs)
        offset = len(lorf_nt) - (label.right() - label.left() + 1)

        upstream_label = gene_key_to_upstream_label[create_gene_key_from_label(label)]
        upstream_left = upstream_right = -1
        upstream_strand = ""

        if upstream_label is not None:
            upstream_left = upstream_label.left() + 1
            upstream_right = upstream_label.right() + 1
            upstream_strand = upstream_label.strand()

        fasta_header = str(i)
        if func_fasta_header_creator is not None:
            if kwargs_fasta_header_creator is not None:
                fasta_header = func_fasta_header_creator(label, offset=offset, lorf_nt=lorf_nt,
                                                         upstream_left=upstream_left,
                                                         upstream_right=upstream_right,
                                                         upstream_strand=upstream_strand,
                                                         **kwargs_fasta_header_creator)
            else:
                fasta_header = func_fasta_header_creator(label, offset=offset, lorf_nt=lorf_nt)

        dict_labeled_sequences[fasta_header] = labeled_sequence

    return dict_labeled_sequences


def extract_labeled_sequences_for_genome(env, gi, **kwargs):
    # type: (Environment, GenomeInfo, Dict[str, Any]) -> Dict[str, Seq]

    pf_sequences = get_pf_sequences_for_genome(env, gi)
    pf_labels = get_pf_labels_for_genome(env, gi, **kwargs)

    try:
        sequences = read_fasta_into_hash(pf_sequences)
        labels = read_labels_from_file(pf_labels, **kwargs)
    except IOError as e:
        log.warning("Could not read sequence/labels files for genome: {}".format(gi.name))
        raise e

    return extract_labeled_sequences(sequences, labels, **kwargs)


def translate_sequences_to_aa(sequences_nt):
    # type: (Dict[str, Seq]) -> Dict[str, Seq]

    sequences_aa = dict()
    for k, v in sequences_nt.items():
        try:
            v_aa = v.translate()
            sequences_aa[k] = v_aa
        except ValueError:
            log.warning("Could not translate sequence:\n{}".format(v))

    return sequences_aa


def dict_intersection_by_key(dict_a, dict_b):
    # type: (Dict, Dict) -> None
    """Removes elements from both dictionaries with keys not in both"""
    keys_a = set(dict_a.keys())
    keys_b = set(dict_b.keys())
    keys_intersection = keys_a.intersection(keys_b)

    keys_a_unique = keys_a.difference(keys_intersection)
    keys_b_unique = keys_b.difference(keys_intersection)

    def remove_keys_from_dict(d, keys):
        # type: (Dict, Iterable) -> None
        for k in keys:
            del d[k]

    remove_keys_from_dict(dict_a, keys_a_unique)
    remove_keys_from_dict(dict_b, keys_b_unique)


def extract_labeled_sequences_for_genomes(env, gil, pf_output, **kwargs):
    # type: (Environment, GenomeInfoList, str, Dict[str, Any]) -> str

    # open file for writing
    try:
        f_aa = open(pf_output, "w")

        for gi in gil:
            func_fasta_header_creator = pack_fasta_header
            kwargs_fasta_header_creator = {"gi": gi}
            try:
                sequences_nt = extract_labeled_sequences_for_genome(
                    env, gi,
                    func_fhc=func_fasta_header_creator,
                    kwargs_fhc=kwargs_fasta_header_creator,
                    **kwargs
                )
                sequences_aa = translate_sequences_to_aa(sequences_nt)

                # only keep sequences that have been translated
                dict_intersection_by_key(sequences_nt, sequences_aa)

                append_sequences_to_file(sequences_aa, f_aa)
            except IOError:
                pass

        f_aa.close()
    except OSError:
        log.warning("Could not open file for writing sequences:\n{}".format(pf_output))

    return pf_output


def run_blast_on_sequence_file(env, pf_q_aa, pf_db, pf_blast_output, **kwargs):
    # type: (Environment, str, str, str, Dict[str, Any]) -> None
    run_blast_alignment(pf_q_aa, pf_db, pf_blast_output, use_diamond=True, **kwargs)


def get_orthologs_from_files(env, pf_q_list, pf_t_list, pf_output, **kwargs):
    # type: (Environment, str, str, str, Dict[str, Any]) -> str

    sbsp_options = get_value(kwargs, "sbsp_options", SBSPOptions(env))  # type: SBSPOptions

    fn_q_labels = get_value(kwargs, "fn_q_labels", "ncbi.gff")
    fn_t_labels = get_value(kwargs, "fn_t_labels", "ncbi.gff")

    q_gil = GenomeInfoList.init_from_file(pf_q_list)
    t_gil = GenomeInfoList.init_from_file(pf_t_list)

    pd_work = env["pd-work"]

    # Extract data for blast run
    pf_q_aa = os.path.join(pd_work, "q.faa")
    pf_q_nt = os.path.join(pd_work, "q.fnt")
    pf_t_aa = os.path.join(pd_work, "t.faa")
    pf_t_nt = os.path.join(pd_work, "t.fnt")

    custom = {
        "reverse_complement": True,
        "ignore_frameshifted": True,
        "ignore_partial": True
    }

    extract_labeled_sequences_for_genomes(env, q_gil, pf_q_aa, fn_labels=fn_q_labels, **custom, **kwargs)
    extract_labeled_sequences_for_genomes(env, t_gil, pf_t_aa, fn_labels=fn_t_labels, **custom, **kwargs)

    pf_blast_db = os.path.join(pd_work, "blast.db")
    create_blast_database(pf_t_aa, pf_blast_db, seq_type="prot", use_diamond=True)  # FIXME: cleanup

    # Run blast
    pf_blast_results = os.path.join(pd_work, "blast.xml")
    run_blast_on_sequence_file(env, pf_q_aa, pf_blast_db, pf_blast_results, **kwargs)

    # Parse data, filter, and write to CSV
    parse_filter_and_convert_to_csv(pf_blast_results, pf_output,
                                    pf_q_original_nt=pf_q_nt,
                                    pf_t_original_nt=pf_t_nt,
                                    pf_q_original_aa=pf_q_aa,
                                    pf_t_original_aa=pf_t_aa,
                                    distance_min=sbsp_options.safe_get("filter-min-distance"),
                                    distance_max=sbsp_options.safe_get("filter-max-distance"),
                                    **kwargs)

    return pf_output
