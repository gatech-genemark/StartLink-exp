import os
from Bio import SeqIO
import Bio.Seq
import sbsp_general.general
from typing import *

from sbsp_container.genome_list import GenomeInfoList
from sbsp_general import Environment
from sbsp_general.general import get_value, except_if_not_in_set
import sbsp_general.labels
from sbsp_general.general import create_gene_key
from sbsp_io.general import write_string_to_file


def read_fasta_into_hash(fnsequences):
    # type: (str) -> dict[str, Bio.Seq.Seq]
    # read all sequences
    sequences = {}  # will contain all protein sequences from file
    for record in SeqIO.parse(fnsequences, "fasta"):
        sequences[record.name.strip().split()[0]] = record.seq

    return sequences


def write_fasta_hash_to_file(fasta, pf_output):
    # type: (Dict[str, Any], str) -> None

    output = ""
    for header in fasta.keys():
        output += ">{}\n{}\n".format(
            header, fasta[header]
        )

    write_string_to_file(output, pf_output)


def extract_genes_for_multiple_genomes(env, genomes, fn_labels, pf_output, seq_type="prot"):
    # type: (Environment, GenomeInfoList, str, str, str) -> None

    # FIXME: remove biogem dependency - do inhouse

    def write_genome_names_to_file(genome_names, pf_out):
        # type: (GenomeInfoList, str) -> None
        with open(pf_out, "w") as f:
            for g in genome_names:
                f.write("{} {}\n".format(g.name, '11'))

    pf_genomes = os.path.join(env["pd-work"], "tmp_genomes.list")

    write_genome_names_to_file(genomes, pf_genomes)

    cmd = "source {}/config.sh".format(env["pd-base"])
    cmd += "\n" + "for f in {}/*.sh; do source $f; done".format(env["pd-bash-lib"])
    cmd += "\nextract_labeled_seqs_for_blastdb {} {} {} {} {} {}\n".format(
        pf_genomes,
        env["pd-data"],
        seq_type,
        env["pd-work"],
        os.path.basename(pf_output),
        fn_labels
    )

    sbsp_general.general.run_shell_cmd(cmd)


def extract_genes_for_multiple_genomes_deprecated(env, genomes, fn_labels, pf_output, seq_type="prot"):
    # type: (iter, list[str], str, str, str) -> None

    def write_genome_names_to_file(genome_names, pf_out):
        with open(pf_out, "w") as f:
            for g in genome_names:
                f.write("{} {}\n".format(g, '11'))

    pf_genomes = os.path.join(env["pd-work"], "tmp_genomes.list")

    write_genome_names_to_file(genomes, pf_genomes)

    cmd = "source {}/config.sh".format(env["pd-base"])
    cmd += "\n" + "for f in {}/*.sh; do source $f; done".format(env["pd-bash-lib"])
    cmd += "\nextract_labeled_seqs_for_blastdb {} {} {} {} {} {}\n".format(
        pf_genomes,
        env["pd-data"],
        seq_type,
        env["pd-work"],
        os.path.basename(pf_output),
        fn_labels
    )

    sbsp_general.general.run_shell_cmd(cmd)


def read_sequences_from_genome_labels_pairs(env, genomes_labels_pairs, **kwargs):
    # type: (Environment, dict[str, sbsp_general.labels.Labels], **str) -> dict[str, dict[str, str]]

    sequence_type = get_value(kwargs, "sequence_type", "both")
    upstream_length_nt = get_value(kwargs, "upstream_length_nt", None)
    downstream_length_nt = get_value(kwargs, "downstream_length_nt", None)
    limit_stop_in_upstream = get_value(kwargs, "limit_stop_in_upstream", True)
    leave_out_gene_stop = get_value(kwargs, "leave_out_gene_stop", False)
    limit_upstream_to_first_candidate = get_value(kwargs, "limit_upstream_to_first_candidate", False)

    force_multiple_of_3 = get_value(kwargs, "force_multiple_of_3", True)

    except_if_not_in_set(sequence_type, ["both", "nucl", "prot"])

    upstream_length_nt = 180

    # if upstream_length_nt is not None:
    #     if downstream_length_nt > 0 and upstream_length_nt % 3 != 0:
    #         raise ValueError("Upstream length nt ({}) is not a multiple of 3".format(upstream_length_nt))

    # if downstream_length_nt is not None:
    #     if downstream_length_nt > 0 and downstream_length_nt % 3 != 0:
    #         raise ValueError("Downstream length nt ({}) is not a multiple of 3".format(downstream_length_nt))

    key_to_sequence = dict()

    pd_data = env["pd-data"]

    # do it genome after genome, so as to only read each sequence file once
    for genome_name in genomes_labels_pairs.keys():

        # read in sequence file
        pf_seq = os.path.join(pd_data, genome_name, "sequence.fasta")

        genome_seqs = read_fasta_into_hash(pf_seq)

        # for each label for that genome
        for label in genomes_labels_pairs[genome_name]:
            left = label.coordinates().left
            right = label.coordinates().right
            strand = label.coordinates().strand
            seqname = label.seqname()

            if seqname not in genome_seqs.keys():
                raise ValueError("Couldn't find accession number ({}) for genome ({})".format(seqname, genome_name))

            # get coordinates of fragment (in case upstream/downstream lengths are set)
            frag_left = left
            frag_right = right

            if leave_out_gene_stop:
                if strand == "+":
                    frag_right -= 3
                else:
                    frag_left += 3

            pos_5prime_in_frag_nt = 0

            if upstream_length_nt is not None or downstream_length_nt is not None:
                if strand == "+":
                    if upstream_length_nt is not None:
                        frag_left = max(0, frag_left - upstream_length_nt)

                        # adjust frag_left to make sure it is is multiple of 3 nt from left
                        #                        rem = (left - frag_left) % 3
                        #                        if rem != 0:
                        #                            frag_left += rem
                        #                            assert ((left - frag_left) % 3 == 0)

                        pos_5prime_in_frag_nt = left - frag_left

                    if downstream_length_nt is not None:
                        if downstream_length_nt > 0:
                            frag_right = min(len(genome_seqs[seqname]) - 1, left + 2 + downstream_length_nt)
                        else:
                            frag_right = min(len(genome_seqs[seqname]) - 1, left + downstream_length_nt)
                #
                #                        if downstream_length_nt > 0:
                #                            rem = (frag_right - left) % 3
                #                            if rem != 0:
                #                                frag_right -= rem
                #                                assert((frag_right - left) % 3 == 0)
                else:
                    if upstream_length_nt is not None:
                        frag_right = min(len(genome_seqs[seqname]) - 1, frag_right + upstream_length_nt)

                        # rem = (frag_right - right) % 3
                        # if rem != 0:
                        #    frag_right = frag_right -  rem
                        #    assert ((frag_right - right) % 3 == 0)

                        pos_5prime_in_frag_nt = frag_right - right

                    if downstream_length_nt is not None:
                        if downstream_length_nt > 0:
                            frag_left = max(0, right - 2 - downstream_length_nt)
                        else:
                            frag_left = max(0, right - downstream_length_nt)

                        # if downstream_length_nt >= 0:

                        #    rem = (right - frag_left) % 3
                        #    if rem != 0:
                        #        frag_left += rem

            frag_nt = genome_seqs[seqname][frag_left:frag_right + 1]

            if strand == "-":
                frag_nt = frag_nt.reverse_complement()

            if limit_stop_in_upstream:

                # find stop codon in upstream region
                actual_upstream_length = pos_5prime_in_frag_nt

                pos_of_stop = None
                for i in range(0, pos_5prime_in_frag_nt, 3):
                    if frag_nt[i:i + 3] in ["TAA", "TGA", "TAG"]:
                        pos_of_stop = i

                if pos_of_stop is not None:
                    frag_nt = frag_nt[pos_of_stop + 3:]
                    pos_5prime_in_frag_nt = pos_5prime_in_frag_nt - pos_of_stop - 3

            if limit_upstream_to_first_candidate:

                actual_upstream_length = pos_5prime_in_frag_nt

                pos_of_candidate = None
                for i in range(0, pos_5prime_in_frag_nt + 1, 3):
                    if frag_nt[i:i + 3] in ["ATG", "GTG", "TTG"]:
                        pos_of_candidate = i
                        break

                if pos_of_candidate is not None:
                    frag_nt = frag_nt[pos_of_candidate:]
                    pos_5prime_in_frag_nt = pos_5prime_in_frag_nt - pos_of_candidate

            pos_5prime_in_frag_aa = pos_5prime_in_frag_nt / 3
            key = create_gene_key(genome_name, seqname, left + 1, right + 1, strand)

            if force_multiple_of_3:
                if len(frag_nt) % 3 != 0:
                    continue

            if key not in key_to_sequence.keys():
                key_to_sequence[key] = dict()



            if sequence_type == "nucl" or sequence_type == "both":
                key_to_sequence[key]["nucl"] = frag_nt._data

                # add position of original 5prime in frag
                key_to_sequence[key]["nucl-pos-5prime-in-frag"] = pos_5prime_in_frag_nt  # FIXME: change nucl to nt

            if sequence_type == "prot" or sequence_type == "both":
                frag_aa = frag_nt.translate()
                key_to_sequence[key]["prot"] = frag_aa._data
                key_to_sequence[key]["prot-pos-5prime-in-frag"] = pos_5prime_in_frag_aa  # FIXME: change prot to aa

    return key_to_sequence
