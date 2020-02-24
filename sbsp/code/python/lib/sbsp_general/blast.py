from __future__ import print_function
import os
from subprocess import CalledProcessError
from typing import *

import numpy as np
import sbsp_general.general

# Generating commands
from sbsp_container.genome_list import GenomeInfoList
from sbsp_general import Environment
from sbsp_io.sequences import extract_genes_for_multiple_genomes


def gen_cmd_run_blastp(pf_q_sequences, pf_blast_db, pf_blast_out, use_diamond, **kwargs):
    # type: (str, str, str, bool, **str) -> str

    max_evalue = sbsp_general.general.get_value(kwargs, "max_evalue", None)

    if use_diamond:
        cmd = "diamond blastp -d {} -q {} -o {} --outfmt 5 --quiet -k 0 --subject-cover 80 --query-cover 80 ".format(
            pf_blast_db, pf_q_sequences, pf_blast_out
        )

        if max_evalue is not None:
            cmd += " --evalue {}".format(max_evalue)

    else:
        # TODO: add NCBI blast support
        raise NotImplemented("Please use diamond blast. NCBI blast not yet supported")

    return cmd


def gen_cmd_create_blast_database(pf_input, pf_blast_db, seq_type, use_diamond):
    # type: (str, str, str, bool) -> str

    if use_diamond:
        cmd = "diamond makedb --in " + pf_input + " -d " + pf_blast_db + " --quiet \n"

    else:
        # TODO: add NCBI blast support
        raise NotImplemented("Please use diamond blast. NCBI blast not yet supported")

    return cmd


# Running commands

def run_blast_alignment(pf_q_sequences, pf_blast_db, pf_blast_out, use_diamond, **kwargs):
    # type: (str, str, str, bool, **str) -> None

    cmd = gen_cmd_run_blastp(pf_q_sequences, pf_blast_db, pf_blast_out, use_diamond, **kwargs)
    try:
        sbsp_general.general.run_shell_cmd(cmd)
    except CalledProcessError:
        raise ValueError("Couldn't run blast")


def create_blast_database(pf_input, pf_blast_db, seq_type="nucl", use_diamond=True):
    # type: (str, str, str, bool) -> None
    cmd = gen_cmd_create_blast_database(pf_input, pf_blast_db, seq_type, use_diamond)
    sbsp_general.general.run_shell_cmd(cmd)


def run_blast(env, pf_q_list, pf_t_list, pf_blast_output, **kwargs):
    # type: (Environment, str, str, str, Dict[str, Any]) -> None

    fn_q_labels = sbsp_general.general.get_value(kwargs, "fn_q_labels", "ncbi.gff")
    fn_t_labels = sbsp_general.general.get_value(kwargs, "fn_t_labels", "ncbi.gff")
    fn_q_proteins = sbsp_general.general.get_value(kwargs, "fn_q_proteins", "q_proteins.faa")
    fn_t_proteins = sbsp_general.general.get_value(kwargs, "fn_t_proteins", "t_proteins.faa")

    fn_q_nucl = "q.fnt"
    fn_t_nucl = "t.fnt"

    fn_blast_db = sbsp_general.general.get_value(kwargs, "fn_blast_db", "db")
    clean = sbsp_general.general.get_value(kwargs, "clean", False)
    max_evalue = sbsp_general.general.get_value(kwargs, "max_evalue", None)

    pd_work = env["pd-work"]

    # get paths to files
    pf_q_proteins = os.path.join(pd_work, fn_q_proteins)
    pf_t_proteins = os.path.join(pd_work, fn_t_proteins)
    pf_blast_db = os.path.join(pd_work, fn_blast_db)

    # extract all proteins
    extract_genes_for_multiple_genomes(env, GenomeInfoList.init_from_file(pf_q_list), fn_q_labels, pf_q_proteins)
    extract_genes_for_multiple_genomes(env, GenomeInfoList.init_from_file(pf_t_list), fn_t_labels, pf_t_proteins)

    # create blast data base from target proteins
    create_blast_database(pf_t_proteins, pf_blast_db, seq_type="prot", use_diamond=True)

    # run blastp
    run_blast_alignment(pf_q_proteins, pf_blast_db, pf_blast_output, use_diamond=True, max_evalue=max_evalue)

    if clean:
        try:
            os.remove(pf_blast_db + ".dmnd")
        except OSError:
            pass

        try:
            os.remove(pf_q_proteins)
        except OSError:
            pass

        try:
            os.remove(pf_t_proteins)
        except OSError:
            pass


def run_blast_from_genome_names(env, q_genomes, t_genomes, fn_q_labels, fn_t_labels, pf_blast_out, **kwargs):
    # type: (dict, list[str], list[str], str, str, str, **str) -> None

    from sbsp_io.sequences import extract_genes_for_multiple_genomes

    fn_q_proteins = sbsp_general.general.get_value(kwargs, "fn_q_proteins", "q_proteins.faa")
    fn_t_proteins = sbsp_general.general.get_value(kwargs, "fn_t_proteins", "t_proteins.faa")
    fn_blast_db = sbsp_general.general.get_value(kwargs, "fn_blast_db", "db")
    clean = sbsp_general.general.get_value(kwargs, "clean", False)
    blast_on_query_database = sbsp_general.general.get_value(kwargs, "blast_on_query_database", False)

    pd_work = env["pd-work"]
    pd_data = env["pd-data"]

    pf_q_proteins = os.path.join(pd_work, fn_q_proteins)
    pf_t_proteins = os.path.join(pd_work, fn_t_proteins)

    pf_blast_db = os.path.join(pd_work, fn_blast_db)

    # extract all proteins
    extract_genes_for_multiple_genomes(env, q_genomes, fn_q_labels, pf_q_proteins)
    extract_genes_for_multiple_genomes(env, t_genomes, fn_t_labels, pf_t_proteins)

    if not blast_on_query_database:
        # create blast data base from target proteins
        create_blast_database(pf_t_proteins, pf_blast_db, seq_type="prot", use_diamond=True)

        # run blastp
        run_blast_alignment(pf_q_proteins, pf_blast_db, pf_blast_out, use_diamond=True)
    else:
        # create blast data base from query proteins
        create_blast_database(pf_q_proteins, pf_blast_db, seq_type="prot", use_diamond=True)

        # run blastp
        run_blast_alignment(pf_t_proteins, pf_blast_db, pf_blast_out, use_diamond=True)

    if clean:
        try:
            os.remove(pf_blast_db + ".dmnd")
        except OSError:
            pass

        try:
            os.remove(pf_q_proteins)
        except OSError:
            pass

        try:
            os.remove(pf_t_proteins)
        except OSError:
            pass


def add_to_dict(original, to_add, key_prefix=""):
    # type: (dict, dict, str) ->  None
    for key in to_add.keys():
        original["{}{}".format(key_prefix, key)] = to_add[key]


def convert_blast_output_to_csv(pf_blast_output, pf_csv, select_best_alignment_per_qt_pair=True, delimiter=",", **kwargs):
    # type: (str, str, bool) -> None

    blast_on_query_database = sbsp_general.general.get_value(kwargs, "blast_on_query_database", False)


    def get_lowest_evalue_of_hsps(alignment):
        return np.min([hsp.expect for hsp in alignment.hsps])

    def set_alignment_if_lowest_evalue(data, genome_name, alignment):
        # type: (dict, str, dict) -> None

        if genome_name not in data:
            data[genome_name] = alignment
        # otherwise, check if it has a lower E-value than what's already there
        else:
            evalue_existing = get_lowest_evalue_of_hsps(data[genome_name])
            evalue_new = get_lowest_evalue_of_hsps(alignment)

            if evalue_new < evalue_existing:
                data[genome_name] = alignment

    # open output file for writing
    with open(pf_csv, "w") as f_csv:

        header = None

        import sbsp_io.blast
        records = sbsp_io.blast.read_hits(pf_blast_output)

        for r in records:

            q_def = r.query  # query definition line
            q_genome = sbsp_general.general.get_genome_name_from_defition_line(q_def)
            if select_best_alignment_per_qt_pair:
                # Structure:
                #   Key     : target genome name
                #   Value   : alignment with the lowest evalue for that target
                best_alignment_per_genome = {}

                # for each target genome, only select the best (lowest e-value) alignment
                for alignment in r.alignments:
                    t_def = alignment.hit_id
                    t_genome = sbsp_general.general.get_genome_name_from_defition_line(t_def)

                    # keep track of the best alignment for this genome
                    set_alignment_if_lowest_evalue(best_alignment_per_genome, t_genome, alignment)

                for t_genome in best_alignment_per_genome.keys():
                    t_def = best_alignment_per_genome[t_genome].hit_id
                    curr_alignment = best_alignment_per_genome[t_genome]
                    evalue = get_lowest_evalue_of_hsps(curr_alignment)


                    if not blast_on_query_database:
                        # now print
                        info = {
                            "q-def": q_def, "t-def" : t_def
                        }

                        q_info = sbsp_general.general.expand_definition_line(q_def, key_prefix="q-")
                        t_info = sbsp_general.general.expand_definition_line(t_def, key_prefix="t-")
                    else:
                        info = {
                            "q-def": t_def, "t-def": q_def
                        }

                        q_info = sbsp_general.general.expand_definition_line(q_def, key_prefix="t-")
                        t_info = sbsp_general.general.expand_definition_line(t_def, key_prefix="q-")

                    def is_frame_shifted(tmp_info, key_prefix):
                        # type: (Dict[str, Any], str) -> bool
                        return (int(tmp_info[key_prefix+"right"]) - int(tmp_info[key_prefix+"left"]) + 1) % 3 != 0

                    if is_frame_shifted(q_info, "q-") or is_frame_shifted(t_info, "t-"):
                        continue


                    info.update(q_info)
                    info.update(t_info)

                    info["evalue"] = evalue

                    # if first line, write header
                    if header is None:
                        header = sorted(info.keys())
                        line = "{}".format(delimiter).join(header) + "\n"
                        f_csv.write(line)
                    # if not first line, make sure header is consistent
                    else:
                        curr_header = sorted(info.keys())
                        if curr_header != header:
                            raise ValueError("CSV rows don't have the same columns")

                    # write data
                    line = "{}".format(delimiter).join([str(info[h]) for h in header]) + "\n"
                    f_csv.write(line)










