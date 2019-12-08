#!/usr/bin/env bash

# Karl Gemayel
# Georgia Institute of Technology


function download_genomic_data_from_assembly_summary() {
    # Parameters:
    #   1: file containing assembly summary (with ftp links)
    #   2: root directory in which data will be download
    fnassembly="$1"
    p_data="$2"
    fnassembly_out="$3"

    # if assembly file doesn't exist, quit
    if [ ! -f $fnassembly ]; then
        echo "Assembly file doesn't exist. Cannot download data..."
        return
    fi

    grep "accession" $fnassembly | head -n 1 > $fnassembly_out

    # create data directory
    mkdir -p ${p_data}

    # get genomic data for tax ids
    cat $fnassembly | while read -r line; do

        ftplink=$(echo "$line" | perl -n -e  '/(ftp:[^\s]+)/ && print "$1\n"')


        # if ftp link not found, move to next line
        if [ -z "$ftplink" ]; then
            continue 
        fi

        gcf=$(echo "$line" | cut -f1)
        acc=$(echo "$line" | cut -f16 | tr " " "_")
        
        gcfacc="${gcf}_${acc}"

        # if directory exists, skip
        if [ -d "${p_data}/$gcfacc" ]; then

            rm -r "${p_data}/$gcfacc"
            #continue
        fi

        echo "$ftplink"

        fngenome="${gcfacc}_genomic.fna"
        fngff="${gcfacc}_genomic.gff"

        # construct links to fasta and gff files
        ftpgenome="${ftplink}/${fngenome}.gz"
        ftpgff="${ftplink}/${fngff}.gz"

        curr_dir=${p_data}/$gcfacc
        mkdir -p ${curr_dir}
        
        wget --quiet $ftpgenome -P ${curr_dir}

        if [ $? -ne 0 ]; then
            rm -r $curr_dir
            continue
        fi

        wget --quiet $ftpgff -P ${curr_dir}
        if [ $? -ne 0 ]; then
            rm -r $curr_dir
            continue
        fi
        
        gunzip ${curr_dir}/${fngenome}.gz
        gunzip ${curr_dir}/${fngff}.gz
        
        mv ${curr_dir}/$fngenome ${curr_dir}/sequence.fasta
        mv ${curr_dir}/$fngff ${curr_dir}/ncbi.gff
    
        echo "$line" >> $fnassembly_out

        
    done
}


function create_list_from_assembly_summary() {
    fnassembly="$1"
    fnlist="$2"

    #FIXME: scripts assumes everything is genetic code 11

    # if assembly file doesn't exist, quit
    if [ ! -f $fnassembly ]; then
        echo "Assembly file doesn't exist. Cannot download data..."
        return
    fi

    if [ -f $fnlist ]; then
        rm $fnlist
    fi

    touch $fnlist


    # get genomic data for tax ids
    cat $fnassembly | while read -r line; do

        ftplink=$(echo "$line" | perl -n -e  '/(ftp:[^\s]+)/ && print "$1\n"')

        # if ftp link not found, move to next line
        if [ -z "$ftplink" ]; then
            continue 
        fi

        gcf=$(echo "$line" | cut -f1)
        acc=$(echo "$line" | cut -f16 | tr " " "_")
        
        gcfacc="${gcf}_${acc}"

        echo -e "$gcfacc\t11" >> $fnlist
        
    done
}


function get_human_readable_names() {
    fnassembly="$1"
    fnlist="$2"

    awk '{print $1}' $fnlist | while read -r line; do
        gcfacc=$line

        gcf=$(echo "$gcfacc" | cut -f-2 -d "_");
        acc=$(echo "$gcfacc" | cut -f3- -d "_" | tr "_" " ");

        echo -e "Here: $gcfacc\t$gcf\t$acc"

        grep "^$gcf" $fnassembly
    done
}

function ncbi_setup_dataset_from_assembly() {
    # Download and setup a dataset from NCBI
    local fn_filtered_assembly="$1"
    local fn_list="$2"
    local path_data="$3"
    local pf_assembly_out="$4"
    
    # download data
    download_genomic_data_from_assembly_summary ${fn_filtered_assembly} ${path_data} ${pf_assembly_out}
    
    # make list
    create_list_from_assembly_summary ${pf_assembly_out} ${fn_list}
    
    # fix fasta header label in sequence files
    set_fasta_header_with_accession ${path_data} ${fn_list}
    
    # extract ribosomal labels
    extract_ribosomal_labels ${path_data} ${fn_list} ncbi.gff ncbi_ribosomal.gff
}



function ncbi_setup_dataset_by_name() {
    # Download and setup a dataset from NCBI
    # E.g. Escherichia 
    local fn_assembly="$1"
    local name_contains="$2"
    local fn_filtered_assembly="$3"
    local fn_list="$4"
    local path_data="$5"
    local pf_assembly_out="$6"
    
    # extract species from assembly
    $bin/extract-from-assembly-summary_py.sh --summary ${fn_assembly} --out $fn_filtered_assembly --key organism_name --key-value ${name_contains} --full-genome --filter-type contains --max-records 20

    # download data
    download_genomic_data_from_assembly_summary ${fn_filtered_assembly} ${path_data} ${pf_assembly_out}
    
    # make list
    # FIXME: use pf_assembly_out
    create_list_from_assembly_summary ${fn_filtered_assembly} ${fn_list}
    
    # fix fasta header label in sequence files
    set_fasta_header_with_accession ${path_data} ${fn_list}
    
    # extract ribosomal labels
    extract_ribosomal_labels ${path_data} ${fn_list} ncbi.gff ncbi_ribosomal.gff
}


function ncbi_get_assembly_entries_from_names() {
    # Get assembly entries for each name 
    local pf_names="$1"
    local pf_assembly="$2"
    local num_per_name="$3"
}



































