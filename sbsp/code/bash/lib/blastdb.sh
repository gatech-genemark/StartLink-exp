#!/bin/bash

# Karl Gemayel
# Georgia Institute of Technology
#
# Created 24/10/2018

function create_blast_db_from_sequence_file() {
    # Creates a blast database from a file containing fasta sequences.
    # Parameters
    #   1: sequence file
    #   2: Path to database
    #   3: Type: prot, nucl
    fnsequences=$1
    pathDB=$2
    type=$3

    makeblastdb -in $fnsequences -dbtype $type -out $pathDB
}

function create_blast_db_from_annotation() {
    # Create a blast database from the annotations of a list of genomes
    # Parameters
    #   1: list of genomes: Format: Name gcode
    #   2: path to genome data
    #   3: Type: prot, nucl
    #   4: Path Work
    #   5: Database name
    create_blast_db $1 $2 $3 $4 $5 $6;
}

function create_blast_db_from_gms2() {
    # Create a blast database from the annotations of a list of genomes
    # Parameters
    #   1: list of genomes: Format: Name gcode
    #   2: path to genome data
    #   3: directory containing gms2
    #   4: Type: prot, nucl
    #   5: Path Work
    #   6: Database name
    if [[ $# < 5 ]]; then
        echo "Not enough input arguments";
        return;
    fi

    fngenomes=$1;
    pathGenomesData=$2;
    dirGMS2=$3;
    type=$4;
    pathWork=$5;
    dbname=$6;

    # check if upstream length inputted
    upstrLen=0
    if [[ $# == 7 ]]; then
        upstrLen=$7;
    fi



    # set option for amino acids if prot type is set
    aaOption=""
    if [[ $type -eq "prot" ]]; then
        aaOption="--aa";
    fi

    # temporary file for storing ORFs
    fntmpORFs=$pathWork/orfs.fasta;
    rm -f $fntmpORFs;

    while read -r line; do
        gcfid=$(echo "$line" | awk '{print $1}');
        gcode=$(echo "$line" | awk '{print $2}');

        pathGenome=$pathGenomesData/$gcfid;
        fnseq=$pathGenome/sequence.fasta;
        fnlabels=$pathGenome/runs/${dirGMS2}/gms2.gff;

        gc="";#$(probuild --stat --gc --seq $fnseq | awk '{print $3}');
        tag="$gcfid;$gcode;$gc;gms2";

        # extract ORFs (amino acid sequences)
        biogem utilities extract-orf --sequences $fnseq --labels $fnlabels --gcode $gcode $aaOption --tag $tag --upstream-len $upstrLen --exclude-no-start >> $fntmpORFs;

    done < "$fngenomes"

    pathDB=$pathWork/$dbname
    create_blast_db_from_sequence_file $fntmpORFs $pathDB $type;
}


function extract_labeled_seqs_for_blastdb {
    # Create a blast database from the annotations of a list of genomes
    # Parameters
    #   1: list of genomes: Format: Name gcode
    #   2: path to genome data
    #   3: Type: prot, nucl
    #   4: Path Work
    #   5: output filename
    #   6: labels filename
    if [[ $# < 6 ]]; then 
        echo "Not enough input arguments";
        return;
    fi

    fngenomes=$1;
    pathGenomesData=$2;
    type=$3;
    pathWork=$4;
    fntmpORFs=$5;
    fnlabelsbase=$6

    echo "$pathGenomesData"

    # check if upstream length inputted
    upstrLen=0
    if [[ $# == 7 ]]; then
        upstrLen=$7;
    fi



    # set option for amino acids if prot type is set
    aaOption=""
    if [[ $type -eq "prot" ]]; then
        aaOption="--aa";
    fi

    # file for storing ORFs
    rm -f $pathWork/$fntmpORFs;

    while read -r line; do
        gcfid=$(echo "$line" | awk '{print $1}');
        gcode=$(echo "$line" | awk '{print $2}');

        pathGenome=$pathGenomesData/$gcfid;
        fnseq=$pathGenome/sequence.fasta;
        fnlabels=$pathGenome/$fnlabelsbase;

        gc="";#$(probuild --stat --gc --seq $fnseq | awk '{print $3}');
        tag="$gcfid;$gcode;$gc";

        # extract ORFs (amino acid sequences)
        # echo "biogem utilities extract-orf --sequences $fnseq --labels $fnlabels --gcode $gcode $aaOption --tag \"$tag\" --upstream-len $upstrLen --exclude-no-start >> $pathWork/$fntmpORFs;"
        biogem utilities extract-orf --sequences $fnseq --labels $fnlabels --gcode $gcode $aaOption --tag "$tag" --upstream-len $upstrLen >> $pathWork/$fntmpORFs;

    done < "$fngenomes"
}

function create_blast_db() {
    # Create a blast database from the annotations of a list of genomes
    # Parameters
    #   1: list of genomes: Format: Name gcode
    #   2: path to genome data
    #   3: Type: prot, nucl
    #   4: Path Work
    #   5: Database name
    if [[ $# < 5 ]]; then 
        echo "Not enough input arguments";
        return;
    fi

    fngenomes=$1;
    pathGenomesData=$2;
    type=$3;
    pathWork=$4;
    dbname=$5;

    # check if upstream length inputted
    upstrLen=0
    if [[ $# == 6 ]]; then
        upstrLen=$6;
    fi



    # set option for amino acids if prot type is set
    aaOption=""
    if [[ $type -eq "prot" ]]; then
        aaOption="--aa";
    fi

    # temporary file for storing ORFs
    fntmpORFs=$pathWork/orfs.fasta;
    rm -f $fntmpORFs;

    while read -r line; do
        gcfid=$(echo "$line" | awk '{print $1}');
        gcode=$(echo "$line" | awk '{print $2}');

        pathGenome=$pathGenomesData/$gcfid;
        fnseq=$pathGenome/sequence.fasta;
        fnlabels=$pathGenome/ncbi.gff;

        gc="";#$(probuild --stat --gc --seq $fnseq | awk '{print $3}');
        tag="$gcfid;$gcode;$gc";

        # extract ORFs (amino acid sequences)
        biogem utilities extract-orf --sequences $fnseq --labels $fnlabels --gcode $gcode $aaOption --tag $tag --upstream-len $upstrLen >> $fntmpORFs;

    done < "$fngenomes"

    pathDB=$pathWork/$dbname
    create_blast_db_from_sequence_file $fntmpORFs $pathDB $type;
}
