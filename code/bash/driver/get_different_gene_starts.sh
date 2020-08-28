#!/bin/bash

# Option default
VERBOSE=0
DEBUG=2
INFO=1
fn_q_labels=verified.gff
fn_q_labels_true=verified.gff
fn_t_labels=ncbi.gff
pbs_conf=defaults
steps=""

# getopts string
# This string needs to be updated with the single character options (e.g. -f)
opts="fvo:"

# Gets the command name without path
cmd(){ echo `basename $0`; }

function verbose() {
    [ ! $# -eq 1 ] && return 1

    level="$1"

    if [ $VERBOSE -ge $level ]; then
        return 0
    else
        return 1
    fi    
}

# Help command output
usage(){
echo "\
`cmd` -qts [OPTION...]
-l, --pf-genomes; Name of list containing query genomes (without .list suffix)
-v, --verbose; Enable verbose output (include multiple times for more
             ; verbosity, e.g. -vvv)
" | column -t -s ";"
}

# Error message
error(){
    echo "`cmd`: invalid option -- '$1'";
    echo "Try '`cmd` -h' for more information.";
    exit 1;
}

# There's two passes here. The first pass handles the long options and
# any short option that is already in canonical form. The second pass
# uses `getopt` to canonicalize any remaining short options and handle
# them
for pass in 1 2; do
    while [ -n "$1" ]; do
        case $1 in
            --) shift; break;;
            -*) case $1 in
                -l|--pf-genomes)    pf_genomes=$2; shift;;
                -v|--verbose)       VERBOSE=$(($VERBOSE + 1));;
                --*)                error $1;;
                -*)                 if [ $pass -eq 1 ]; then ARGS="$ARGS $1";
                                    else error $1; fi;;
                esac;;
            *)  if [ $pass -eq 1 ]; then ARGS="$ARGS $1";
                else error $1; fi;;
        esac
        shift
    done
    if [ $pass -eq 1 ]; then ARGS=`getopt $opts $ARGS`
        if [ $? != 0 ]; then usage; exit 2; fi; set -- $ARGS
    fi
done


# Check options
[ -z "$pf_genomes" ] && usage


# check that files exist
[ ! -f "${pf_genomes}" ] && { echo "File doesn't exist: $pf_genomes"; exit 1; };

cat ${pf_genomes} | awk '{if (NR > 1) print}' | while read -r line; do
    gcfid=$(echo "$line"| cut -f1 -d,)
    pf_ncbi=/storage4/karl/sbsp/similarity-based-start-prediction/data/large/$gcfid/ncbi.gff
    pf_gms2=/storage4/karl/sbsp/similarity-based-start-prediction/data/large/$gcfid/runs/gms2/gms2.gff
    pf_diff=/storage4/karl/sbsp/similarity-based-start-prediction/data/large/$gcfid/diff_5prime.gff
    
    [ ! -f "${pf_gms2}" ] && { echo "GMS2 run for '$gcfid' doesn't exist."; continue; };
    
    compp -a $pf_ncbi -b $pf_gms2 -q -n -S -L | grep -v "#" > $pf_diff
done

