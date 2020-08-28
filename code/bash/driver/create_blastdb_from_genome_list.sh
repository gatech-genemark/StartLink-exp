#!/bin/bash

# Option default
VERBOSE=0
DEBUG=2
INFO=1
pbs_conf=$config/pbs_defaults.conf

# getopts string
# This string needs to be updated with the single character options (e.g. -f)
opts="vp:"

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
--pf-genome-list; Name of list containing query genomes (without .list suffix)
--pf-blast-db; Name of output blast file
-t; Name of list containing target genomes (without .list suffix)
-p, --pf-pbs-conf;  Configuration number of PBS options (found in pbs_NUM.conf) (default: ${pbs_conf})
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
                --pf-genome-list)   pf_genome_list=$2; shift;;
                --pf-blast-db)      pf_blast_db=$2; shift;;
                -p|--pf-pbs-conf)   pf_pbs_conf=$2; shift;;
                
                --tag)              tag=$2; shift;;
                -v|--verbose)       VERBOSE=$(($VERBOSE + 1));;
                -e)                 email_on_complete=1; shift;;
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
([ -z "$pf_genome_list" ]) && usage

# check that files exist
[ ! -f "${pf_genome_list}" ] && { echo "File doesn't exist: $pf_q_list"; exit 1; };
[ ! -f "${pf_pbs_conf}" ] && { echo "File doesn't exist: $pf_pbs_conf"; exit 1; };

pd_work=$(pwd)
pf_sequences_tmp=$(mktemp --tmpdir=${pd_work})

$bin/extract_annotated_sequences_py.sh --pf-genome-list $pf_genome_list --pf-output $pf_sequences_tmp  --pf-pbs-options $pf_pbs_conf
diamond makedb --in $pf_sequences_tmp --db $pf_blast_db
rm $pf_sequences_tmp
