#!/bin/bash

# Option default
VERBOSE=0
DEBUG=2
INFO=1

database=
ancestor=
num_per_taxid=1

# getopts string
# This string needs to be updated with the single character options (e.g. -f)
opts="a:d:vn:"

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
-a, --ancestor; Name of list containing query genomes (without .list suffix)
-d, --database; Database type (options: genbank, refseq)
-n, --num-per-taxid; Number of targets per taxid
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
                -a|--ancestor)      ancestor=$2; shift;;
                -d|--database)      database=$2; shift;;
                -n|--num-per-taxid) num_per_taxid=$2; shift;;
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
([ -z "$ancestor" ] || [ -z "$database" ]) && usage

function generate_pbs_string() {
    ancestor="$1"
    database="$2"
    ancestor_valid="$3"
    
    echo "#PBS -N $ancestor_valid"
    echo "#PBS -j oe"
    echo "#PBS -l nodes=1:ppn=1"
    echo "#PBS -l walltime=07:00:00:00"
    echo "#PBS -V"
    echo "#PBS -W umask=002"
    echo "export PATH=\"/home/karl/anaconda/envs/sbsp/bin:\$PATH\""
    echo "PBS_O_WORKDIR=/storage4/karl/sbsp/biogem/sbsp/tmp/downloading"
    echo "cd \$PBS_O_WORKDIR"

    echo "$bin/download-data-by-ancestor_py.sh  --ancestor-id \"$ancestor\" --ancestor-id-type name_txt --pf-taxonomy-tree tree --pf-assembly-summary $metadata/assembly_summary_${database}.txt --pd-output $data --pf-output-list $lists/${database}_${ancestor_valid}.list --number-per-taxid ${num_per_taxid} --favor-assembly-level-order --valid-assembly-levels \"Complete Genome\" Scaffold Contig"
}

ancestor_valid=$(echo "$ancestor" | tr '[:upper:]' '[:lower:]' | tr '/' "_" | tr " " "_" | tr "-" "_")
pf_pbs=${database}_${ancestor_valid}.pbs

generate_pbs_string "$ancestor" "$database" "${ancestor_valid}" > $pf_pbs

qsub $pf_pbs

