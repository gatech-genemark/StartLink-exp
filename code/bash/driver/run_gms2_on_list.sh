#!/bin/bash

# Option default
VERBOSE=0
DEBUG=2
INFO=1
pf_list=
dn_output=gms2
domain=

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
-l, --list; File containing list of genomes
-o, --dn-output; Name of output directory under genome/runs  (default: ${dn_output})
-d, --domain; Domain of species (options: bacteria, archaea)
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
                -l|--list)          pf_list=$2; shift;;
                -o|--dn-output)     dn_output=$2; shift;;
                -d|--domain)        domain=$2; shift;;
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


# # Check options
([ -z "${pf_list}" ] || [ -z "${dn_output}" ] || [ -z "$domain" ]) && usage
([ "$domain" != "bacteria" ] && [ "$domain" != "archaea" ]) && { echo "Domain must be archaea or bacteria"; exit 1; }
 
# check that files exist
[ ! -f "${pf_list}" ] && { echo "File doesn't exist: $pf_list"; exit 1; };

# run command
function create_pbs_header() {
    gcfid="$1"
    pd_output="$2"
    echo "\
#PBS -N $gcfid
#PBS -o ${pd_output}/gms2.oe
#PBS -j oe
#PBS -l nodes=1:ppn=1
#PBS -W umask=002

PBS_O_WORKDIR=$pd_output
cd \$PBS_O_WORKDIR

ls | grep -v ".pbs" | xargs rm

"
}

total=$(wc -l $pf_list | awk '{print $1 - 1}')
current=0

awk '{if (NR > 1) print}' $pf_list | while read -r line; do
    gcfid=$(echo "$line" | awk -F "," '{print $1}')
    gencode=$(echo "$line" | awk -F "," '{print $2}')

    current=$((current + 1))
    echo -ne "Progress: $current / $total\r"

    pd_gcfid=$data/$gcfid
    pd_gcfid_runs=$pd_gcfid/runs
    pd_output=$pd_gcfid_runs/$dn_output
    pf_sequence=$pd_gcfid/sequence.fasta

    mkdir -p $pd_output

    [ ! -d ${pd_gcfid} ] && { echo "Warning: could not find directory for $gcfid."; continue; };
    [ ! -d ${pd_output} ] && { echo "Warning: Could not create output directory: $pd_output"; continue; };
    [ ! -f ${pf_sequence} ] && { echo "Warning: Could not find sequence file for $gcfid"; continue; };

    # create PBS file
    pf_pbs=$pd_output/run.pbs
    create_pbs_header $gcfid $pd_output > $pf_pbs
    echo "${bin_external}/gms2/gms2.pl --gcode $gencode --format gff --out gms2.gff --seq $pf_sequence  --v --genome-type $domain --fgio-dist-thresh 25;" >> $pf_pbs

    qsub -z $pf_pbs

done

