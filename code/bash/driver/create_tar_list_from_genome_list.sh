#!/bin/bash

# Option default
VERBOSE=0
DEBUG=2
INFO=1
pf_genome_list=
pf_tar=
email_on_complete=

# getopts string
# This string needs to be updated with the single character options (e.g. -f)
opts="v"

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
`cmd` --pf-genome-list PF --pf-tr PF [OPTION...]
--pf-genome-list; Genome list file
--pf-tar; Output tar file
-e; Sends email when job completes
-v, --verbose; Enable verbose output (include multiple times for more
             ; verbosity, e.g. -vvv)
" | column -t -s ";"
exit 1; 
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
                --pf-tar)           pf_tar=$2; shift;;
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
([ -z "${pf_genome_list}" ] || [ -z "${pf_tar}" ]) && usage


awk -F "," '{if (NR > 1) print $1}' $pf_genome_list | while read -r gcfid; do echo -e "$gcfid/sequence.fasta\n$gcfid/ncbi.gff" ; done > x.tmp.list

FROMSIZE=`wc -l x.tmp.list | awk '{print $1}'`
echo "$FROMSIZE"
CHECKPOINT=`echo "$FROMSIZE*40" | bc`;
echo "Estimated: [==================================================]";
echo -n "Progess:   [";
tar  --record-size=1024 --checkpoint="${CHECKPOINT}" --checkpoint-action="ttyout=>" -c  -f $pf_tar --directory=$data -T x.tmp.list 
echo "]"

rm x.tmp.list
