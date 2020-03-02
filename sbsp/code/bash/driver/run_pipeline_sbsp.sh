#!/bin/bash

# Option default
VERBOSE=0
DEBUG=2
INFO=1
fn_q_labels=verified.gff
fn_q_labels_true=verified.gff
fn_t_labels=ncbi.gff
pbs_conf=defaults
tag=
steps=""
email_on_complete=

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
-q; Name of list containing query genomes (without .list suffix)
-t; Name of list containing target genomes (without .list suffix)
-s, --sbsp-conf;  Configuration number of SBSP options (found in file sbsp_NUM.conf)
-p, --pbs-conf;  Configuration number of PBS options (found in pbs_NUM.conf) (default: ${pbs_conf})
--q-labels; Name of labels file for query genomes (default: ${fn_q_labels})
--t-labels; Name of labels file for target genomes (default: ${fn_t_labels})
-e; Sends email when job completes
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
                -q)                 query=$2; shift;;
                -t)                 target=$2; shift;;
                -s|--sbsp-conf)     sbsp_conf=$2; shift;;
                -p|--pbs-conf)      pbs_conf=$2; shift;;
                --steps)            steps="--steps $2"; shift;;
                --q-labels-true)    fn_q_labels_true=$2; shift;;
                --q-labels)         fn_q_labels=$2; shift;;
                --t-labels)         fn_t_labels=$2; shift;;
                --tag)              tag="$2"; shift;;
                -v|--verbose)       VERBOSE=$(($VERBOSE + 1));;
                -e)                 email_on_complete=1;;
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
([ -z "$query" ] || [ -z "$target" ] || [ -z "${sbsp_conf}" ]) && usage

pf_q_list=$lists/${query}.list;
pf_t_list=/scratch/karl/sbsp_data/${target}.dmnd;
pf_sbsp_conf=$config/sbsp_${sbsp_conf}.conf
pf_pbs_conf=$config/pbs_${pbs_conf}.conf

# check that files exist
[ ! -f "${pf_q_list}" ] && { echo "File doesn't exist: $pf_q_list"; exit 1; };
[ ! -f "${pf_sbsp_conf}" ] && { echo "File doesn't exist: $pf_sbsp_conf"; exit 1; };
[ ! -f "${pf_pbs_conf}" ] && { echo "File doesn't exist: $pf_pbs_conf"; exit 1; };

# run command
fn_q_labels_no_ext="${fn_q_labels%.*}"
fn_t_labels_no_ext="${fn_t_labels%.*}"

if [ ! -z "$tag" ]; then
    tag="${tag}_"
fi

pd_run=$exp/${tag}q_${query}_t_${target}_sbsp_${sbsp_conf}_ql_${fn_q_labels_no_ext}_tl_${fn_t_labels_no_ext}
dn_compute=pbs_${tag}q_${query}_t_${target}_sbsp_${sbsp_conf}_ql_${fn_q_labels_no_ext}_tl_${fn_t_labels_no_ext}
pf_output=${pd_run}/output.csv

mkdir -p $pd_run

if verbose $INFO; then
    echo "Working directory: $pd_run"
fi

$bin/pipeline_msa_py.sh --pf-q-list $pf_q_list --pf-t-db $pf_t_list --fn-q-labels $fn_q_labels --pf-sbsp-options $pf_sbsp_conf --pf-pbs-options $pf_pbs_conf --pf-output $pf_output --pd-work $pd_run --fn-q-labels-true $fn_q_labels_true --dn-compute ${dn_compute} $steps

if [ ! -z ${email_on_complete} ]; then
    mail -s "Done" <<< "Job at ${pd_run}"  2> /dev/null
fi

