#!/bin/bash

# Option default
VERBOSE=0
DEBUG=2
INFO=1
fn_q_labels=ncbi.gff
fn_q_labels_true=ncbi.gff
fn_t_labels=ncbi.gff
pbs_conf=defaults
steps=""
tag=
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
                --pf-genome-list)   pf_genome_list=$2; shift;;
                -s|--sbsp-conf)     sbsp_conf=$2; shift;;
                -p|--pbs-conf)      pbs_conf=$2; shift;;
                --steps)            steps="--steps $2"; shift;;
                --q-labels-true)    fn_q_labels_true=$2; shift;;
                --q-labels)         fn_q_labels=$2; shift;;
                --t-labels)         fn_t_labels=$2; shift;;
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
([ -z "$pf_genome_list" ] || [ -z "${sbsp_conf}" ]) && usage

pf_sbsp_conf=$config/sbsp_${sbsp_conf}.conf
pf_pbs_conf=$config/pbs_${pbs_conf}.conf

# check that files exist
[ ! -f "${pf_genome_list}" ] && { echo "File doesn't exist: $pf_genome_list"; exit 1; };
[ ! -f "${pf_sbsp_conf}" ] && { echo "File doesn't exist: $pf_sbsp_conf"; exit 1; };
[ ! -f "${pf_pbs_conf}" ] && { echo "File doesn't exist: $pf_pbs_conf"; exit 1; };

# run command
fn_q_labels_no_ext="${fn_q_labels%.*}"
fn_t_labels_no_ext="${fn_t_labels%.*}"


function parse_entry_for_gcfid() {
    echo "$1" | awk -F "," '{print $1}'
}

function parse_entry_for_scientific_name() {
    echo "$1" | sed -E 's/^.*name=([^;]+).*$/\1/g'
}

function parse_entry_for_ancestor_name() {
    echo "$1" | sed -E 's/^.*ancestor=([^;]+).*$/\1/g'
}

function convert_string_to_valid_directory_name() {
    line="$1"
    echo "$line" | tr '[:upper:]' '[:lower:]' | tr '/' "_" | tr " " "_" | tr "-" "_"
}

function run_sbsp_on_genome_list_entry() {
    local entry="$1"
    local gcfid=$(parse_entry_for_gcfid "$entry")
    local ancestor=$(parse_entry_for_ancestor_name "$entry")

    # create temporary list for gcfid
    local pf_list=$lists/tmp_${gcfid}.list

    echo "gcfid,genetic-code,attributes" > $pf_list
    echo "$entry" >> $pf_list
    
    local ancestor_valid=$(convert_string_to_valid_directory_name "$ancestor")
    local tag_option=""
    if [ ! -z "$tag" ]; then
        tag_option="--tag $tag"
    fi

    $bin/run_pipeline_sbsp_sh.sh  -q tmp_$gcfid -t genbank_${ancestor_valid}  -s ${sbsp_conf} -p ${pbs_conf} --q-labels ncbi.gff --q-labels-true ${fn_q_labels_true} ${tag_option} $steps
    rm $pf_list
}


function run_sbsp_on_genome_list() {
    local pf_genome_list="$1"
    [ "$#" -ne 1 ] && { echo "Missing argument 'pf_genome_list'"; return ; }
    
    local total=$(wc -l $pf_genome_list | awk '{print $1-1}')
    local current=1
    awk -F "," '{if (NR > 1) print }' $pf_genome_list | while read -r line; do
        local gcfid=$(parse_entry_for_gcfid "$line")
        echo -e "Progress: $current/$total. GCFID=$gcfid"
        
        run_sbsp_on_genome_list_entry "$line" &
        sleep 5
        current=$((current + 1))
    done
}

run_sbsp_on_genome_list $pf_genome_list
