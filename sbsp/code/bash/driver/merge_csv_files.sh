#!/bin/bash

# Option default
VERBOSE=0
DEBUG=2
INFO=1
list_pf_csv=
pf_output=

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
--in; List of CSV files
--pf-output; Output file
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
                --in)               list_pf_csv="${list_pf_csv} $2"; while [ "$#" -ge 3 ] && [ "$3" != "--pf-output" ]; do shift; list_pf_csv="${list_pf_csv} $2"; done; shift;;
                --pf-output)        pf_output="$2"; shift;;
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
([ -z "$list_pf_csv" ] || [ -z "$pf_output" ]) && usage

list_pf_csv=$(echo "$list_pf_csv" | while read -r line; do ls $line; done)

first="yes"

echo "$list_pf_csv" | while read -r pf_input; do
    # if file empty or doesn't exist
    if [ ! -f $pf_input ] || [ ! -s $pf_input ]; then
        continue
    fi

    if [ "$first" == "yes" ]; then
        cp $pf_input $pf_output
        first="no"
    else
        awk '{if (NR > 1) print}' "$pf_input" >> $pf_output
    fi
done








