


pf_list=$lists/genbank_corynebacteriales.list
pf_list=$lists/genbank_actinobacteria.list
#pf_list=$lists/verified.list

function get_up_nodes() {
    pbsnodes -l  up | awk '{print $1}'
}

function get_numbered_nodes() {
    get_up_nodes | nl -s , -n rn -w 1
}


# run command
function create_pbs_header() {
    local node="$1"
    local pf_list="$2"
    local pd_download="/scratch/karl/sbsp_data"

    rm -f output_${node}
    
    echo "\
#PBS -N $node
#PBS -o output_${node}
#PBS -j oe
#PBS -l nodes=$node
#PBS -l walltime=07:00:00:00
#PBS -V
#PBS -W umask=002

# rm -rf /scratch/karl/sbsp_data/*
mkdir -p $pd_download

PBS_O_WORKDIR=${pd_download}
cd \$PBS_O_WORKDIR


time { pf_list=$pf_list; total=\$(awk 'END{print NR - 1}' \$pf_list); counter=1; awk -F \",\" '{if (NR > 1) print \$1}' \$pf_list | while read -r gcfid; do echo -ne \"\$counter / \$total\\r\"; counter=\$((counter + 1)); rsync   -r --exclude=runs  $data/\$gcfid  /scratch/karl/sbsp_data; done }
"
}

pf_job_ids=job_id.txt
rm -r $pf_job_ids

get_up_nodes | while read -r node; do 
    create_pbs_header $node $pf_list > ${node}.pbs
    qsub ${node}.pbs >> $pf_job_ids
done

# wait until all jobs are done

job_ids=$(cat $pf_job_ids)

echo "$job_ids" | while read -r jobname; do

    while [ $(qstat -u karl -a | grep " R\|Q\|H " | grep "$jobname"  | wc -l) != 0 ]; do 
        sleep 30
    done
done


