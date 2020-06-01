##### Similarity-Based Start Prediction #####

### Step 0: Initialization

# ALWAYS source this configuration file to get the environment variables
# set up. This is especially important for running all the bash executables, since
# they assume that the environment variables have been initialized.
source config.sh        # set up env variables for bin, data, etc..
.
# If this is the first time running SBSP, set up the environment variables and 
# create executables for each driver program.
source install.sh       # create executables

### Step 1: Data download
# WARNING: This step uses PBS, takes a long time, and a lot of disk space.
# SBSP only comes with genomes that have verified genes. Download the remaining genomes
# This step downloads genomes (one per taxa id) under the given ancestors.
# The genomes are located in $data
# For each ancestor, the list of downloaded genomes will be located at: $lists/${database}_${ancestor}.list
# NOTE: If this loop fails (due to resource allocation), run each command separately, waiting until its done
# before running the next (hint: put an "echo" around the inner command to get all the download commands)
for database in refseq genbank; do
    for ancestor in Enterobacterales "Cyanobacteria/Melainabacteria group" Corynebacteriales Actinobacteria Alphaproteobacteria Archaea; do
        $bin/download_data_under_ancestor_sh.sh -a "Enterobacterales" -d $database
    done
done

# Download data for representative genomes
$bin/download_data_under_ancestor_sh.sh -a "Bacteria" -d refseq_representative
$bin/download_data_under_ancestor_sh.sh -a "Archaea" -d refseq_representative


### Step 2: Running SBSP
# Here are the individual commands to run all genomes with their corresponding ancestor using the 
# refseq database. Running on genbank requires replaced "refseq" with "genbank"
# These will also take a while to run and make use of PBS.
# The runs use sbsp_1.conf and pbs_1.conf configuration files, located in $config 
$bin/run_pipeline_sbsp_sh.sh  -q verified_ecoli -t genbank_enterobacterales -s 1 -p 1
$bin/run_pipeline_sbsp_sh.sh  -q verified_mtuberculosis -t genbank_corynebacteriales -s 1 -p 1
$bin/run_pipeline_sbsp_sh.sh  -q verified_hsalinarum -t genbank_archaea -s 1 -p 1
$bin/run_pipeline_sbsp_sh.sh  -q verified_npharaonis -t genbank_archaea -s 1 -p 1
$bin/run_pipeline_sbsp_sh.sh  -q verified_rdenitrificans -t genbank_alphaproteobacteria -s 1 -p 1
$bin/run_pipeline_sbsp_sh.sh  -q verified_synechocystis -t genbank_cyanobacteria_melainabacteria_group -s 1 -p 1

# You can run "$bin/run_pipeline_sbsp_sh.sh -h" for more information.

# Each run should create a directory located in $experiments. The directory name will have the form:
# q_verified_X_t_genbank_Y_ql_verified_tl_ncbi
# where X is the query genome and Y is the target ancestor
# In each directory is a subdirectory called "accuracy", and in it is a file called accuracy.csv which contains
# the accuracy of SBSP compared to the set of verified genes.
# Also in that directory are subdirectories msa_false and msa_true which contain the MSA for the true and false predictions

### Stats: Date of annotation for representative set
function get_hist_of_annotation_dates() {
    local pf_list="$1"

    echo "Date,Num Genomes"
    awk -F "," '{if (NR > 1) print $1}' $pf_list | while read -r gcfid; do
        curr_date=$(grep annotation-date $data/$gcfid/ncbi.gff | sed -E 's/^#[^ ]+\s([0-9]+)\/([0-9]+)\/([^ ]+).+/\1\/\3/g')

        echo "$curr_date"
    done | sort | uniq -c | awk '{print $2","$1}'
}

get_hist_of_annotation_dates $lists/refseq_representative_archaea.list

### Stats: GMS2 vs NCBI on representative set of genomes
$bin/run_gms2_on_list_sh.sh -l $lists/refseq_representative_bacteria.list --type bacteria
$bin/run_gms2_on_list_sh.sh -l $lists/refseq_representative_archaea.list --type archaea

function collect_gms2_vs_ncbi_stats() {
    local pf_list="$1"              # list of genomes
    local pf_output="$2"            # Path to output file 

    echo "GCFID,Group,GC,Found,Identical" > $pf_output
    dn_gms2=gms2

    # for each genome
    awk -F "," '{if (NR > 1) print $1}' $pf_list | while read -r gcfid; do

        pd_genome=$data/$gcfid

        group=$(grep group $pd_genome/runs/${dn_gms2}/GMS2.mod | cut -f 2 -d "-")        
        gc=$(probuild --stat --gc --seq $pd_genome/sequence.fasta | awk '{print $3}')

        v=$(compp -q -a $pd_genome/ncbi.gff -b $pd_genome/runs/${dn_gms2}/gms2.gff -v)

        found=$(echo "$v" | grep "found_in_A" | awk '{print $2}')
        ident=$(echo "$v" | grep "identical_in_A" | awk '{print $2}')

        echo -e "$gcfid,$group,$gc,$found,$ident"   
    done >> $pf_output
}

mkdir -p $exp/gms2_vs_ncbi_representative
for domain in bacteria archaea; do
    collect_gms2_vs_ncbi_stats $lists/refseq_representative_${domain}.list $exp/gms2_vs_ncbi_representative/stats_${domain}.csv
done




function get_2020_genomes() {
    local pf_list="$1"

    head -n 1 $pf_list
    awk -F "," '{if (NR > 1) print }' $pf_list | while read -r line ; do
        gcfid=$(echo "$line" | awk -F "," '{print $1}')
        year=$(grep annotation-date $data/$gcfid/ncbi.gff | sed -E 's/^#[^ ]+\s([0-9]+)\/([0-9]+)\/([^ ]+).+/\3/g')

        if [ "$year" == "2020" ]; then
            echo $line
        fi

        # echo "$curr_date"
    done 
}


#### Validate: On single-candidate genes

function collect_genes_with_ssc_from_list() {
    local pf_list="$1"
    local ancestor_name=$(basename $pf_list .list)

    local gcfid=ssc_${ancestor_name}
    local pd_gcfid=$data/$gcfid

    mkdir -p $pd_gcfid

    $bin/collect_genes_with_single_start_candidate_py.sh --pf-genome-list $pf_list --pd-output $pd_gcfid

    # create list
    pf_gcfid_list=$lists/{$gcfid}.list
    echo "gcfid,genetic-code,attributes" > $pf_gcfid_list
    echo "$gcfid,11,." >> $pf_gcfid_list

    echo $gcfid
}

# Example: corynebacteriales

target_id=genbank_corynebacteriales
pf_list=$lists/${target_id}.list

gcfid=$(collect_genes_with_ssc_from_list $pf_list)
pf_gcfid_list=$lists/${gcfid}.list

$bin/run_pipeline_sbsp_sh.sh  -q $gcfid -t $target_id -s 1 -p 1


##### Gathering results
















