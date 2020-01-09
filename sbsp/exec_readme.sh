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
for database in refseq genbank; do
    for ancestor in Enterobacterales "Cyanobacteria/Melainabacteria group" Corynebacteriales Actinobacteria Alphaproteobacteria Archaea; do
        $bin/download_data_under_ancestor_sh.sh -a "Enterobacterales" -d $database
    done
done

# Download data for representative genomes
$bin/download_data_under_ancestor_sh.sh -a "root" -d refseq_representative


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

