
function setup_busco_dataset() {
    # Download and setup a dataset from busco
    # E.g. enterobacteriales 9 data lists
    local name=$1
    local odbversion=$2
    local pd_data=$3
    local pf_list=$4
    local w=$5
    
    name_and_version="${name}_odb${odbversion}"
    dlname=${name_and_version}.tar.gz
    
    # download data (if doesn't exist)
    if [ ! -f $w/$dlname ]; then 
    	wget http://busco.ezlab.org/v2/datasets/$dlname -P $w
    fi
    
    # unzip 
    tar -zxf $w/$dlname -C $w
    
    pf_assembly=$metadata/assembly_summary.txt
    fn_filtered_assembly=$w/filtered_assembly_summary.txt
    pf_species_info=$w/${name_and_version}/info/species.info
    
    # extract species from assembly
    $bin/extract-from-assembly-summary_py.sh --summary ${pf_assembly} --out $fn_filtered_assembly --key taxid --pf-key-values ${pf_species_info}  --full-genome
    
    # download data
    download_genomic_data_from_assembly_summary ${fn_filtered_assembly} ${pd_data}
    
    # make list
    create_list_from_assembly_summary ${fn_filtered_assembly} ${pf_list}
    
    # fix fasta header label in sequence files
    set_fasta_header_with_accession ${pd_data} ${pf_list}
    
    # extract ribosomal labels
    extract_ribosomal_labels ${pd_data} ${pf_list} ncbi.gff ncbi_ribosomal.gff
}




