


function construct_training_set_loop_propagation() {

    pf_q_list="$1"
    pf_t_list="$2"
    fn_q_labels="$3"
    fn_t_labels="$4"
    pf_output="$5"
    pd_work="$6"
    pd_data="$7"
    pd_bin="$8"
    pbs_opts=""

    if [[ "$#" -eq 9 ]]; then
        pbs_opts="$9"
    fi

    if [[ "$#" -lt 8 ]]; then
        echo "Incorrect number of arguments"
        return
    fi

    q_genome_name=$(head -n 1 $pf_q_list | awk '{print $1}')

    # Construct samples from verified
    pf_verified_final_data=$pd_work/verified_final_data.tab
    pf_verified_inter_data=$pd_work/verified_inter_data.tab

    pf_verified_final_data_labeled=$pd_work/verified_final_data_labeled.tab
    pf_verified_inter_data_labeled=$pd_work/verified_inter_data_labeled.tab

    pf_verified_final_data_labeled_with_features=$pd_work/verified_final_data_labeled_with_features.tab
    pf_verified_inter_data_labeled_with_features=$pd_work/verified_inter_data_labeled_with_features.tab

    # get orthologs
    
    $bin/loop-propagate_py.sh --pf-q-list $pf_q_list --pf-t-list $pf_t_list --fn-q-labels $fn_q_labels --fn-t-labels $fn_t_labels --prefix "v_" --pf-final-data $pf_verified_final_data --pf-intermediate-data $pf_verified_inter_data --path-work $pd_work --path-data $pd_data ${pbs_opts}

    # label
    $bin/generate-labeled-data-from-alignments-file_py.sh --propagation $pf_verified_final_data --out $pf_verified_final_data_labeled from-labels --labels $pd_data/$q_genome_name/verified.gff --gt-source target

    # add features    
    $bin/add-features-to-alignments_py.sh --data $pf_verified_final_data_labeled --p-dir-data $pd_data --out $pf_verified_final_data_labeled_with_features --num-processors 8 -l INFO --write-partial ${pbs_opts} --pd-work $pd_work

    cp $pf_verified_final_data_labeled_with_features $pf_output
}


function construct_training_set_negative() {

    pf_q_list="$1"
    pf_t_list="$2"
    fn_q_labels_begin="$3"
    fn_t_labels_begin="$4"
    pf_output="$5"
    search_method="$6"
    offset="$7"
    region_length="$8"
    pd_work="$9"
    pd_data="${10}"
    pd_bin="${11}"

    ignore_stops=""

    if [[ "$#" -gt 11 ]]; then
        ignore_stops="${12}"
    fi

    max_sample_size="1"

    if [[ "$#" -gt 12 ]]; then
        max_sample_size="${13}"
    fi

    pbs_opts=""

    if [[ "$#" -eq 14 ]]; then
        pbs_opts="${14}"
    fi
    
    q_genome_name=$(head -n 1 $pf_q_list | awk '{print $1}')
    pd_q_genome=$pd_data/$q_genome_name
    pf_q_sequence=$pd_q_genome/sequence.fasta
    pf_q_verified=$pd_q_genome/${fn_q_labels_begin}

    echo "${q_genome_name}"


    prefix=${search_method}_of_${offset}_${region_length}_

    tag=$(echo "${search_method}" | cut -b1)


    # Construct samples from downstream of verified
    fn_q_labels=${prefix}_${fn_q_labels_begin}
    fn_t_labels=${fn_t_labels_begin}

    pf_final_data=$pd_work/${prefix}verified_final_data.tab
    pf_inter_data=$pd_work/${prefix}verified_inter_data.tab

    pf_final_data_labeled=$pd_work/${prefix}verified_final_data_labeled.tab
    pf_inter_data_labeled=$pd_work/${prefix}verified_inter_data_labeled.tab

    pf_final_data_labeled_with_features=$pd_work/${prefix}verified_final_data_labeled_with_features.tab
    pf_inter_data_labeled_with_features=$pd_work/${prefix}verified_inter_data_labeled_with_features.tab

    echo "$bin/get-candidate-starts_py.sh  --p-sequences ${pf_q_sequence} --p-labels ${pf_q_verified} --search-method ${search_method} --offset $offset --region-length ${region_length} --max-sample-size ${max_sample_size} ${ignore_stops} > ${pd_q_genome}/${fn_q_labels}"
    $bin/get-candidate-starts_py.sh  --p-sequences ${pf_q_sequence} --p-labels ${pf_q_verified} --search-method ${search_method} --offset $offset --region-length ${region_length} --max-sample-size ${max_sample_size} ${ignore_stops} > $pd_q_genome/${fn_q_labels}

    echo "$bin/loop-propagate_py.sh --pf-q-list $pf_q_list --pf-t-list $pf_t_list --fn-q-labels $fn_q_labels --fn-t-labels $fn_t_labels --prefix \"${tag}_\" --pf-final-data $pf_final_data --pf-intermediate-data $pf_inter_data --path-work $pd_work --path-data $pd_data --no-loop --verbose ${pbs_opts}"
    $bin/loop-propagate_py.sh --pf-q-list $pf_q_list --pf-t-list $pf_t_list --fn-q-labels $fn_q_labels --fn-t-labels $fn_t_labels --prefix "${tag}_" --pf-final-data $pf_final_data --pf-intermediate-data $pf_inter_data --path-work $pd_work --path-data $pd_data --no-loop --verbose ${pbs_opts}

    
    $bin/generate-labeled-data-from-alignments-file_py.sh --propagation $pf_inter_data --exclude-positively-labeled --out $pf_inter_data_labeled from-labels --labels $pf_q_verified --gt-source query 


    $bin/add-features-to-alignments_py.sh --data $pf_inter_data_labeled --p-dir-data $pd_data --out $pf_inter_data_labeled_with_features --num-processors 8 -l INFO --write-partial ${pbs_opts} --pd-work $pd_work 

    cp $pf_inter_data_labeled_with_features $pf_output
}


function construct_testing_set() {

    pf_q_list="$1"
    pf_t_list="$2"
    fn_q_labels="$3"
    fn_t_labels="$4"
    pf_output="$5"
    pd_work="$6"
    pd_data="$7"
    pd_bin="$8"

    pbs_opts=""

    if [[ "$#" -eq 9 ]]; then
        pbs_opts="$9"
    fi


    if [[ "$#" -ne 8 ]]; then
        echo "Illegal number of parameters"
    fi

    if [ -f $pf_output ] ; then
        rm $pf_output
    fi
    
    while read -r q_line; do
        
        q_genome_name=$(echo "$q_line" | awk '{print $1}')
        pd_q_genome=$pd_data/${q_genome_name}
        
        while read -r t_line; do
            t_genome_name=$(echo "$t_line" | awk '{print $1}')
            
            echo "$t_genome_name"
            
            pref=${q_genome_name}_and_${t_genome_name}
            
            pf_alignments=$pd_work/${pref}_alignments
            
            # get orthologs and locally align starts
            $bin/locally-align-start-region_py.sh --q-name ${q_genome_name} --t-name ${t_genome_name} --fn-q-labels ${fn_q_labels} --fn-t-labels ncbi.gff --pf-out $pf_alignments --pd-data $pd_data --pd-work $pd_work

            if [[ -f $pf_output ]]; then
                grep -v "strand" $pf_alignments >> $pf_output
            else
                cp $pf_alignments $pf_output
            fi

            rm $pf_alignments
            
        done < "$pf_t_list"
        
    done < "$pf_q_list"


}




