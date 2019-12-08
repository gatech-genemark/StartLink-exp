

function set_fasta_header_with_accession() {
	local datapath="$1";
	local genomeList="$2";
	awk '{print $1}' $genomeList | while read -r line; do
		sed -E -i   's/^(>[^ ]+).*/\1/' $datapath/$line/sequence.fasta 
	done
}
