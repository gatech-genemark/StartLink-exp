#!/usr/bin/env bash

# Karl Gemayel
# Georgia Institute of Technology

function extract_ribosomal_labels() {
	local datapath="$1";
	local genomeList="$2";
	local nameannot="$3";
	local nameout="$4";
	awk '{print $1}' $genomeList | while read -r line; do
	echo "$line"
		grep -E "\s+ribosomal\s+" $datapath/$line/$nameannot > $datapath/$line/$nameout
	done
}