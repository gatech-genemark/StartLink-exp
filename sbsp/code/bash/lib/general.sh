#!/usr/bin/env bash

function get_column_names_indeces() {
	fname="$1"

	head -n 1 $fname | awk -F "," '{ for (i = 1; i <= NF; i++) {print i"\t"$i}}'
}


function print_alignments_in_hr_format() {
	fname="$1"
	filter="$2"

	cat $fname | grep "${filter}$" |  awk -F "," '{print $29"\n"$30"\n"$31"\n"$33"\n"$34"\n"$32"\n"}'
}