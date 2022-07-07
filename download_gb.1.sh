#!/bin/bash
while getopts ":i:o:a:" arg; do
		case "${arg}" in
				i)
						in_list=${OPTARG}
						;;
				o)
						out_dir=${OPTARG}
						;;
				a)
						index_in=${OPTARG}
						;;
		esac
done

mkdir -p ${out_dir} ; cd ${out_dir}

while IFS=$'\t' read -r -a my_array
do
	reference=${my_array[${index_in}]}
	echo $reference | awk 'BEGIN{OFS=FS="/"}{print $0,$NF"_genomic.gbff.gz"}' | xargs -n1 wget
done < ${in_list}
