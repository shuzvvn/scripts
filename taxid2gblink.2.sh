#!/bin/bash
while getopts ":i:l:o:a:g:" arg; do
		case "${arg}" in
				i)
						in_file=${OPTARG}
						;;
				l)
						in_list=${OPTARG}
						;;
				o)
						out_file=${OPTARG}
						;;
				a)
						index_in=${OPTARG}
						;;
				b)
						complete_genome=${OPTARG}
						;;
		esac
done

# initialize outfile
if [ -e ${out_file} ]; then
	rm ${out_file}
fi;

# read line to array

while IFS=$'\t' read -r -a my_array
do
	taxid=${my_array[${index_in}]}
	if [ ${complete_genome} ]; then
	awk -v TAXID="$taxid" 'BEGIN{FS=OFS="\t"}{if($6==TAXID && $12=="Complete Genome" && $14=="Full") {print $5,$6,$8,$9,$15,$20}}' ${in_list}
	else
	awk -v TAXID="$taxid" 'BEGIN{FS=OFS="\t"}{if($6==TAXID) {print $5,$6,$8,$9,$12,$14,$15,$20}}' ${in_list}
fi
done < ${in_file} > ${out_file}
