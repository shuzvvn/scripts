#!/bin/bash
while getopts ":i:l:o:a:" arg; do
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
	awk -v TAXID="$taxid" 'BEGIN{FS=OFS="\t"}{if($6==TAXID && $12=="Complete Genome" && $14=="Full") {print $5,$6,$8,$9,$15,$20}}' ${in_list} >> ${out_file}
done < ${in_file}
