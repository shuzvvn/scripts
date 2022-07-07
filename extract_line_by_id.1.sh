#!/bin/bash
while getopts ":i:l:o:a:b:" arg; do
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
			index_out=${OPTARG}
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
	awk 'BEGIN {FS="\t"} $"'"${index_out}"'"=="'"${my_array[${index_in}]}"'"' ${in_file} >> ${out_file}
done < ${in_list}