#!/bin/bash
while getopts ":i:l:o:a:b:" arg; do
	case "${arg}" in
		i)
			in_file=${OPTARG}
			;;
		l)
			rename_table=${OPTARG}
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
cp ${in_file} ${out_file}

# read line to array
while IFS=$'\t' read -r -a my_array
do
	sed -i "s/${my_array[${index_in}]}/${my_array[${index_out}]}/g" ${out_file};
done < ${rename_table}