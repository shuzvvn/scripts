#!/bin/bash
while getopts ":i:l:o:" arg; do
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
	esac
done

# time used
start=`date +%s`

# initialize outfile
cp ${in_file} ${out_file}

# read line to array
while IFS=$'\t' read -r -a my_array
do
	sed -i "s/${my_array[0]}/${my_array[3]}/g" ${out_file};
done < ${rename_table}

# print time used
end=`date +%s`
runtime=$((end-start))
echo "Time used ${runtime}s"
