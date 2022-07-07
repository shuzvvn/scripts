#!/bin/bash
while getopts ":i:g:o:" arg; do
	case "${arg}" in
		i)
			in_file=${OPTARG}
			;;
		g)
			gi_taxid=${OPTARG}
			;;
		o)
			out_file=${OPTARG}
			;;
	esac
done

# check if out_file already exist => delete the file
if [ -e ${out_file} ]; then
	rm ${out_file}
fi

# time used
start=`date +%s`

# processing in_file line by line
while IFS=$'\t' read -r -a my_array # read line to array
do
taxid=`grep --max-count=1 "^${my_array[2]}"$'\t' ${gi_taxid} | cut -f2` # find the taxon id of the gi (=${my_array[2]})
echo "${my_array[*]}"$'\t'"${taxid}" >> ${out_file} # add the taxon id to the end of the line, write to out_file
done < ${in_file}

# print time used
end=`date +%s`
runtime=$((end-start))
echo "Time used ${runtime}s"