#!/bin/bash
while getopts ":i:g:a:d:o:" arg; do
	case "${arg}" in
		i)
			in_file=${OPTARG}
			;;
		g)
			gi_taxid=${OPTARG}
			;;
		a)
			names=${OPTARG}
			;;
		d)
			nodes=${OPTARG}
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

# Obtain the name corresponding to a taxid or the taxid of the parent taxa
get_name_or_taxid()
{
	grep --max-count=1 "^${1}"$'\t' "${2}" | cut --fields="${3}"
}

while IFS=$'\t' read -r -a my_array # read line to array
do
	gi=${my_array[2]} # read gi
	taxonomy="" # initialize taxonomy as NA

	# Get the taxid corresponding to the GI number
	taxid=$(get_name_or_taxid "${gi}" "${gi_taxid}" "2")
	init_taxid=${taxid}

	# Loop until you reach the root of the taxonomy (i.e. taxid = 1)
	if [[ ${taxid} =~ ^[0-9]+$ ]]; then
		while [ "${taxid}" -gt 1 ]; do
			# Obtain the scientific name corresponding to a taxid
			name=$(get_name_or_taxid "${taxid}" "${names}" "3")
			# Obtain the parent taxa taxid
			parent=$(get_name_or_taxid "${taxid}" "${nodes}" "3")
			# Build the taxonomy path
			taxonomy="${name}\t${taxonomy}"
			taxid="${parent}"
		done
	else
		init_taxid='NA'
		taxonomy='NA'
	fi
	( IFS=$'\t'; echo -e "${my_array[*]}\t${init_taxid}\t${taxonomy}") >> ${out_file} # add the taxon id to the end of the line, write to out_file
done < ${in_file}

# print time used
end=`date +%s`
runtime=$((end-start))
echo "Time used ${runtime}s"