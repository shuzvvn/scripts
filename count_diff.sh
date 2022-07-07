#!/bin/bash
# START of shell script to run predictnls
while getopts ":i:j:" arg; do
	case "${arg}" in
		i)
			inlist1=${OPTARG}
			;;
		j)
			inlist2=${OPTARG}
			;;
	esac
done

# check if item in list1 also in list2
list1=`cat ${inlist1} | tr '\n' ' '`
count_in=0
count_out=0
for i in ${list1}; do
count_in=$((count_in+1))
if grep -Fq "${i}" ${inlist2} ; then
count_out=$((count_out+1))
fi
done
# print report
echo count_in = ${count_in}, count_out = ${count_out} 
# END of shell script
