#!/bin/bash
# START of shell script to run predictnls
while getopts ":i:o:" arg; do
	case "${arg}" in
		i)
			inDir=${OPTARG}
			;;
		o)
			outDir=${OPTARG}
			;;
	esac
done
# prepare directory
mkdir -p ${outDir}; cd ${outDir}
# run PredictNLS
for locus_tag in $(echo `ls ${inDir}` | sed 's/.fasta / /g' | sed 's/.fasta//g')
do
	predictnls fileIn=${inDir}/${locus_tag}.fasta \
	fileOut=${outDir}/${locus_tag}.nls \
	fileSummary=${outDir}/${locus_tag}.sum;
	grep -A3 "List of NLS's found in sequence" ${locus_tag}.nls | tail -1 | sed 's/^| *//g' | sed 's/|$//g' | sed 's/| */\t/g' > ${locus_tag}.nls.1
done
# check result
grep '[0-1]' *.sum | sed 's/.sum:/\t/g' > predictnls_sum
grep '.' *.nls.1 | sed 's/.nls.1:/\t/g' > predictnls_seq
# clean up
rm *.nls*
rm *.sum
# print report
inSeqNo=`ls ${inDir} | grep -c 'fasta'`
nonnls=`awk '$2==0 {count++} END {print count}' predictnls_sum`
nlsNo=$((inSeqNo-nonnls))
echo count_in = ${inSeqNo}, count_out = ${nlsNo}
# END of shell script to run predictnls
