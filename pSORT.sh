#!/bin/bash
# START of shell script to Run pSORT for NLS prediction
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
# run pSORT on webserver http://psort1.hgc.jp/cgi-bin/okumura.pl
for locus_tag in $(echo `ls ${inDir}` | sed 's/.fasta / /g' | sed 's/.fasta//g')
do
	seq=`tail -1 ${inDir}/${locus_tag}.fasta`
	curl "http://psort1.hgc.jp/cgi-bin/okumura.pl?origin=plant&title=${locus_tag}&sequence=${seq}" -H 'Accept-Encoding: gzip, deflate, sdch' -H 'Accept-Language: en-US,en;q=0.8,zh-TW;q=0.6,zh;q=0.4' -H 'Upgrade-Insecure-Requests: 1' -H 'User-Agent: Mozilla/5.0 (Macintosh; Intel Mac OS X 10_11_2) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/48.0.2564.116 Safari/537.36' -H 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8' -H 'Referer: http://psort1.hgc.jp/form.html' -H 'Connection: keep-alive' -H 'AlexaToolbar-ALX_NS_PH: AlexaToolbar/alxg-3.3' --compressed > ${locus_tag}.txt
done
# nucleus score
grep 'nucleus --- Certainty=' *.txt | sed 's/.txt:                          nucleus --- Certainty= /\t/g' | sed 's/(Affirmative) < succ>//g' > pSORT.nuc &
# chloroplast thylakoid membrane score
grep 'chloroplast thylakoid membrane --- Certainty=' *.txt | sed 's/.txt:   chloroplast thylakoid membrane --- Certainty= /\t/g' | sed 's/(Affirmative) < succ>//g' > pSORT.chltm &
# chloroplast stroma score
grep 'chloroplast stroma --- Certainty= ' *.txt | sed 's/.txt:               chloroplast stroma --- Certainty= /\t/g' | sed 's/(Affirmative) < succ>//g' > pSORT.chlstr &
# Nuclear Signal
grep 'Nuclear Signal' *.txt | sed 's/.txt:/\t/g' | sed 's/Nuclear Signal   Status: //g' > pSORT1
grep -L 'Nuclear Signal' *.txt | sed 's/.txt/\tNA/g' > pSORT0
sort pSORT1 pSORT0 > pSORT
rm pSORT1 pSORT0
# Final Results
grep -A2 'Final Results' *.txt | grep 'Affirmative' | sed 's/.txt- */\t/g' | sed 's/^ *//g' | sed 's/ --- Certainty= / /g' | sed 's/(Affirmative) < succ>//g' > pSORT.result1
grep -L 'Final Results' *.txt | sed 's/.txt/\tNA/g' > pSORT.result0
sort pSORT.result1 pSORT.result0 > pSORT.result
rm pSORT.result1 pSORT.result0
join -t $'\t' pSORT pSORT.result > pSORT_sum
rm pSORT pSORT.result
# print report
inSeqNo=`ls ${inDir} | grep -c 'fasta'` 
nlsNo=`awk '$2~/positive/ {count++} END {print count}' pSORT_sum`
echo count_in = ${inSeqNo}, count_out = ${nlsNo} 
# END of shell script to run pSORT
