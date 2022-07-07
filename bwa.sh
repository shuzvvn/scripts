#!/bin/bash
# set default for AS
if [ "${run_PE_AS}" == "true" ]; then
	run_PE_AS=200
fi
if [ "${run_flash_AS}" == "true" ]; then
	run_flash_AS=420
fi
thread="${thread:-32}"

# bwa index
bwa index -a is ${fasta_file} ;

# index reference sequence in the FASTA format
samtools faidx ${fasta_file} ;

# bwa mem
mkdir -p ${out_dir} ;
cd ${out_dir} ;

if ${run_PE} ; then
	# report job start time
	now=$(date +"%r")
	echo " bwa for ${project_ID} PE start : $now"

	# map PE reads
	bwa mem -t ${thread} ${fasta_file} ${PE1} ${PE2} 2> ${out_dir}/PE.log | samtools view -Su - | samtools sort -m 30000000000 - ${out_dir}/PE ;

	# report job finished time
	now=$(date +"%r")
	echo " bwa for ${project_ID} PE end : $now"
	
	# filter by alignment score 30
	samtools view ${out_dir}/PE.bam | perl -ne '/AS\:i\:(\d+)/; print if $1 >= 30' | samtools view -bt ${fasta_file}.fai - > ${out_dir}/PE_AS30.bam ;
	samtools index ${out_dir}/PE_AS30.bam ;
	# clean up
	rm ${out_dir}/PE.bam

	# filter by setting alignment score
	if [ ${run_PE_AS} ] ; then
	samtools view ${out_dir}/PE_AS30.bam | perl -ne '/AS\:i\:(\d+)/; print if $1 >= '"${run_PE_AS}"'' | samtools view -bt ${fasta_file}.fai - > ${out_dir}/PE_AS${run_PE_AS}.bam ;
	samtools index ${out_dir}/PE_AS${run_PE_AS}.bam ;

	# report job finished time
	now=$(date +"%r")
	echo " AS${run_PE_AS} filter for ${project_ID} PE end : $now"
	fi
fi

if ${run_flash} ; then
	# report job start time
	now=$(date +"%r")
	echo " bwa for ${project_ID} flash start : $now"
	
	# map flash reads
	bwa mem -p -t ${thread} ${fasta_file} ${flash_file} 2> ${out_dir}/flash.log | samtools view -Su - | samtools sort -m 30000000000 - ${out_dir}/flash ;
	samtools index ${out_dir}/flash.bam ; 

	# report job finished time
	now=$(date +"%r")
	echo " bwa for ${project_ID} flash end : $now"

	# filter by alignment score 
	if [ ${run_flash_AS} ] ; then
	samtools view ${out_dir}/flash.bam | perl -ne '/AS\:i\:(\d+)/; print if $1 >= '"${run_flash_AS}"'' | samtools view -bt ${fasta_file}.fai - > ${out_dir}/flash_AS${run_flash_AS}.bam ;
	samtools index ${out_dir}/flash_AS${run_flash_AS}.bam ;

	# report job finished time
	now=$(date +"%r")
	echo " AS filter for ${project_ID} flash end : $now"
	fi
fi
### end of script ###
