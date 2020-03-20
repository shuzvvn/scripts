#!/bin/bash
mkdir -p ${out_dir}

# get regionX
read -p "contigID:start-end: " regionX
echo Getting region.sam from .bam ...
samtools view ${bam_file} "${regionX}" | samtools view -Shu -t ${ref_fasta}.fai - | samtools sort -m 30000000000 -o ${out_dir}regionX.bam ;
echo Getting region.id from region.sam ...
samtools view ${out_dir}regionX.bam | cut -f 1 | sort -u > ${out_dir}regionX.id ;
echo Extracting raw reads by region.id ...

# extract reads
perl /home/chkuo/plscript/filter_fastq_by_id.1.pl --list_file=${out_dir}regionX.id --regex_list="^(\S+)" --in_file=${paired_fastq} --regex_in="^@(\S+)" --out_file=${out_dir}regionX.fastq --mode="include" --verbose=1

# covert fastq to fastq and qual
python3 /home/shutingcho/pyscript/convert_fastq_to_fasta.1.py --in_fastq=${out_dir}regionX.fastq --out_fasta=${out_dir}regionX.fasta --out_qual=${out_dir}regionX.qual

echo Running phrap ...
phrap ${out_dir}regionX.fasta &> ${out_dir}regionX.phrapout ;
echo Renaming contigs from phrap results...
perl /home/chkuo/plscript/seq_rename_regex.6.pl --in_file=${out_dir}regionX.fasta.contigs --out_file=${out_dir}regionX.fasta.contigs.rename --regex="\/regionX.fasta.(Contig\d+)";
cat ${out_dir}regionX.fasta.contigs.rename
