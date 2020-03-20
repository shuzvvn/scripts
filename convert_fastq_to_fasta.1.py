#!/usr/bin/python3

# convert_fastq_to_fasta.1.py

# Shu-Ting Cho <vivianlily6@hotmail.com>
# convert fastq to fasta and qual
# requires:  biopython

# v1 2020/03/20

# Usage: python3 /home/shutingcho/pyscript/convert_fastq_to_fasta.1.py --in_fastq=/scratch/shutingcho/phyto57/phyto57.08/phrap/regionX.fastq --out_fasta=/scratch/shutingcho/phyto57/phyto57.08/phrap/regionX.fasta --out_qual=/scratch/shutingcho/phyto57/phyto57.08/phrap/regionX.qual

# for reading fasta file
from Bio import SeqIO

# for making directory
import os

# get options
import sys, getopt
opts, args = getopt.getopt(sys.argv[1:], '', longopts=[
	'in_fastq=',
	'out_fasta=',
	'out_qual='])

# get variables from opts
for opt, arg in opts:
	if opt == "--in_fastq":
		try:
			infastq = list(SeqIO.parse(str(arg), "fastq")) # check file exist and format
			infastq = str(arg)
		except IOError:
			print("Could not read file:", str(arg))
	elif opt == "--out_fasta":
		out_fasta = str(arg)
		out_dir = os.path.dirname(out_fasta)
		if not os.path.exists(out_dir): # make dir if not exist
			os.makedirs(out_dir)
	elif opt == "--out_qual":
		out_qual = str(arg)
		out_dir = os.path.dirname(out_qual)
		if not os.path.exists(out_dir): # make dir if not exist
			os.makedirs(out_dir)
	else:
		assert False, "unhandled option"

# convert(in_file, in_format, out_file, out_format, alphabet=None)
SeqIO.convert(infastq, "fastq", out_fasta, "fasta")
SeqIO.convert(infastq, "fastq", out_qual, "qual")

