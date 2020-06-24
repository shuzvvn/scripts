#!/usr/bin/python3

# prot_mw.1.py

# Shu-Ting Cho <vivianlily6@hotmail.com>
# Calculates the molecular weight of a protein.
# v1 2020/06/24

# Usage: python3 /home/shutingcho/py/prot_mw.1.py --in_file=/scratch/shutingcho/phyto21/phyto21.56/run6/signalp/phyto21.aa.fasta.1 --out_file=/scratch/shutingcho/phyto21/phyto21.56/run6/signalp/phyto21.aa.mw

import sys, getopt

# for reading fasta file
from Bio import SeqIO

# Analyzing protein sequences with the ProtParam module.
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# get options
opts, args = getopt.getopt(sys.argv[1:], '', longopts=[
	'in_file=',
	'out_file='])

# get variables from opts
for opt, arg in opts:
	if opt == "--in_file":
		aafasta_records = list(SeqIO.parse(str(arg), 'fasta'))
	elif opt == "--out_file":
		out_file = str(arg)
	else:
		assert False, "unhandled option"


# Calculates the molecular weight of all protein seq in the aa fasta file
count_seq = 0

out_file_h = open(out_file, 'w')

for i in aafasta_records:
	locus_tag = str(i.id)
	my_seq = str(i.seq)
	analysed_seq = ProteinAnalysis(my_seq)
	kDa = analysed_seq.molecular_weight() / 1000
	out_file_h.write(locus_tag + '\t' + str(round(kDa, 2)) + '\n')
	count_seq += 1

out_file_h.close()

# print report
print('count_seq =', count_seq)
# end of script