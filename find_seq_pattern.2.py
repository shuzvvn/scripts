#!/usr/bin/python3

# find_seq_pattern.2.py

# Shu-Ting Cho <vivianlily6@hotmail.com>
# find specified seq pattern in fasta, report seq and position
# Seq or a MutableSeq of DNA are acceptable
# For tDNA border finding or regulator binding motifs

# v1 
# v2 2019/05/15



# Usage: python3 /home/shutingcho/py/find_seq_pattern.2.py --nt_fasta=/home/shutingcho/Downloads/db.fasta --pattern="ATTTCGGTTGTAGC"

import sys, getopt, re

# for reading fasta file
from Bio import SeqIO
# for reverse complement
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

# get options
opts, args = getopt.getopt(sys.argv[1:], '', longopts=[
	'nt_fasta=',
	'pattern='])

# get variables from opts
for opt, arg in opts:
	if opt == "--nt_fasta":
		record_dict = SeqIO.to_dict(SeqIO.parse(str(arg), "fasta"))
	elif opt == "--pattern":
		pattern_seq = str(arg)
	else:
		assert False, "unhandled option"

# reverse complement seq
def reverse_complement(seq):
	my_dna = Seq(seq, IUPAC.ambiguous_dna)
	seq_RC = str(my_dna.reverse_complement())
	return seq_RC

# mutable DNA seq to re and capital
def mut2RE(seq):
	seq = seq.upper()
	for ori, REreplacement in IUPAC_dict.items():
		seq = seq.replace(ori, REreplacement)
	return seq

# IUPAC dict
IUPAC_dict = {'A':'A', 'T':'T', 'C':'C', 'G':'G', 'M':'[AC]', 'R':'[AG]', 'W':'[AT]', 'S':'[CG]', 'Y':'[CT]', 'K':'[GT]', 'V':'[ACG]', 'H':'[ACT]', 'D':'[AGT]', 'B':'[CGT]', 'X':'[ACGT]', 'N':'[ACGT]'}

# check inputs and report
# check if pattern_seq valid
for character in pattern_seq:
	if character not in IUPAC_dict:
		print('!!ERROR: pattern contains NOT DNA characters: ' + character)
		exit()

print('input_seq = %i\npattern_lgth = %i' % (len(record_dict), len(pattern_seq)))

# init result summary
sum_dict = {}

# REpattern dict 
pattern_dict = {'1':mut2RE(pattern_seq), '-1':mut2RE(reverse_complement(pattern_seq))}

# start looking for pattern
for seq_id in record_dict:
	# init hit
	hit_num = 0
	for strand in pattern_dict:
		regex = re.compile(pattern_dict[strand])
		for match in re.finditer(regex, str(record_dict[seq_id].seq)):
			print('%s\t%i\t%i\t%s\t%s' % (seq_id , match.start() , match.end(), strand, match.group()))
			hit_num = hit_num + 1

	sum_dict[seq_id] = hit_num
	
# print report
for seq_id in sum_dict:
	print(seq_id + ' = ' + str(sum_dict[seq_id]))