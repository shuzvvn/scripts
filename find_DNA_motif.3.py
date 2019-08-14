#!/usr/bin/python3

# find_DNA_motif.3.py

# Shu-Ting Cho <vivianlily6@hotmail.com>
# find specified seq motif in fasta, report seq and position
# this version will always treat seq as circular

# v3 2019/08/14
#    use 'max_mismatch'
#    modify output format 'rep_sum'
# v2 2019/05/15
#    allow IUPAC ambiguous DNA as motif
#    add check reverse complement option
#    modify output format
# v1 2019/05/13

# Usage: python3 /home/shutingcho/py/find_DNA_motif.3.py --nt_fasta=/home/shutingcho/Downloads/db.fasta --motif="RTTDCAWWTGHAAY" --max_mismatch=0 --check_RC=0 --rep_sum=0

import sys, getopt

# for reading fasta file
from Bio import SeqIO
# for reverse complement
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


# set defalt
check_RC = 0
rep_sum = 0

# get options
opts, args = getopt.getopt(sys.argv[1:], '', longopts=[
	'nt_fasta=',
	'motif=',
	'max_mismatch=',
	'check_RC=',
	'rep_sum='])

# get variables from opts
for opt, arg in opts:
	if opt == "--nt_fasta":
		record_dict = SeqIO.to_dict(SeqIO.parse(str(arg), "fasta"))
	elif opt == "--motif":
		motif_seq = str(arg).upper().replace('I', 'N') # Inosine is basically "N"
		# set defalt for min match to match the whole
		min_match_num = len(motif_seq)
	elif opt == "--max_mismatch":
		min_match_num = len(motif_seq) - int(arg)
		if min_match_num < 1:
			assert False, "!!ERROR: max_mismatch is too long!"
	elif opt == "--check_RC":
		check_RC = bool(arg)
	elif opt == "--rep_sum":
		rep_sum = bool(arg)
	else:
		assert False, "unhandled option"

# function to compare seq and count score
def score_seqs(seq1, seq2):
	site = 0
	score = 0
	while site < len(seq2):
		if seq1[site] in IUPAC_dict[seq2[site]]:
			score = score + 1
		site = site + 1
	return score

# reverse complement seq
def reverse_complement(seq):
	my_dna = Seq(seq, IUPAC.ambiguous_dna)
	seq_RC = str(my_dna.reverse_complement())
	return seq_RC

# IUPAC dict
IUPAC_dict = {'A':['A'], 'T':['T'], 'C':['C'], 'G':['G'], 'M':['A','C'], 'R':['A','G'], 'W':['A','T'], 'S':['C','G'], 'Y':['C','T'], 'K':['G','T'], 'V':['A','C','G'], 'H':['A','C','T'], 'D':['A','G','T'], 'B':['C','G','T'], 'X':['A','C','G','T'], 'N':['A','C','G','T']}


# check inputs and report
if rep_sum:
	print('input_seq = %i\nmotif_lgth = %i' % (len(record_dict), len(motif_seq)))

sum_dict = {}
for seq_id in record_dict.keys():
	# init hit
	hit_num = 0
	# set sliding window word
	# always treat seq as circular
	# init window_end
	window_end = 1
	while window_end <= len(record_dict[seq_id].seq):
		# set window_start
		window_start = window_end - len(motif_seq)
		# get window word
		if window_start < 0:
			sliding_word = str(record_dict[seq_id].seq)[window_start:] + str(record_dict[seq_id].seq)[:window_end]
		else:
			sliding_word = str(record_dict[seq_id].seq)[window_start:window_end]
		# compare seqs, same strand
		if score_seqs(sliding_word, motif_seq) >= min_match_num:
			if window_start < 0:
				print('%s\t%i\t%i\t%i\t%i\t%s' % (seq_id, window_start, window_end, score_seqs(sliding_word, motif_seq), 1, sliding_word))
			else:
				print('%s\t%i\t%i\t%i\t%i\t%s' % (seq_id, window_start + 1, window_end, score_seqs(sliding_word, motif_seq), 1, sliding_word))
			hit_num = hit_num + 1
		elif (check_RC and score_seqs(reverse_complement(sliding_word), motif_seq) >= min_match_num):
			if window_start < 0:
				print('%s\t%i\t%i\t%i\t%i\t%s' % (seq_id, window_start, window_end, score_seqs(reverse_complement(sliding_word), motif_seq), -1, reverse_complement(sliding_word)))
			else:
				print('%s\t%i\t%i\t%i\t%i\t%s' % (seq_id, window_start + 1, window_end, score_seqs(reverse_complement(sliding_word), motif_seq), -1, reverse_complement(sliding_word)))
			hit_num = hit_num + 1
		window_end = window_end + 1
	sum_dict[seq_id] = hit_num

# print summary report
if rep_sum:
	for seq_id in sum_dict:
		print(seq_id + ' = ' + str(sum_dict[seq_id]))
