#!/usr/bin/python3

# find_DNA_motif.1.py

# Shu-Ting Cho <vivianlily6@hotmail.com>
# find specified seq motif in fasta, report seq and position
# this version will always treat seq as circular

# v1 2019/05/13


# Usage: python3 /home/shutingcho/py/find_DNA_motif.1.py --nt_fasta=/home/shutingcho/Downloads/db.fasta --motif="ATTTCGGTTGTAGC" --min_match=12

import sys, getopt

# for reading fasta file
from Bio import SeqIO
# for reverse complement
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

# get options
opts, args = getopt.getopt(sys.argv[1:], '', longopts=[
	'nt_fasta=',
	'motif=',
	'min_match='])

# get variables from opts
for opt, arg in opts:
	if opt == "--nt_fasta":
		record_dict = SeqIO.to_dict(SeqIO.parse(str(arg), "fasta"))
	elif opt == "--motif":
		motif_seq = str(arg)
	elif opt == "--min_match":
		min_match_num = int(arg)
	else:
		assert False, "unhandled option"

# function to compare seq and count score
def score_seqs(seq1, seq2):
	site = 0
	score = 0
	while site < len(seq2):
		if seq1[site] == seq2[site]:
			score = score + 1
		site = site + 1
	return score

# reverse complement seq
def reverse_complement(seq):
	my_dna = Seq(seq, generic_dna)
	seq_RC = str(my_dna.reverse_complement())
	return seq_RC

# check inputs and report
# min_match_num cant be longer than motif length
if min_match_num < len(motif_seq):
	sum_dict = {}
	print('input_seq = %i\nmotif_lgth = %i' % (len(record_dict), len(motif_seq)))
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
			if window_start < 0:
				sliding_word = str(record_dict[seq_id].seq)[window_start:] + str(record_dict[seq_id].seq)[:window_end]
			else:
				sliding_word = str(record_dict[seq_id].seq)[window_start:window_end]
			# compare seqs, check both strands
			if score_seqs(sliding_word, motif_seq) >= min_match_num:
				print(seq_id + '\t' + str(window_start + 1) + '\t' + str(window_end) + '\t' + '+' + '\t' + sliding_word + '\t' + str(score_seqs(sliding_word, motif_seq)))
				hit_num = hit_num + 1
			elif score_seqs(reverse_complement(sliding_word), motif_seq) >= min_match_num:
				print(seq_id + '\t' + str(window_start + 1) + '\t' + str(window_end) + '\t' + '-' + '\t' + reverse_complement(sliding_word) + '\t' + str(score_seqs(reverse_complement(sliding_word), motif_seq)))
				hit_num = hit_num + 1
			window_end = window_end + 1
		sum_dict[seq_id] = hit_num
		print('\n')
	# print report
	for seq_id in sum_dict:
		print(seq_id + ' = ' + str(sum_dict[seq_id]))
else:
	print('!!ERROR: min_match cannot be longer than motif length!')
