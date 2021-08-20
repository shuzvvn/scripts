#!/usr/bin/python3

# pairwise_ident_average.2.py

# Shu-Ting Cho <vivianlily6@hotmail.com>
# calculate average of pairwise identities for each fasta in the input dir and output table
# requires:  biopython
# v1 2021/07/21
# v2 2021/08/19 mask gaps

# Usage:
# python3 /home/shu-ting/pyscript/pairwise_ident_average.2.py --in_dir=/home/shu-ting/project/ecoli01/ecoli01.35/aa_alignment/ --out_file=/home/shu-ting/project/ecoli01/ecoli01.35/aa_alignment.ident


## import modules
import argparse
import os

from Bio import SeqIO
from itertools import combinations
from statistics import mean

# parsing arguments
def parse_args():
	parser = argparse.ArgumentParser(description='Calculate average of pairwise identities for each fasta in the input dir and output table.')
	parser.add_argument('-i', '--in_dir', help='input directory')
	parser.add_argument('-e', '--in_ext', default="fasta", help='read files with extension (default: fasta)')
	parser.add_argument('-o', '--out_file', default="./", help='output table (default: current working directory)')
	return parser.parse_args()

# make dir if not exist
def make_outdir(out_file):
	out_dir = os.path.dirname(out_file)
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)

# get all combinations of pairs from list of seq IDs
def pairs(IDs):
	pairs_list = combinations(IDs, 2)
	return pairs_list

# get identity of two seqs
def identity(seq_pair):
	seq1 = seq_pair[0]
	seq2 = seq_pair[1]
	if seq1 == seq2:
		identity = 1
	else:
		bases = len(seq1)
		matches = 0
		for pos in range(bases):
			if seq1[pos] == seq2[pos]:
				matches += 1
		identity = matches / bases
	return identity

# main
def main():
	args = parse_args()

	# write output
	out_file = args.out_file
	make_outdir(out_file)
	out_file_h = open(out_file, 'w')
	out_file_h.write("geneID,ave_ident\n") # title

	# read files with ext under dir
	in_dir = args.in_dir.rstrip('/')
	in_ext = args.in_ext
	for file in os.listdir(in_dir):
		if file.endswith(in_ext):
			in_file = in_dir + '/' + file
			geneID = file[:-(len(in_ext)+1)]

			# read fasta
			record_dict = SeqIO.index(in_file, "fasta")
			IDs = [*record_dict] # list of IDs
			ID_pair_list = pairs(IDs) # list of paired IDs

			# get the list of all gap positions
			gap_pos_l = []
			for ID in IDs:
				gap_pos_l += [pos for pos, char in enumerate(str(record_dict[ID].seq)) if char == '-']
			gap_pos_l = list(set(gap_pos_l))

			if len(gap_pos_l) != 0:
				# remove gap pos
				seq_len = len(str(record_dict[IDs[0]].seq))
				nogap_pos_l = [pos for pos in list(range(seq_len)) if pos not in gap_pos_l]
				# get identities of masked pairwise seqs
				identity_list = []
				for ID_pair in ID_pair_list: # get identity for each pair of seqs
					seq1 = record_dict[ID_pair[0]].seq
					seq2 = record_dict[ID_pair[1]].seq
					nseq1 = ''
					nseq2 = ''
					for pos in range(seq_len):
						if pos in nogap_pos_l:
							nseq1 += str(record_dict[ID_pair[0]].seq)[pos]
							nseq2 += str(record_dict[ID_pair[1]].seq)[pos]
					seq_pair = [nseq1, nseq2]
					identity_list.append(identity(seq_pair))
			else: # if no gaps
				# get identities of original pairwise seqs
				identity_list = []
				for ID_pair in ID_pair_list: # get identity for each pair of seqs
					seq_pair = [record_dict[ID_pair[0]].seq, record_dict[ID_pair[1]].seq]
					identity_list.append(identity(seq_pair))

			# average of identity_list
			ave_ident = mean(identity_list)
			out_file_h.write(geneID + "," + str(ave_ident) + "\n")

			# close record_dict
			record_dict.close()

	# close output
	out_file_h.close()

if __name__ == '__main__':
	main()