#!/usr/bin/python3

# get_prior_genomegamap.2.py

# Shu-Ting Cho <vivianlily6@hotmail.com>
# calculate the codon frequencies of sequences under a given directory for genomegaMap analyses
# requires:  biopython
# v1 2021/07/21
# v2 2021/08/19 handle alignment with gaps, by excluding all "-", ".", and "N" from the concat seq, for calculating codon usage

# Usage:
# python3 /home/shu-ting/pyscript/get_prior_genomegamap.2.py --in_dir=/home/shu-ting/project/ecoli01/ecoli01.35/pal2nal.rename/ --out_file=/home/shu-ting/project/ecoli01/ecoli01.35/genomegamap/random_align_334.fasta


## import modules
import argparse
import os

from Bio import SeqIO
import random

# parsing arguments
def parse_args():
	parser = argparse.ArgumentParser(description='Calculate the codon frequencies of sequences under a given directory for genomegaMap analyses.')
	parser.add_argument('-i', '--in_dir', help='input directory')
	parser.add_argument('-e', '--in_ext', default="fasta", help='read files with extension (default: fasta)')
	parser.add_argument('-o', '--out_file', default="./out.fasta", help='write random aligned codons to file (default: ./out.fasta)')
	parser.add_argument('-n', '--n_codons', type=int, default=334, help='number of random codons to choose from inputs (default: 334)')
	return parser.parse_args()


# main
def main():
	args = parse_args()

	# make dir if not exist
	out_file = args.out_file
	out_dir = os.path.dirname(out_file)
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)

	# store all seqs as one string
	seq = ''
	seq_dict = {}

	# read files with ext under dir
	print("Reading input files ...")
	in_dir = args.in_dir.rstrip('/')
	count_in_files = 0
	for file in os.listdir(in_dir):
		count_in_files += 1
		if file.endswith(args.in_ext):
			in_file = in_dir + '/' + file

			# read fasta
			for record in SeqIO.parse(in_file, "fasta"):
				ID_h = str(record.id)
				seq_h = str(record.seq)
				
				# concat all seqs for calculation codon usage, exclude gaps
				seq = seq + seq_h.replace('-', '').replace('.', '').replace('N', '')

				# concat seqs for each gene
				try:
					seq_dict[ID_h] = seq_dict[ID_h] + seq_h
				except:
					seq_dict[ID_h] = seq_h


	# length of concat alignment
	len_align = len(seq_dict[list(seq_dict.keys())[0]])
	print("count_genomes = " + str(len(seq_dict)) + ", count_genes = " + str(count_in_files) + ", align_length = " + str(len_align))

	n_codon_align = int(len_align/3)
	
	# get random number for n_codon_align
	n_codons = int(args.n_codons)
	print("Choose " + str(n_codons) + " random codons from " + str(n_codon_align) + " aligned codons ...")
	random_codons_pos = random.choices(range(n_codon_align), k=n_codons)

	# generate alignment base on position chosen
	# write output
	out_file_h = open(out_file, 'w')
	for ID_h in seq_dict:
		out_file_h.write(">" + ID_h + "\n")
		random_seq_h = ''
		for pos in random_codons_pos:
			random_seq_h = random_seq_h + seq_dict[ID_h][pos*3: pos*3+3]
		out_file_h.write(random_seq_h + "\n")
	out_file_h.close()


	# Split concated all seq into list of codons
	print("Calculate codon frequencies ...")
	seq_codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
	total = len(seq_codons)

	# codon order list
	codons = ['TTT', 'TTC', 'TTA', 'TTG', 'TCT', 'TCC', 'TCA', 'TCG', 'TAT', 'TAC', 'TGT', 'TGC', 'TGG', 'CTT', 'CTC', 'CTA', 'CTG', 'CCT', 'CCC', 'CCA', 'CCG', 'CAT', 'CAC', 'CAA', 'CAG', 'CGT', 'CGC', 'CGA', 'CGG', 'ATT', 'ATC', 'ATA', 'ATG', 'ACT', 'ACC', 'ACA', 'ACG', 'AAT', 'AAC', 'AAA', 'AAG', 'AGT', 'AGC', 'AGA', 'AGG', 'GTT', 'GTC', 'GTA', 'GTG', 'GCT', 'GCC', 'GCA', 'GCG', 'GAT', 'GAC', 'GAA', 'GAG', 'GGT', 'GGC', 'GGA', 'GGG']
	codon_freqs = []

	# calculate codon frequency
	# get count and relative frequencies (count/total) for each codon in codons list
	for codon in codons:
		count = seq_codons.count(codon)
		freq = count/total
		codon_freqs.append(freq)

	# print result to screen
	print(*codon_freqs)

if __name__ == '__main__':
	main()