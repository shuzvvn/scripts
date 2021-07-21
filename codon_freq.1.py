#!/usr/bin/python3

# codon_freq.1.py

# Shu-Ting Cho <vivianlily6@hotmail.com>
# calculate the codon frequencies of sequences under a given directory for genomegaMap analyses
# requires:  biopython
# v1 2021/07/21

# Usage:
# python3 /home/shu-ting/pyscript/codon_freq.1.py --in_dir=/home/shu-ting/project/ecoli01/ecoli01.35/pal2nal/


## import modules
import argparse
import os

from Bio import SeqIO


# parsing arguments
def parse_args():
	parser = argparse.ArgumentParser(description='Calculate the codon frequencies of sequences under a given directory for genomegaMap analyses.')
	parser.add_argument('-i', '--in_dir', help='input directory')
	parser.add_argument('-e', '--in_ext', default="fasta", help='read files with extension (default: fasta)')
	return parser.parse_args()

# calculate codon frequency
def codon_freq(seq):
	Seq_Cai = CodonUsage.CodonAdaptationIndex()
	return Seq_Cai.cai_for_gene(seq)


# main
def main():
	args = parse_args()

	# codon order list
	codons = ['TTT', 'TTC', 'TTA', 'TTG', 'TCT', 'TCC', 'TCA', 'TCG', 'TAT', 'TAC', 'TGT', 'TGC', 'TGG', 'CTT', 'CTC', 'CTA', 'CTG', 'CCT', 'CCC', 'CCA', 'CCG', 'CAT', 'CAC', 'CAA', 'CAG', 'CGT', 'CGC', 'CGA', 'CGG', 'ATT', 'ATC', 'ATA', 'ATG', 'ACT', 'ACC', 'ACA', 'ACG', 'AAT', 'AAC', 'AAA', 'AAG', 'AGT', 'AGC', 'AGA', 'AGG', 'GTT', 'GTC', 'GTA', 'GTG', 'GCT', 'GCC', 'GCA', 'GCG', 'GAT', 'GAC', 'GAA', 'GAG', 'GGT', 'GGC', 'GGA', 'GGG']
	codon_freqs = []

	# store all seqs as one string
	seq = ''

	# read files with ext under dir
	in_dir = args.in_dir.rstrip('/')
	for file in os.listdir(in_dir):
		if file.endswith(args.in_ext):
			in_file = in_dir + '/' + file

			# read fasta
			for record in SeqIO.parse(in_file, "fasta"):
				seq = seq + str(record.seq)

	# Split seq into list of codons
	seq_codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
	total = len(seq_codons)

	# get count and relative frequencies (count/total) for each codon in codons list
	for codon in codons:
		count = seq_codons.count(codon)
		freq = count/total
		codon_freqs.append(freq)

	# print result to screen
	print(*codon_freqs)

if __name__ == '__main__':
	main()