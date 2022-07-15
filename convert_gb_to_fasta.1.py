#!/usr/bin/env python3

# convert_gb_to_fasta.1.py

# Shu-Ting Cho <vivianlily6@hotmail.com>
# parse whole genome fasta from gb
# requires:  biopython
# v1 2022/07/14

# Usage:
# ~/script/convert_gb_to_fasta.1.py --in_dir=/Users/stc/project/DVT_ecoli01/gb/ --out_dir=/Users/stc/project/DVT_ecoli01/fasta/

## import modules
from Bio import SeqIO
import argparse
import os

# parsing arguments
def parse_args():
	parser = argparse.ArgumentParser(description='Parse whole genome fasta from gbk.')
	parser.add_argument('-i', '--in_dir', help='input dir')
	parser.add_argument('-e', '--in_ext', default="gbk", help='read files with extension (default: gbk)')
	parser.add_argument('-o', '--out_dir', default="./", help='output dir (default: current working directory)')
	parser.add_argument('-v', '--verbose', default=True, help='report counts')
	return parser.parse_args()

# make dir if not exist
def make_outdir(out_dir):
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)

# main
def main():
	args = parse_args()

	# write output
	out_dir = args.out_dir.rstrip('/')
	make_outdir(out_dir)

	# read files with ext under dir
	in_dir = args.in_dir.rstrip('/')
	in_ext = args.in_ext
	file_count = 0
	for file in os.listdir(in_dir):
		if file.endswith(in_ext):
			in_file = in_dir + '/' + file
			genomeID = file[:-(len(in_ext)+1)]
			out_file = out_dir + '/' + genomeID + '.fasta'

			# convert gbk to fasta
			SeqIO.convert(in_file, "genbank", out_file, "fasta")
			file_count+=1
	
	verbose = args.verbose
	if verbose:
		print('file_count = %i' % (file_count))

if __name__ == '__main__':
	main()
	