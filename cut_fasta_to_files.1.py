#!/usr/bin/python3

# cut_fasta_to_files.1.py

# Shu-Ting Cho <vivianlily6@hotmail.com>
# cut seq in input fasta to individual fasta files
# v1 2021/01/11

# Usage: 
# python3 /home/shu-ting/pyscript/cut_fasta_to_files.1.py --in_file=/home/shu-ting/scratch/ecoli01/ecoli01.13/blastn/hits.1.fasta --out_dir=/home/shu-ting/scratch/ecoli01/ecoli01.13/blastn/hits.2.fasta --verbose=1


## import modules
import sys, getopt, os

# for reading fasta file
from Bio import SeqIO


## set default values
verbose = True



## read arguments
# get options
opts, args = getopt.getopt(sys.argv[1:], '', longopts=[
	'in_file=',
	'verbose=',
	'out_dir='])

# get variables from opts
for opt, arg in opts:
	if opt == "--in_file":
		try:
			records = list(SeqIO.parse(str(arg), "fasta"))
		except IOError:
			outmessage = "Could not read file:" + str(arg) + ", not in fasta format or does not exist."
			sys.exit(outmessage)
	elif opt == "--out_dir":
		out_dir = str(arg).rstrip('/')
	elif opt == "--verbose":
		try:
			verbose = bool(arg)
		except ValueError:
			print("value for verbose should be 0 or 1. Proceed with verbose = 1 ...")
	else:
		assert False, "unhandled option"


# make dir if not exist
if not os.path.exists(out_dir):
	os.makedirs(out_dir)


## main ##

# count seqs in input file
count_in = len(records)


for record in records:
	# write output
	out_file = out_dir + '/' + str(record.id) + ".fasta"
	out_file_h = open(out_file, 'w')
	out_file_h.write('>' + str(record.id) + '\n' + str(record.seq))
	out_file_h.close()


# print verbose
if verbose:
	print('count_in = %i' % (count_in))