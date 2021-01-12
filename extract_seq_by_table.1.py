#!/usr/bin/python3

# extract_seq_by_table.1.py

# Shu-Ting Cho <vivianlily6@hotmail.com>
# extract seq (start, end), rename, reorder base on input table
# in_file fasta
# table tsv
# out_file fasta
# v1 2021/01/10

# Usage: 
# python3 /home/shu-ting/pyscript/extract_seq_by_table.1.py --in_file=/home/shu-ting/scratch/ecoli01/ecoli01.13/source_data/all.fasta --table=/home/shu-ting/scratch/ecoli01/ecoli01.13/blastn/blaCTX-M-15.all --out_file=/home/shu-ting/scratch/ecoli01/ecoli01.13/blastn/hits.1.fasta --index_in=2 --index_out=2 --index_start=7 --index_end=8 --index_strand=9 --upstream=5000 --downstream=5000 --verbose=1


## import modules
import sys, getopt, os, csv

# for reading fasta file
from Bio import SeqIO

# for reverse complement
from Bio.Seq import Seq

# for reading tsv file as DataFrame
# import pandas as pd



## set default values
index_in = 0
index_out = None
verbose = True

index_start = None
index_end = None

seq_pos_start = 0
seq_pos_end = -1

upstream = 0
downstream = 0

index_strand = None


## read arguments
# get options
opts, args = getopt.getopt(sys.argv[1:], '', longopts=[
	'in_file=',
	'table=',
	'index_in=',
	'index_out=',
	'index_start=',
	'index_end=',
	'index_strand=',
	'upstream=',
	'downstream=',
	'verbose=',
	'out_file='])

# get variables from opts
for opt, arg in opts:
	if opt == "--in_file":
		try:
			record_dict = SeqIO.to_dict(SeqIO.parse(str(arg), "fasta"))
		except IOError:
			outmessage = "Could not read file:" + str(arg) + ", not in fasta format or does not exist."
			sys.exit(outmessage)
	elif opt == "--table":
		try:
			file_h = open(str(arg), 'r')
		except IOError:
			outmessage = "File does not exist: " + str(arg)
			sys.exit(outmessage)
		try:
			reader = csv.reader(file_h, delimiter = '\t')
			ncol = len(next(reader)) # Read first line and count columns
		except ValueError:
			outmessage = "File not in tsv format: " + str(arg)
			sys.exit(outmessage)
	elif opt == "--index_in":
		try:
			index_in = int(arg)
		except ValueError:
			print("Not valid number for index_in. Proceed with index_in = 0 ...")
		if index_in + 1 > ncol:
			print("index_in" + str(index_in) + "exceed ncol" + str(ncol) + ". Proceed with index_in = 0 ...")
	elif opt == "--index_out":
		try:
			index_out = int(arg)
		except ValueError:
			print("Not valid number for index_out. Proceed with index_out = None ...")
		if index_out + 1 > ncol:
			print("index_out" + str(index_out) + "exceed ncol" + str(ncol) + ". Proceed with index_out = index_in ...")
	elif opt == "--index_start":
		try:
			index_start = int(arg)
		except ValueError:
			print("Not valid number for index_start. Proceed with index_start = 0 ...")
		if index_start + 1 > ncol:
			print("index_start" + str(index_start) + "exceed ncol" + str(ncol) + ". Proceed with index_start = None ...")
	elif opt == "--index_end":
		try:
			index_end = int(arg)
		except ValueError:
			print("Not valid number for index_end. Proceed with index_end = None ...")
		if index_end + 1 > ncol:
			print("index_end" + str(index_end) + "exceed ncol" + str(ncol) + ". Proceed with index_end = None ...")
	elif opt == "--index_strand":
		try:
			index_strand = int(arg)
		except ValueError:
			print("Not valid number for index_strand. Proceed with index_strand = None ...")
	elif opt == "--upstream":
		try:
			upstream = int(arg)
		except ValueError:
			print("Not valid number for upstream. Proceed with upstream = 0 ...")
	elif opt == "--downstream":
		try:
			downstream = int(arg)
		except ValueError:
			print("Not valid number for downstream. Proceed with downstream = 0 ...")
	elif opt == "--out_file":
		out_file = str(arg)
	elif opt == "--verbose":
		try:
			verbose = bool(arg)
		except ValueError:
			print("value for verbose should be 0 or 1. Proceed with verbose = 1 ...")
	else:
		assert False, "unhandled option"


# assign value if not given
if index_out is None:
	index_out = index_in




# make dir if not exist
import os
out_dir = os.path.dirname(out_file)
if not os.path.exists(out_dir):
	os.makedirs(out_dir)



## main ##

# count seqs in input file
count_in = len(record_dict)
# count items in input list
count_table = 0
# count seqs in output file
count_out = 0

# write output
out_file_h = open(out_file, 'w')

for line in file_h:
	count_table += 1
	words = line.strip('\n').split('\t')

	# get seq in input file
	try:
		# get names
		seq_name_in = words[index_in]
		seq_name_out = words[index_out]

		record_h = record_dict[seq_name_in]

		# get strand
		if index_strand != None:
			seq_strand = words[index_strand]
			if seq_strand == "plus" or seq_strand == "+" or seq_strand == "1":
				seq_strand = 1
			elif seq_strand == "minus" or seq_strand == "-" or seq_strand == "-1":
				seq_strand = -1
			else:
				print("strand for " + seq_name_out + "not recognized. Proceed with strand = 1 ..." )
				seq_strand = 1
		else:
			seq_strand = 1

		# get positions
		if index_start != None:
			seq_pos_start = int(words[index_start])
			if seq_strand == 1:
				if seq_pos_start - upstream < 0:
					seq_out_start = 1
				else:
					seq_out_start = seq_pos_start - upstream
			else:
				if seq_pos_start - downstream < 0:
					seq_out_start = 1
				else:
					seq_out_start = seq_pos_start - downstream

		if index_end != None:
			seq_pos_end = int(words[index_end])
			if seq_strand == 1:
				if seq_pos_end + downstream > len(record_h.seq):
					seq_out_end = len(record_h.seq)
				else:
					seq_out_end = seq_pos_end + downstream
			else:
				if seq_pos_end + upstream > len(record_h.seq):
					seq_out_end = len(record_h.seq)
				else:
					seq_out_end = seq_pos_end + upstream

		# get seq
		if seq_strand == 1:
			seq_h = record_h.seq[seq_out_start-1:seq_out_end]
		else:
			seq_h = record_h.seq[seq_out_start-1:seq_out_end].reverse_complement()

		# write output		
		if index_start == None and index_end == None and index_strand == None:
			out_file_h.write('>' + str(seq_name_out) + '\n' + str(seq_h) + '\n')
		else:
			out_file_h.write('>' + str(seq_name_out) + "_" + str(seq_out_start) + "_" + str(seq_out_end) + "_" + str(seq_strand) + '\n' + str(seq_h) + '\n')
		# count output
		count_out += 1

	except:
		print(str(seq_name_in) + ' not found!!!')


# close files
out_file_h.close()
file_h.close()


# print verbose
if verbose:
	print('count_in = %i, count_table = %i, count_out = %i' % (count_in, count_table, count_out))
