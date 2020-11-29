#!/usr/bin/python3

# extract_fasta_by_id.1.py

# Shu-Ting Cho <vivianlily6@hotmail.com>
# extract seq, rename, reorder base on input list
# in_file fasta
# list_file tsv
# out_file fasta
# v1 2020/11/29

# Usage: 
# python3 /home/shu-ting/pyscript/extract_fasta_by_id.1.py --in_file=/home/shu-ting/scratch/ecoli01/ecoli01.04/source_data/nt_fasta/ecoli01.04.fasta --list_file=/home/shu-ting/project/ecoli01/ecoli01.04/gyrA.tsv --out_file=/home/shu-ting/project/ecoli01/ecoli01.04/unaligned/gyrA.nt.fasta --index_in=0 --index_out=0 --verbose=1


## import modules
import sys, getopt, os

# for reading fasta file
from Bio import SeqIO

# for reading tsv file as DataFrame
# import pandas as pd



## set default values
index_in = 0
index_out = None
verbose = True



## read arguments
# get options
opts, args = getopt.getopt(sys.argv[1:], '', longopts=[
	'in_file=',
	'list_file=',
	'index_in=',
	'index_out=',
	'verbose=',
	'out_file='])

# get variables from opts
for opt, arg in opts:
	if opt == "--in_file":
		try:
			record_dict = SeqIO.to_dict(SeqIO.parse(str(arg), "fasta"))
		except IOError:
			print("Could not read file:", str(arg), ", not in fasta format or does not exist.")
	elif opt == "--list_file":
		try:
			list_file_h = open(str(arg), 'r')
		except IOError:
			print('list_file not found!')
	elif opt == "--index_in":
		try:
			index_in = int(arg)
		except ValueError:
			print("Oops!  That was no valid number for index_in. Try again...")
	elif opt == "--index_out":
		try:
			index_out = int(arg)
		except ValueError:
			print("Oops!  That was no valid number for index_out. Try again...")
	elif opt == "--out_file":
		out_file = str(arg)
	elif opt == "--verbose":
		verbose = bool(arg)
	else:
		assert False, "unhandled option"


# assign value if not given
if index_out is None:
	index_out = index_in


# tsv to df
# try:
# 	# Read a tsv file into DataFrame.
# 	tsv_read = pd.read_csv(list_file, sep='\t', header=None, usecols=[index_in,index_out])
# except IOError:
# 	print("Something wrong with list_file, index_in, or index_out.")


# df to dict
# my_dict = tsv_read.to_dict()



# make dir if not exist
import os
out_dir = os.path.dirname(out_file)
if not os.path.exists(out_dir):
	os.makedirs(out_dir)



## main ##

# count seqs in input file
count_in = len(record_dict)
# count items in input list
count_list = 0
# count seqs in output file
count_out = 0

# write output
out_file_h = open(out_file, 'w')

for line in list_file_h:
	count_list += 1
	words = line.strip('\n').split('\t')

	# get names
	seq_name_in = words[index_in]
	seq_name_out = words[index_out]

	# get seq in input file
	try:
		record_h = record_dict[seq_name_in]
		seq_h = record_h.seq
		out_file_h.write('>' + str(seq_name_out) + '\n' + str(seq_h) + '\n')
		# count output
		count_out += 1
	except:
		print(str(seq_name_in) + ' not found!!!')

out_file_h.close()
list_file_h.close()

# print verbose
if verbose:
	print('count_in = %i, count_list = %i, count_out = %i' % (count_in, count_list, count_out))