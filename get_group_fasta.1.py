#!/usr/bin/python3

# get_group_fasta.1.py

# Shu-Ting Cho <vivianlily6@hotmail.com>
# extract seq and/or rename base on input group
# in_fasta
# in_group tsv or csv
# out_dir
# v1 2021/03/05

# Usage: 
# python3 /home/shu-ting/pyscript/get_group_fasta.1.py --in_fasta=/scratch/shutingcho/phyto56/phyto56.11/source_data/aa_fasta/phyto56.11.fasta --in_group=/scratch/shutingcho/phyto56/phyto56.11/orthomcl/group/all+1.group --out_dir=/scratch/shutingcho/phyto56/phyto56.11/aa_fasta/ --sep="," --rename=1 --verbose=1 


## import modules
import sys, getopt, os

# for reading fasta file
from Bio import SeqIO


## set default values
verbose = 1
rename = 1
sep = ","


## read arguments
# get options
opts, args = getopt.getopt(sys.argv[1:], '', longopts=[
	'in_fasta=',
	'in_group=',
	'sep=',
	'rename=',
	'out_dir=',
	'verbose='])

# get variables from opts
for opt, arg in opts:
	if opt == "--in_fasta":
		in_fasta = str(arg)
		try:
			record_dict = SeqIO.to_dict(SeqIO.parse(in_fasta, "fasta"))
		except IOError:
			print("Could not read file:", in_fasta, ", not in fasta format or does not exist.")
	elif opt == "--in_group":
		try:
			in_group_h = open(str(arg), 'r')
		except IOError:
			print('in_group not found!')
	elif opt == "--sep":
		sep = str(arg)
	elif opt == "--out_dir":
		out_dir = str(arg).rstrip("/")
	elif opt == "--rename":
		rename = bool(arg)
	elif opt == "--verbose":
		verbose = bool(arg)
	else:
		assert False, "unhandled option"



# make dir if not exist
import os
if not os.path.exists(out_dir):
	os.makedirs(out_dir)



## main ##

# count seqs in input file
count_in = len(record_dict)

# count seqs in output file
count_out = 0

# read groups to list
lines = in_group_h.readlines()

# count items in input list
count_group = 0

# first line is title (Genome ID)
genomeIDs = lines[0].strip('\n').split(sep)


# process from the second line
for line in lines[1:]:
	count_group += 1
	words = line.strip('\n').split(sep)
	groupID = str(words[0])
	seqIDs = words[1:]
	# write output
	out_file = out_dir + "/" + groupID + ".fasta"
	out_file_h = open(out_file, 'w')

	index_h = 0
	for seqID in seqIDs:
		index_h += 1
		seq_name_in = str(seqID)
		if rename:
			seq_name_out = str(genomeIDs[index_h])
		else:
			seq_name_out = seq_name_in
	
		# get seq in input file
		try:
			record_h = record_dict[seq_name_in]
			seq_h = record_h.seq
			out_file_h.write('>' + str(seq_name_out) + '\n' + str(seq_h) + '\n')
			# count output
			count_out += 1
		except:
			print('Cannot find ' + str(seq_name_in) + ' in group ' + groupID)

	out_file_h.close()

in_group_h.close()

# print verbose
if verbose:
	print('in_fasta = %s, %i seqs' % (in_fasta, count_in))
	print('Processed %i groups and %i seqs' % (count_group, count_out))