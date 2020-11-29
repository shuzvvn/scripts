#!/usr/bin/python3

# cut_seq_by_region.1.py

# Shu-Ting Cho <vivianlily6@hotmail.com>
# cut sequences in fasta by region
# in_file fasta
# start end postion
# out_file fasta
# include or exclude
# v1 2020/11/29

# Usage: 
# python3 /home/shu-ting/pyscript/cut_seq_by_region.1.py --in_file=/home/shu-ting/scratch/ecoli01/ecoli01.04/source_data/nt_fasta/ecoli01.04.fasta --out_file=/home/shu-ting/project/ecoli01/ecoli01.04/unaligned/gyrA.nt.fasta --start=80 --end=84 --include=1 --pattern=GDSAVYD --verbose=1


## import modules
import sys, getopt, re

# for reading fasta file
from Bio import SeqIO



## set default values
start = 0
end = -1
include = True
verbose = True
pattern = None
out_match = None



## read arguments
# get options
opts, args = getopt.getopt(sys.argv[1:], '', longopts=[
	'in_file=',
	'start=',
	'end=',
	'include=',
	'pattern=',
	'verbose=',
	'out_file=',
	'out_match='])

# get variables from opts
for opt, arg in opts:
	if opt == "--in_file":
		try:
			record_dict = SeqIO.to_dict(SeqIO.parse(str(arg), "fasta"))
		except IOError:
			print("Could not read file:", str(arg), ", not in fasta format or does not exist.")
	elif opt == "--start":
		try:
			start = int(arg)
		except ValueError:
			print("Oops!  That was no valid number for start. Try again...")
	elif opt == "--end":
		try:
			end = int(arg)
		except ValueError:
			print("Oops!  That was no valid number for end. Try again...")
	elif opt == "--out_file":
		out_file = str(arg)
	elif opt == "--out_match":
		out_match = str(arg)
	elif opt == "--pattern":
		pattern = re.compile(str(arg))
	elif opt == "--include":
		try:
			include = bool(arg)
		except ValueError:
			print("value for include should be a boolean")
	elif opt == "--verbose":
		try:
			verbose = bool(arg)
		except ValueError:
			print("value for verbose should be a boolean")
	else:
		assert False, "unhandled option"


# make dir if not exist
import os
out_dir = os.path.dirname(out_file)
if not os.path.exists(out_dir):
	os.makedirs(out_dir)

if out_match is not None:
	out_dir = os.path.dirname(out_match)
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	out_match_h = open(out_match, 'w')



## main ##

# count seqs in input file
count_in = len(record_dict)
# count seqs in output file
count_out = 0
# count match
count_match = 0



# write output
out_file_h = open(out_file, 'w')

# get seq in input file
for record in record_dict:
	try:
		record_h = record_dict[record]
		seq_h = record_h.seq
		cut_seq = str(seq_h[start-1:end])
		out_file_h.write('>' + str(record) + '\n' + cut_seq + '\n')
		
		# if pattern is given, check if output seq match pattern
		if pattern is not None:
			result = pattern.fullmatch(cut_seq)
			if result is None:
				out_match_h.write(str(record) + '\t' + cut_seq + '\t0\n')
			else:
				out_match_h.write(str(record) + '\t' + result.group(0) + '\t1\n')
				count_match += 1

		# count output
		count_out += 1
	except:
		print(str(record) + ' has no seq in region!!!')

out_file_h.close()

if out_match is not None:
	out_match_h.close()

# print verbose
if verbose:
	print('count_in = %i, count_out = %i' % (count_in, count_out))
if pattern is not None:
	print('count_match = %i' % (count_match))