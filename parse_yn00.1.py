#!/usr/bin/python3

# parse_yn00.1.py

# Shu-Ting Cho <vivianlily6@hotmail.com>
# parse results from yn00 (PAML)
# v1 2021/03/08

# Usage: 
# python3 /home/shu-ting/pyscript/parse_yn00.1.py --in_dir=/home/shu-ting/project/ecoli01/ecoli01.24/yn00/ --out_dir=/home/shu-ting/project/ecoli01/ecoli01.24/yn00_parsed/ --verbose=1

## import modules
import sys, getopt, os

## set default values
verbose = 1


## read arguments
# get options
opts, args = getopt.getopt(sys.argv[1:], '', longopts=[
	'in_dir=',
	'out_dir=',
	'verbose='])

# get variables from opts
for opt, arg in opts:
	if opt == "--in_dir":
		in_dir = str(arg).rstrip("/")
	elif opt == "--out_dir":
		out_dir = str(arg).rstrip("/")
	elif opt == "--verbose":
		verbose = int(arg)
	else:
		assert False, "unhandled option"

verbose = bool(verbose)

# make dir if not exist
if not os.path.exists(out_dir):
	os.makedirs(out_dir)


## main ##
count_in = 0

# List all subdirectories using scandir()
with os.scandir(in_dir) as entries:
	for entry in entries:
		if entry.is_dir():
			count_in += 1
			file_id = entry.name
			# read file to list
			dN_file = in_dir + "/" + file_id + "/2YN.dN"
			dN_file_h = open(dN_file, 'r')
			dN_lines = dN_file_h.readlines()
			dS_file = in_dir + "/" + file_id + "/2YN.dS"
			dS_file_h = open(dS_file, 'r')
			dS_lines = dS_file_h.readlines()

			# write output
			out_file = out_dir + "/" + file_id + ".rate"
			out_file_h = open(out_file, 'w')

			# add first seq_id to list
			seq_list = []
			dN_line = dN_lines[1]
			dN_words = dN_line.strip('\n').split()
			seq_list.append(dN_words[0])

			for line_index in range(2,len(dN_lines)-1):
				dN_line = dN_lines[line_index]
				dS_line = dS_lines[line_index]
				dN_words = dN_line.strip('\n').split()
				dS_words = dS_line.strip('\n').split()
				
				seq_1 = dN_words[0]
				seq_list.append(seq_1)
				for pairs_index in range(len(seq_list) - 1):
					seq_2 = seq_list[pairs_index]
					dN = float(dN_words[pairs_index + 1])
					dS = float(dS_words[pairs_index + 1])
					if dS == 0:
						dNdS = "NA"
					else:
						dNdS = round(dN / dS, 4)

					# write to output
					out_file_h.write(str(seq_1) + '\t' + str(seq_2) + '\t' + str(dNdS) + '\t' + str(dN) + '\t' + str(dS) + '\n')

			# close files
			out_file_h.close()
			dN_file_h.close()
			dS_file_h.close()


# print verbose
if verbose:
	print('count_in = %i' % (count_in))