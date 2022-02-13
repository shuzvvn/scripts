#!/usr/bin/python3

# seq_len.2.py

# Shu-Ting Cho <vivianlily6@hotmail.com>
# seq_len for fastq/fasta in a dir

# v2 2022/02/12 input accept both fastq and fasta
# v1 2020/12/13

# Usage: 
# python3 /home/shu-ting/pyscript/seq_len.1.py --in_dir=/home/shu-ting/scratch/raw_reads/Nanopore_20201213/CD07748/ --print_each=0


import sys, getopt

from Bio import SeqIO
import os


## set default values
# verbose = True
print_each = False
in_format = 'fasta'


# get options
opts, args = getopt.getopt(sys.argv[1:], '', longopts=[
	'in_format=',
	'print_each=',
	'in_dir='])



# get variables from opts
for opt, arg in opts:
	if opt == "--in_format":
		in_format = str(arg)
	elif opt == "--in_dir":
		in_dir = str(arg).rstrip("/")
		try:
			list_files = []
			for file in os.listdir(in_dir):
				if file.endswith("." + in_format):
					list_files.append(in_dir + '/' + file)
		except IOError:
			print('Cannot read in_dir!')
	elif opt == "--print_each":
		try:
			print_each = bool(int(arg))
		except ValueError:
			print("value for print_each should be a boolean")
	else:
		assert False, "unhandled option"


# Calculate N50 for a sequence of numbers.
def calculate_N50(list_of_lengths):
    tmp = []
    for tmp_number in set(list_of_lengths):
            tmp += [tmp_number] * list_of_lengths.count(tmp_number) * tmp_number
    tmp.sort()
 
    if (len(tmp) % 2) == 0:
        median = (tmp[int(len(tmp) / 2) - 1] + tmp[int(len(tmp) / 2)]) / 2
    else:
        median = tmp[int(len(tmp) / 2)]
 
    return median


## main ##

dir_length = 0
dir_count = 0
list_dir_length = []


for in_file in list_files:
	file_name = os.path.basename(in_file)
	records = SeqIO.index(in_file, in_format)
	file_count = len(records)
	file_length = 0
	list_file_length = []
	for record in records:
		seq_length = len(records[record])
		file_length += seq_length
		list_file_length.append(seq_length)

	# output for each file
	if print_each is True:
		print(file_name + '\t' + str(file_count) + '\t' + str(file_length) + '\t' + str(calculate_N50(list_file_length)))
	
	# count for whole dir
	dir_length += file_length
	dir_count += file_count
	list_dir_length = list_dir_length + list_file_length

	records.close()

print(os.path.basename(in_dir) + '\t' + str(dir_count) + '\t' + str(dir_length) + '\t' + str(calculate_N50(list_dir_length)))