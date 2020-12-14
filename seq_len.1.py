#!/usr/bin/python3

# seq_len.1.py

# Shu-Ting Cho <vivianlily6@hotmail.com>
# seq_len for fastq in a dir
# v1 2020/12/13

# Usage: 
# python3 /home/shu-ting/pyscript/seq_len.1.py --in_dir=/home/shu-ting/scratch/raw_reads/Nanopore_20201213/CD07748/ --print_each=0


import sys, getopt

from Bio import SeqIO
import os


## set default values
# verbose = True
print_each = False


# get options
opts, args = getopt.getopt(sys.argv[1:], '', longopts=[
	# 'in_file_1=',
	# 'in_file_2=',
	# 'index_1=',
	# 'index_2=',
	# 'col_1=',
	# 'col_2=',
	# 'out_file=',
	'print_each=',
	'in_dir='])



# get variables from opts
for opt, arg in opts:
	if opt == "--in_dir":
		in_dir = str(arg).rstrip("/")
		try:
			list_files = []
			for file in os.listdir(in_dir):
				if file.endswith(".fastq"):
					list_files.append(in_dir + '/' + file)
		except IOError:
			print('Cannot read in_dir!')
	# elif opt == "--in_file_2":
	# 	try:
	# 		# read csv file as a list of lists
	# 		with open(str(arg), 'r') as read_obj:
	# 			csv_reader = csv.reader(read_obj, delimiter="\t")
	# 			list_of_in_file_2_h = list(csv_reader)
	# 	except IOError:
	# 		print('Cannot read in_file_2!')
	# elif opt == "--index_1":
	# 	try:
	# 		index_1 = int(arg)
	# 	except ValueError:
	# 		print("Oops!  That was no valid number for index_1.  Try again...")
	# elif opt == "--index_2":
	# 	try:
	# 		index_2 = int(arg)
	# 	except ValueError:
	# 		print("Oops!  That was no valid number for index_2.  Try again...")
	# elif opt == "--col_1":
	# 	col_1 = str(arg)
	# elif opt == "--col_2":
	# 	col_2 = str(arg)
	# elif opt == "--out_file":
	# 	out_file = str(arg)
	elif opt == "--print_each":
		try:
			print_each = bool(int(arg))
		except ValueError:
			print("value for print_each should be a boolean")
	# elif opt == "--verbose":
	# 	try:
	# 		verbose = bool(arg)
	# 	except ValueError:
	# 		print("value for verbose should be a boolean")
	else:
		assert False, "unhandled option"


# Calculate N50 for a sequence of numbers.
    # Args:
    #     list_of_lengths (list): List of numbers.
 
    # Returns:
    #     float: N50 value.
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


for in_fastq in list_files:
	file_name = os.path.basename(in_fastq)
	records = SeqIO.index(in_fastq, "fastq")
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