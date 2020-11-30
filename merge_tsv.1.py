#!/usr/bin/python3

# merge_tsv.1.py

# Shu-Ting Cho <vivianlily6@hotmail.com>
# merge_tsv base on columne
# v1 2020/11/30

# Usage: 
#python3 /home/shu-ting/pyscript/merge_tsv.1.py --in_file_1=/scratch/shutingcho/phyto21/phyto21.21/source_data/cds/BGWL.info --in_file_2=/scratch/shutingcho/phyto21/phyto21.21/annotation/ortho_info/BGWL.txt --index_1=0 --index_2=2 --col_1="0,1,2,3,4,5,6,7,8" --col_2="0" --out_file=/scratch/shutingcho/phyto21/phyto21.21/annotation/phyto21.info.1 --print_empty="NA" --verbose=1


import sys, getopt

# for read tsv to list of lists
import csv


## set default values
verbose = True
print_empty = None


# get options
opts, args = getopt.getopt(sys.argv[1:], '', longopts=[
	'in_file_1=',
	'in_file_2=',
	'index_1=',
	'index_2=',
	'col_1=',
	'col_2=',
	'out_file=',
	'print_empty=',
	'verbose='])



# get variables from opts
for opt, arg in opts:
	if opt == "--in_file_1":
		try:
			# read csv file as a list of lists
			with open(str(arg), 'r') as read_obj:
				csv_reader = csv.reader(read_obj, delimiter="\t")
				list_of_in_file_1_h = list(csv_reader)
		except IOError:
			print('Cannot read in_file_1!')
	elif opt == "--in_file_2":
		try:
			# read csv file as a list of lists
			with open(str(arg), 'r') as read_obj:
				csv_reader = csv.reader(read_obj, delimiter="\t")
				list_of_in_file_2_h = list(csv_reader)
		except IOError:
			print('Cannot read in_file_2!')
	elif opt == "--index_1":
		try:
			index_1 = int(arg)
		except ValueError:
			print("Oops!  That was no valid number for index_1.  Try again...")
	elif opt == "--index_2":
		try:
			index_2 = int(arg)
		except ValueError:
			print("Oops!  That was no valid number for index_2.  Try again...")
	elif opt == "--col_1":
		col_1 = str(arg)
	elif opt == "--col_2":
		col_2 = str(arg)
	elif opt == "--out_file":
		out_file = str(arg)
	elif opt == "--print_empty":
		print_empty = str(arg)
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



# function to make dict using index content as key
def list2dict_index(in_list, ID, in_index = 0, col = "all"):
	out_dict = {}

	if col != "all":
		try:
			col_list = col.split(',')
		except:
			print("Cannot read col for " + ID)
	else:
		col_list = list(range(len(in_list[0])))#.remove(in_index)
	# read line to dict
	for line in in_list:
		# set key
		try:
			key_h = line[in_index]
		except:
			print("in_index out of range for " + ID)
		# set value
		value_h = []
		try:
			for i in col_list:
				i = int(i)
				value_h.append(line[i])
		except:
			print("col out of range for " + ID)
		out_dict[key_h] = value_h

	return out_dict



## main ##

# count lines in input file
count_in_1 = len(list_of_in_file_1_h)
count_in_2 = len(list_of_in_file_2_h)
# count lines in output file
count_out = 0


# read list to dict
dict_1 = list2dict_index(in_list = list_of_in_file_1_h, ID = "in_file_1", in_index = index_1, col = col_1)
dict_2 = list2dict_index(in_list = list_of_in_file_2_h, ID = "in_file_2", in_index = index_2, col = col_2)



# write output
out_file_h = open(out_file, 'w')

for key_h in dict_1:
	value_1_h = dict_1[key_h]
	# find key in dict_2
	if key_h in dict_2:
		value_2_h = dict_2[key_h]
	elif print_empty is not None:
		value_2_h = [print_empty] * len(col_2)
	else:
		value_2_h = [""] * len(col_2)
	# combine value and write to output
	out_file_h.write("\t".join(value_1_h) + "\t" + "\t".join(value_2_h) + "\n")
	# count output
	count_out += 1

out_file_h.close()


# print verbose
if verbose:
	print('count_in_1 = %i, count_in_2 = %i, count_out = %i' % (count_in_1, count_in_2, count_out))