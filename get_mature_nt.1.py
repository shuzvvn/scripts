#!/usr/bin/python3

# get_mature_nt.1.py

# Shu-Ting Cho <vivianlily6@hotmail.com>
# get mature nt seq according to list
# v1 2018/03/13

# Usage: python3 /home/shutingcho/pyscript/get_mature_nt.1.py --in_fasta=/scratch/shutingcho/phyto30/phyto30.12/source_data/cds/PLY.nt.fasta --in_list=/home/shutingcho/project/phyto30/phyto30.14/run3/effector/run3.list.2 --id_index=0 --pos_index=1 --out_file=/scratch/shutingcho/phyto21r/phyto21r.02/effector/mature.nt.1.fasta

import sys, getopt

# for reading fasta file
from Bio import SeqIO

# for reading csv
import csv

# get options
opts, args = getopt.getopt(sys.argv[1:], '', longopts=[
	'in_fasta=',
	'in_list=',
	'id_index=',
	'pos_index=',
	'out_file='])

# get variables from opts
for opt, arg in opts:
	if opt == "--in_fasta":
		ntfasta_dict = SeqIO.index(str(arg), 'fasta')
	elif opt == "--in_list":
		list_file_h = open(str(arg), 'r')
		in_list = list(csv.reader(list_file_h, delimiter='\t'))
		list_file_h.close()
	elif opt == "--id_index":
		id_index = int(arg)
	elif opt == "--pos_index":
		pos_index = int(arg)
	elif opt == "--out_file":
		out_file = str(arg)
	else:
		assert False, "unhandled option"

# make dir if not exist
import os
if out_file:
	out_dir = os.path.dirname(out_file)
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)

# function to cut nt seq base on SP cleavage cite position
count_in, count_out = 0, 0
print('#locus_tag\tcut_after:\tmature_nt')
with open(out_file, 'w') as out_file_h:
	for i in in_list:
		count_in += 1
		id_h = i[id_index]
		pos_h = int(i[pos_index])
		mature_start = pos_h * 3 - 2
		try:
			in_seq_h = str(ntfasta_dict[id_h].seq)
			out_seq_h = in_seq_h[mature_start-1:]
			out_file_h.write('>%s\n%s\n' % (id_h, out_seq_h))
			print('%s\t%i\t%i' % (id_h, mature_start-1, len(out_seq_h)))
			count_out += 1
		except:
			print(str(id_h) + ' was not found!!!')
print('count_in = %i, count_out = %i, count_diff = %i' % (count_in, count_out, count_in - count_out))
