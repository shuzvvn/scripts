#!/usr/bin/python3

# SNP_position.1.py

# Shu-Ting Cho <vivianlily6@hotmail.com>
# report SNP position, in inter/intragenic region
# v1 2018/01/11

# Usage: python3 /home/shutingcho/pyscript/SNP_position.1.py python3 /home/shutingcho/pyscript/SNP_position.1.py --in_list=/home/shutingcho/project/phyto39/phyto39.03/SNP.2.list --in_list_index=0 --cds_info=/scratch/shutingcho/phyto39/phyto39.03/source_data/cds/CP000061.info --pseudo_info=/scratch/shutingcho/phyto39/phyto39.03/source_data/pseudo/CP000061.info --other_info=/scratch/shutingcho/phyto39/phyto39.03/source_data/RNA/CP000061.RNA --out_file=/home/shutingcho/project/phyto39/phyto39.03/SNP.CP000061.out

import sys, getopt, re

# get options
opts, args = getopt.getopt(sys.argv[1:], '', longopts=[
	'in_list=',
	'in_list_index=',
	'cds_info=',
	'pseudo_info=',
	'other_info=',
	'out_file='])

# get variables from opts
for opt, arg in opts:
	if opt == "--in_list":
		in_list = str(arg)
	elif opt == "--in_list_index":
		in_list_index = int(arg)
	elif opt == "--cds_info":
		cds_info = str(arg)
	elif opt == "--pseudo_info":
		pseudo_info = str(arg)
	elif opt == "--other_info":
		other_info = str(arg)
	elif opt == "--out_file":
		out_file = str(arg)
	else:
		assert False, "unhandled option"

# function to get list
def get_list_from_table(in_file, in_index = 0):
	out_list = []
	in_index = int(in_index)
	with open(in_file, 'r') as in_file_h:
		for line in in_file_h:
			words = line.strip('\n').split('\t')
			out_list.append(words[in_index])
	return out_list

# function to get dict
def get_range_from_info(in_file, start_index, name_index):
	out_dict = {}
	with open(in_file, 'r') as in_file_h:
		for line in in_file_h:
			words = line.strip('\n').split('\t')
			out_dict[range(int(words[start_index]), int(words[start_index+1])+1)] = [words[0], words[name_index], words[name_index+1]]
	return out_dict

# get list of SNP
snp_list = get_list_from_table(in_list, in_list_index)

# get list of feature ranges
feature_lists = []
if cds_info:
	cds_ranges = get_range_from_info(cds_info, 3, 7)
	feature_lists.append(cds_ranges)
if pseudo_info:
	pseudo_ranges = get_range_from_info(pseudo_info, 3, 7)
	feature_lists.append(pseudo_ranges)
if other_info:
	other_ranges = get_range_from_info(other_info, 1, 4)
	feature_lists.append(other_ranges)

# store result in a dict
out_dict = {}

# store count for SNP in each feature
count_cds = 0
count_pseudo = 0
count_others = 0
intragenic_list = []

# intragenic SNP
for feature_ranges in feature_lists: # each info dict
	for feature_range in feature_ranges: # each range in dict
		snp_in_range = 0 # init snp in this range
		for i in snp_list:
			i = int(i)
			if i in feature_range:
				snp_in_range += 1
				if feature_ranges == cds_ranges:
					feature_type = 'cds'
					count_cds += 1
				elif feature_ranges == pseudo_ranges:
					feature_type = 'pseudo'
					count_pseudo += 1
				else:
					feature_type = 'others'
					count_others += 1
				out_dict[i] = [feature_type] + feature_ranges[feature_range]
				intragenic_list.append(str(i))
			else:
				if snp_in_range != 0:
					break

# intergenic SNP
intergenic_list = [x for x in snp_list if x not in intragenic_list]
if intergenic_list:
	count_intergenic = len(intergenic_list)
	feature_type = 'intergenic'
	for i in intergenic_list:
		i = int(i)
		out_dict[i] = [feature_type]

# write output
out_file_h = open(out_file, 'w')
for i in snp_list:
	out_file_h.write(i + '\t' + '\t'.join(out_dict[int(i)]) + '\n')
out_file_h.close()

# print report
print('count_SNP = %i\ncount_in_cds = %i\ncount_in_pseudo = %i\ncount_in_others = %i\ncount_in_intergenic = %i' % (len(snp_list), count_cds, count_pseudo, count_others, count_intergenic))

# end of script
