#!/usr/bin/python3

# SNP_position.2.py

# Shu-Ting Cho <vivianlily6@hotmail.com>
# report SNP position, in inter/intragenic region
# v1 2018/01/11
# v2 2018/02/25
# in_list => vcf

# Usage: python /home/shutingcho/pyscript/SNP_position.1.py --in_list=/home/shutingcho/project/phyto39/phyto39.03/SNP.2.list --in_list_index=0 --cds_info=/scratch/shutingcho/phyto39/phyto39.03/source_data/cds/CP000061.info --pseudo_info=/scratch/shutingcho/phyto39/phyto39.03/source_data/pseudo/CP000061.info --other_info=/scratch/shutingcho/phyto39/phyto39.03/source_data/RNA/CP000061.RNA --out_file=/home/shutingcho/project/phyto39/phyto39.03/SNP.CP000061.out

import sys, getopt, re

# get options
opts, args = getopt.getopt(sys.argv[1:], '', longopts=[
	'vcf=',
	'cds_info=',
	'pseudo_info=',
	'other_info=',
	'out_file='])

# get variables from opts
for opt, arg in opts:
	if opt == "--vcf":
		vcf = str(arg)
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

# function to get dict from vcf
def get_dict_from_table(in_file, key_index = 0, value_index = 1):
	out_dict = {}
	key_index = int(key_index)
	value_index = int(value_index)
	with open(in_file, 'r') as in_file_h:
		for line in in_file_h:
			words = line.strip('\n').split('\t')
			try:
				out_dict[words[key_index]]
			except:
				out_dict[words[key_index]] = []
			out_dict[words[key_index]].append(words[value_index])
	return out_dict

# function to get dict
def get_range_from_info(in_file, feature_type, key_index = 1, start_index = 3, name_index = 7):
	out_dict = {}
	range_dict = {}
	with open(in_file, 'r') as in_file_h:
		for line in in_file_h:
			words = line.strip('\n').split('\t')
			try:
				out_dict[words[key_index]]
			except:
				out_dict[words[key_index]] = {} # set each locus as a dict
			out_dict[words[key_index]][range(int(words[start_index]), int(words[start_index+1])+1)] = [feature_type, words[0], words[name_index], words[name_index+1]]
	return out_dict

# get list of SNP
snp_dict = get_dict_from_table(vcf)

# get list of feature ranges
all_feature_dict = {}
if cds_info:
	cds_ranges = get_range_from_info(cds_info, 'CDS', 1, 3, 7)
	all_feature_dict = cds_ranges
if pseudo_info:
	pseudo_ranges = get_range_from_info(pseudo_info, 'psuedo', 1, 3, 7)
	for locus in pseudo_ranges:
		all_feature_dict[locus] = {**all_feature_dict[locus], **pseudo_ranges[locus]}
if other_info:
	other_ranges = get_range_from_info(other_info, 'RNA', 1, 3, 7)
	for locus in pseudo_ranges:
		all_feature_dict[locus] = {**all_feature_dict[locus], **pseudo_ranges[locus]}

'''

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
'''