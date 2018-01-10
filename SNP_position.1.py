#!/usr/bin/python3

import sys, getopt, re

# Usage: python3 /home/shutingcho/pyscript/change_annotation.1.py --in_list=/home/shutingcho/project/phyto30/phyto30.14/run3/effector/run3.list.2 --in_list_index=0 --in_info=/home/shutingcho/project/phyto30/phyto30.12/info/PLY_v1.cds.03.info.ko.cog.desc.merg --in_info_index=0 --match_pattern=hypothetical protein --match_index=8 --change_to=putative effector --out_info=/home/shutingcho/project/phyto30/phyto30.12/info/PLY_v1.cds.03.info.ko.cog.desc.merg.1

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

# make function to get list
def get_list_from_table(in_table, in_index = 0, get_range = False):
	out_list = []
	in_index = int(in_index)
	with open(in_table, 'r') as in_table_h:
		for line in in_table_h:
			words = line.strip('\n').split('\t')
			if get_range:
				out_list.append({
					locus_tag : words[0], 
					feature_range : range(int(words[3]), int(words[4])+1),
					feature_name : ': '.join(words[7:9])})
			else:
				out_list.append(words[in_index])
	return out_list

# get list of SNP
snp_list = get_list_from_table(in_list, in_list_index)

# get list of feature ranges
if cds_info:
	cds_range_list = get_list_from_table(cds_info, get_range = True)
if pseudo_info:
	pseudo_range_list = get_list_from_table(pseudo_info, get_range = True)
if other_info:
	other_range_list = get_list_from_table(other_info, get_range = True)

# write output
out_file_h = open(out_file, 'w')

for i in snp_list:
	if i in 

		if words[in_file_index] in candidates:
			if re.match(match_pattern, words[match_index]): # re.match(pattern, string)
				words[match_index] = change_to
				changed_num += 1
				out_file_h.write('\t'.join(words))
			else:
				unchanged_num += 1
				out_file_h.write(line)
				print(words[in_file_index], words[match_index])
			candidates.remove(words[in_file_index])
		else:
			out_file_h.write(line)
out_file_h.close()

# print candidates not in in_file
if candidates:
	print(str(len(candidates)), 'items not in in_info: ') 
	print('\n'.join(candidates))

# print report
print('count_in = ', str(count_candidate), ', count_changed = ', str(changed_num), ', count_unchanged = ', str(unchanged_num))
# end of script



count_candidate = len(candidates)
