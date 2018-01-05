#!/usr/bin/python3

import sys, getopt, re

# Usage: python3 /home/shutingcho/pyscript/change_annotation.1.py --in_list=/home/shutingcho/project/phyto30/phyto30.14/run3/effector/run3.list.2 --in_list_index=0 --in_info=/home/shutingcho/project/phyto30/phyto30.12/info/PLY_v1.cds.03.info.ko.cog.desc.merg --in_info_index=0 --match_pattern=hypothetical protein --match_index=8 --change_to=putative effector --out_info=/home/shutingcho/project/phyto30/phyto30.12/info/PLY_v1.cds.03.info.ko.cog.desc.merg.1

opts, args = getopt.getopt(sys.argv[1:], '', longopts=[
	'in_list=', 
	'in_list_index=', 
	'in_info=', 
	'in_info_index=', 
	'match_pattern=',
	'match_index=',
	'change_to=',
	'out_info='])

# get variables from opts
for opt, arg in opts:
	if opt == "--in_list":
		in_list = str(arg)
	elif opt == "--in_list_index":
		in_list_index = int(arg)
	elif opt == "--in_info":
		in_file = str(arg)
	elif opt == "--in_info_index":
		in_file_index = int(arg)
	elif opt == "--match_pattern":
		match_pattern = str(arg)
	elif opt == "--match_index":
		match_index = int(arg)
	elif opt == "--change_to":
		change_to = str(arg)
	elif opt == "--out_info":
		out_file = str(arg)
	else:
		assert False, "unhandled option"

# count number for record
changed_num = 0
unchanged_num = 0

# get list of in candidates
candidates = []
with open(in_list) as in_list_h:
	for line in in_list_h:
		words = line.split('\t')
		candidates.append(words[in_list_index])
candidate_num = len(candidates)

out_file_h = open(out_file, 'w')

with open(in_file) as in_file_h:
	for line in in_file_h:
		words = line.split('\t')
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
print('count_in = ', str(candidate_num), ', count_changed = ', str(changed_num), ', count_unchanged = ', str(unchanged_num))
# end of script
