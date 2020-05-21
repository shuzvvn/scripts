#!/usr/bin/python3

# convert_info_to_circos.1.py

# Shu-Ting Cho <vivianlily6@hotmail.com>
# convert info file to circos input: region / color table
# v1 2020/05/21


# Usage: python3 /home/shutingcho/pyscript/convert_info_to_circos.1.py --info=/home/shutingcho/project/phyto57/phyto57.11/info/phyto57.info --index_CHROM=1 --index_start=3 --index_end=4 --index_cog=12 --index_strand=6 --color_table=/home/shutingcho/project/phyto21/phyto21.36/color.txt --out_files=/scratch/shutingcho/phyto57/phyto57.11/circos/phyto57.info


import sys, getopt, os


# get options
opts, args = getopt.getopt(sys.argv[1:], '', longopts=[
	'info=',
	'index_CHROM=',
	'index_start=',
	'index_end=',
	'index_cog=',
	'index_strand=',
	'color_table=',
	'out_files='])

# get variables from opts
for opt, arg in opts:
	if opt == "--info":
		info = str(arg)
	elif opt == "--index_CHROM":
		try:
			index_CHROM = int(arg)
		except ValueError:
			print("Oops!  That was no valid number for index_CHROM.  Try again...")
	elif opt == "--index_start":
		try:
			index_start = int(arg)
		except ValueError:
			print("Oops!  That was no valid number for index_start.  Try again...")
	elif opt == "--index_end":
		try:
			index_end = int(arg)
		except ValueError:
			print("Oops!  That was no valid number for index_end.  Try again...")
	elif opt == "--index_cog":
		try:
			index_cog = int(arg)
		except ValueError:
			print("Oops!  That was no valid number for index_cog.  Try again...")
	elif opt == "--index_strand":
		try:
			index_strand = int(arg)
		except ValueError:
			print("Oops!  That was no valid number for index_strand.  Try again...")
	elif opt == "--color_table":
		color_table = str(arg)
	elif opt == "--out_files":
		out_files = str(arg)
		out_dir = os.path.dirname(out_files)
		if not os.path.exists(out_dir): # make dir if not exist
			os.makedirs(out_dir)
	else:
		assert False, "unhandled option"


# get color_dict from color_table
color_dict = {}
with open(color_table, 'r') as in_file_h:
	for line in in_file_h:
		words = line.strip('\n').split('\t')
		# key is COG, value is color code
		color_dict[str(words[0])] = str(words[1])


# name out file

### read input and write output ###
out_file_f_h = open(out_files + '.+', 'w')
out_file_r_h = open(out_files + '.-', 'w')

# count line
in_line = 0
out_f_line = 0
out_r_line = 0

with open(info, 'r') as in_file_h:
	for line in in_file_h:
		in_line = in_line + 1
		words = line.strip('\n').split('\t')
		chrom_id = str(words[index_CHROM])
		start = str(words[index_start])
		end = str(words[index_end])
		strand = int(words[index_strand])
		cog = str(words[index_cog])
		try:
			out_line = '\t'.join([chrom_id, start, end, 'color=' + color_dict[cog]])
		except KeyError:
			out_line = '\t'.join([chrom_id, start, end, 'color=' + "black"])
		if strand == 1:
			out_file_f_h.write(out_line + '\n')
			out_f_line = out_f_line + 1
		else:
			out_file_r_h.write(out_line + '\n')
			out_r_line = out_r_line + 1

out_file_f_h.close()
out_file_r_h.close()

# print report
print('in_line = %i\nout_forward = %i\nout_reverse = %i\nout_total = %i' % (in_line, out_f_line, out_r_line, out_f_line + out_r_line))
# end of script
