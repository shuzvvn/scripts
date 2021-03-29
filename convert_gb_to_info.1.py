#!/usr/bin/python3

# convert_gb_to_info.1.py

# Shu-Ting Cho <vivianlily6@hotmail.com>
# grab the gene records from a genbank file
# requires:  biopython
# v1 2021/03/28

# Usage:
# python3 /home/shutingcho/pyscript/convert_gb_to_info.1.py --in_file=/home/shutingcho/project/phyto30/phyto30.12/tbl2asn/v2.gb --out_file=/home/shutingcho/project/phyto29/phyto29.03/v2.info

## import modules
from Bio import SeqIO
import sys, getopt

## set default values
verbose = 1


## read arguments
# get options
opts, args = getopt.getopt(sys.argv[1:], '', longopts=[
	'in_file=', 
	'out_file=',
	'verbose='])

# get variables from opts
for opt, arg in opts:
	if opt == "--in_file":
		in_file = str(arg)
	elif opt == "--out_file":
		out_file = str(arg)
	elif opt == "--verbose":
		verbose = int(arg)
	else:
		assert False, "unhandled option"

verbose = bool(verbose)

# make dir if not exist
import os
out_dir = os.path.dirname(out_file)
if not os.path.exists(out_dir):
	os.makedirs(out_dir)




## main ##

# count records
count_CDS = 0
count_rRNA = 0
count_tRNA = 0
count_pseudo = 0
count_others = 0


# write output
out_file_h = open(out_file, 'w')

for record in SeqIO.parse(open(in_file, "r"), "genbank") :
	locus = record.name
	for feature in record.features:
		if feature.type != 'source':
			start = feature.location.start.position + 1
			stop = feature.location.end.position
			length = int(stop) - int(start) + 1
			# locus_tag
			try:
				locus_tag = feature.qualifiers['locus_tag'][0] 
			except:
				locus_tag = "NA"
			# gene name
			try:
				gene = feature.qualifiers['gene'][0]
			except:
				gene = "NA"
			# product description
			try:
				product = feature.qualifiers['product'][0]
			except:
				product = "NA"
			# strand
			strand = feature.strand

			# gene type
			if feature.type == 'tRNA':
				gene_type = 'tRNA'
				count_tRNA += 1
			elif feature.type == 'rRNA':
				gene_type = 'rRNA'
				count_rRNA += 1
			elif feature.type == 'CDS':
				gene_type = 'CDS'
				count_CDS += 1
			else:
				gene_type = False
				count_others += 1

			# pseudo
			try:
				if feature.qualifiers['pseudo']:
					pseudo = 'pseudo'
					count_pseudo += 1
			except:
				pseudo = 'NA'

			# write to output
			if gene_type:
				info_line = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n".format(locus_tag, locus, gene_type, start, stop, length, strand, gene, product, pseudo)
				out_file_h.write(info_line)

out_file_h.close()


# print verbose
if verbose:
	print('CDS = %i, rRNA = %i, tRNA = %i, pseudo = %i, others = %i' % (count_CDS, count_rRNA, count_tRNA, count_pseudo, count_others))