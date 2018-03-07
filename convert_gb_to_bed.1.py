#!/usr/bin/python3

# convert_gb_to_bed.1.py

# Shu-Ting Cho <vivianlily6@hotmail.com>
# grab the gene records from a genbank file
# requires:  biopython

# v1 2018/02/25


from Bio import SeqIO

import sys, getopt

# Usage: python3 /home/shutingcho/pyscript/convert_gb_to_bed.1.py --in_file=/home/shutingcho/project/phyto30/phyto30.12/tbl2asn/v2.gb --out_file=/home/shutingcho/project/phyto29/phyto29.03/v2.bed

opts, args = getopt.getopt(sys.argv[1:], '', longopts=[
	'in_file=', 
	'out_file='])

# get variables from opts
for opt, arg in opts:
	if opt == "--in_file":
		in_file = str(arg)
	elif opt == "--out_file":
		out_file = str(arg)
	else:
		assert False, "unhandled option"

outf = open(out_file, 'w')
for record in SeqIO.parse(open(in_file, "rU"), "genbank") :
	locus = record.name
	for feature in record.features:
		if feature.type != 'source':
			start = feature.location.start.position
			stop = feature.location.end.position
			try:
				name = feature.qualifiers['locus_tag'][0] + ':' + feature.qualifiers['gene'][0]
			except:
				# some features only have a locus tag
				name = feature.qualifiers['locus_tag'][0]
			if feature.strand < 0:
				strand = "-"
			else:
				strand = "+"
			try:
				if feature.qualifiers['pseudo']:
					colorcode = '153,153,153'
			except:
				if feature.type == 'tRNA':
					colorcode = '51,102,255'
				elif feature.type == 'rRNA':
					colorcode = '51,204,51'
				elif feature.type == 'CDS':
					colorcode = '0,51,102'
				else:
					colorcode = False
			if colorcode:
				bed_line = "\t{0}\t{1}\t{2}\t1000\t{3}\t{0}\t{1}\t{4}\n".format(start, stop, name, strand, colorcode)
				outf.write(locus + bed_line)
outf.close()
