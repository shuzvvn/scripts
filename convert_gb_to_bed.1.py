#!/usr/bin/env python
# encoding: utf-8
"""
bed_from_genbank.py
grab the gene records from a genbank file (edit for other record types).
- requires:  biopython
"""

from Bio import SeqIO

import pdb, sys, getopt

# Usage: python3 /home/shutingcho/pyscript/change_annotation.1.py --in_list=/home/shutingcho/project/phyto30/phyto30.14/run3/effector/run3.list.2 --in_list_index=0 --in_info=/home/shutingcho/project/phyto30/phyto30.12/info/PLY_v1.cds.03.info.ko.cog.desc.merg --in_info_index=0 --match_pattern=hypothetical protein --match_index=8 --change_to=putative effector --out_info=/home/shutingcho/project/phyto30/phyto30.12/info/PLY_v1.cds.03.info.ko.cog.desc.merg.1

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
		for feature in record.features:
			if feature.type == 'gene':
				start = feature.location.start.position
				stop = feature.location.end.position
				try:
					name = feature.qualifiers['gene'][0]
				except:
					# some features only have a locus tag
					name = feature.qualifiers['locus_tag'][0]
				if feature.strand < 0:
					strand = "-"
				else:
					strand = "+"
				bed_line = "record.name\t{0}\t{1}\t{2}\t1000\t{3}\t{0}\t{1}\t65,105,225\n".format(start, stop, name, strand)
				outf.write(bed_line)
outf.close()
