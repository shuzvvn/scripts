#!/usr/bin/python3

import sys, getopt

# Usage: python3 make_bed.1.py --in_file=/scratch/shutingcho/phyto38/phyto38.03/phyto38.03.length --out_file=/scratch/shutingcho/phyto38/phyto38.03/phyto38.03.w200.bed --win_size=250

opts, args = getopt.getopt(sys.argv[1:], '', longopts=['in_file=', 'out_file=', 'win_size='])

# get variables from opts
for opt, arg in opts:
	if opt == "--in_file":
		inLgthLs = arg
	elif opt == "--out_file":
		outFile = arg
	elif opt == "--win_size":
		setWinSize = int(arg)
	else:
		assert False, "unhandled option"

outFile_h = open(outFile, 'w')
with open(inLgthLs) as inFile_h:
	lastWinID = 0
	for line in inFile_h:
		line = line.rstrip('\n')
		words = line.split('\t')
		chrID = str(words[0])
		lgth = int(words[1])
		i = 0
		while i < lgth:
			start = i + 1
			windowID = lastWinID + i / setWinSize + 1
			if i < (lgth - setWinSize):
				end = i + setWinSize
			else:
				end = lgth
			winSize = end - start + 1
			outFile_h.write('%s\t%s\t%s\t.\t.\t+\t%s\t%s\n' % (str(chrID), str(start), str(end), str(int(windowID)), str(winSize)))
			i += setWinSize
		lastWinID = windowID
outFile_h.close()
# end of script
