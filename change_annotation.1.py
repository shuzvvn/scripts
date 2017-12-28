#!/usr/bin/python3

import sys, getopt

# Usage: python3 make_bed.1.py --in_list=/scratch/shutingcho/phyto38/phyto38.03/phyto38.03.length --out_file=/scratch/shutingcho/phyto38/phyto38.03/phyto38.03.w200.bed --win_size=250

opts, args = getopt.getopt(sys.argv[1:], '', longopts=['in_list=', 'in_info=', 'out_info='])

# get variables from opts
for opt, arg in opts:
	if opt == "--in_list":
		inList = arg
	elif opt == "--in_info":
		inInfo = arg
	elif opt == "--out_info":
		outFile = arg
	else:
		assert False, "unhandled option"

candidates=[]

with open(inList) as inList_h:
    for line in inList_h:
        words = line.split('\t')
        candidates.append(words[0])

outFile_h = open(outFile, 'w')

with open(inInfo) as inInfo_h:
	for line in inInfo_h:
		line = line.rstrip('\n')
		words = line.split('\t')
		if words[0] in candidates:
            if words[8] == 'hypothetical protein':
                words[8] = 'putative effector'
            else:
                print('\t'.join(words))



                
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
