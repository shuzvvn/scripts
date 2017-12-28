# open and read bpo
with open('/scratch/shutingcho/phyto38/phyto38.01/orthomcl/info/phyto+1others-.list') as file_phyto: 
	data = file_phyto.read() 
	l_phyto = data.split()
	with open('/scratch/shutingcho/phyto38/phyto38.01/bpo/phyto38.01.bpo') as file_bpo: 
		for line in file_bpo: 
			line = line.rstrip('\n')
			words = line.split(';')
			if (words[1] in l_phyto) and (words[3] not in l_phyto):
				print(line.replace(';', '\t'))
# end of script
