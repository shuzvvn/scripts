#!/usr/bin/python3

# SNP_position.3.py

# Shu-Ting Cho <vivianlily6@hotmail.com>
# report SNP position, in inter/intragenic region
# v1 2018/01/11
# v2 2018/02/26: in_list => vcf; include genome(locus) in the record
# v3 2018/03/01: import vcf module; define nonsynonymous SNP

# Usage: python /home/shutingcho/pyscript/SNP_position.2.py --in_vcf=/home/shutingcho/project/phyto39/phyto39.03/SNP.2.list --cds_info=/scratch/shutingcho/phyto39/phyto39.03/source_data/cds/CP000061.info --ntfasta=/--pseudo_info=/scratch/shutingcho/phyto39/phyto39.03/source_data/pseudo/CP000061.info --RNA_info=/scratch/shutingcho/phyto39/phyto39.03/source_data/RNA/CP000061.RNA --out_file=/home/shutingcho/project/phyto39/phyto39.03/SNP.CP000061.out

import sys, getopt, re, vcf

# for dna translation
from Bio.seq import Seq
from Bio.Alphabet import IUPAC

# for reading fasta file
from Bio import SeqIO

# get options
opts, args = getopt.getopt(sys.argv[1:], '', longopts=[
	'in_vcf=',
	'cds_info=',
	'nt_fasta=',
	'pseudo_info=',
	'RNA_info=',
	'out_file='])

# get variables from opts
for opt, arg in opts:
	if opt == "--in_vcf":
		in_vcf = str(arg)
	elif opt == "--cds_info":
		cds_info = str(arg)
	elif opt == "--nt_fasta":
		nt_fasta = str(arg)
	elif opt == "--pseudo_info":
		pseudo_info = str(arg)
	elif opt == "--RNA_info":
		RNA_info = str(arg)
	elif opt == "--out_file":
		out_file = str(arg)
	else:
		assert False, "unhandled option"

'''
# function to get list
def get_list_from_table(in_file, index_0 = 0, index_1 = 1):
	out_list = []
	with open(in_file, 'r') as in_file_h:
		for line in in_file_h:
			words = line.strip('\n').split('\t')
			out_list.append((words[int(index_0)], words[int(index_1)]))
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
			out_dict[words[key_index]].append(int(words[value_index]))
	return out_dict
'''

# function to get dict
def get_range_from_info(in_file, feature_type, key_index = 1, start_index = 3, name_index = 7):
	out_dict = {}
	range_dict = {}
	with open(in_file, 'r') as in_file_h:
		for line in in_file_h:
			words = line.strip('\n').split('\t')
			if not out_dict[words[key_index]]: 
				out_dict[words[key_index]] = {} # set each locus as a dict
			# key is feature range, value is feature info
			# dict[locus][range(feature_start, feature_end)] = [feature_type, locus_tag, gene_name, product]
			out_dict[words[key_index]][range(int(words[start_index]), int(words[start_index+1])+1)] = [feature_type, words[0], words[name_index], words[name_index+1]]
	return out_dict

# class to record the vcf characterization record
class VcfOut:
	def __init__(self, locus, position, REF, ALT):
		self.CHROM = locus
		self.POS = position
		self.REF = REF
		self.ALT = ALT
		# set defalt values
		self.feature_type = 'intergenic'
		self.vcf_type = 'SNP'
		self.synonymous = True

	def add_feature_info(self, feature_type, locus_tag = '', name = '', product = ''): # CDS / pseudo / RNA / intergenic
		self.feature_type = feature_type
		self.locus_tag = locus_tag
		self.name = name
		self.product = product

	def add_vcf_type(self, vcf_type): # SNP / indel
		self.vcf_type = vcf_type

	def nonsynonymous(self): # if SNP, synonymous / non-synonymous
		self.synonymous = False

# get mutating postion in feature
def mutate_pos(feature_range, SNP_position):
	position = int(SNP_position) - feature_range[0] + 1
	return position

# get mut_coding_dna
def SNP_mutate(nt_seq, position, alt): # position start from 1
	new_seq = nt_seq[:int(position)-1] + alt + nt_seq[int(position):]
	return new_seq

# check non-synonymous mutation
def check_mutation(transl_table = 11, nt1, nt2):
	genetic_code = int(transl_table)
	coding_dna_1 = Seq(nt1, IUPAC.unambiguous_dna)
	coding_dna_2 = Seq(nt2, IUPAC.unambiguous_dna)
	if str(coding_dna_1.translate(table = genetic_code)) == str(coding_dna_2.translate(table = genetic_code)):
		return True
	else:
		return False










### read input
# read nt_fasta to dict
record_dict = SeqIO.to_dict(SeqIO.parse(nt_fasta, 'fasta'))





# get list of SNP
snp_list = get_list_from_table(in_vcf)
snp_dict = get_dict_from_table(in_vcf)

# get list of feature ranges
all_feature_dict = {}
if cds_info:
	cds_ranges = get_range_from_info(cds_info, 'CDS', 1, 3, 7)
	all_feature_dict = cds_ranges.copy()
if pseudo_info:
	pseudo_ranges = get_range_from_info(pseudo_info, 'pseudo', 1, 3, 7)
	for locus in pseudo_ranges:
		all_feature_dict[locus].update(pseudo_ranges[locus])
if RNA_info:
	RNA_ranges = get_range_from_info(RNA_info, 'RNAs', 1, 3, 7)
	for locus in RNA_ranges:
		all_feature_dict[locus].update(RNA_ranges[locus])

# store characterized vcf record in list
out_vcf_record = []

# read vcf file
vcf_reader = vcf.Reader(open(in_vcf, 'r'))
for record_i in vcf_reader: # for each vcf record
	locus = str(record_i.CHROM)
	record_o = VcfOut(locus, record_i.POS, record_i.REF, record_i.ALT)
	intragenic_switch = False
	for feature_range in all_feature_dict[locus]: # for each range in the vcf locating locus 
		feature_record = all_feature_dict[locus][feature_range]
		#snp_in_range = 0
		if record_i.POS in feature_range:
			intragenic_switch = True
			#snp_in_range += 1
			#dict[locus][range(feature_start, feature_end)] = [feature_type, locus_tag, gene_name, product]
			record_o.add_feature_info(feature_record[0], feature_record[1], feature_record[2], feature_record[3])
			if len(record_i.REF) == len(record_i.ALT):
				record_o.add_vcf_type('SNP')
			elif len(record_i.REF) > len(record_i.ALT)::
				record_o.add_vcf_type('deletion')
			else:
				record_o.add_vcf_type('insertion')
			############################## test synonymous ################################






			record_o.add_feature(feature_record[0])
			all_feature_dict[locus][feature_range]
	if not intragenic_switch: # this vcf is not in any gene => it's intergenic
		record_o.add_feature_info('intergenic')





	def add_feature_info(self, feature_type, locus_tag = '', name = '', product = ''): # CDS / pseudo / RNA / intergenic
		self.feature_type = feature_type
		self.locus_tag = locus_tag
		self.name = name
		self.product = product

	def add_vcf_type(self, vcf_type): # SNP / indel
		self.vcf_type = vcf_type

	def nonsynonymous(self): # if SNP, synonymous / non-synonymous
		self.synonymous = False







# store result in a dict
out_dict = {}
# intragenic SNP
intragenic_dict = {}
for locus in all_feature_dict:
	out_dict[locus] = {}
	intragenic_dict[locus] = []
	for feature_range in all_feature_dict[locus]:
		snp_in_range = 0
		for i in snp_dict[locus]:
			if int(i) in feature_range:
				snp_in_range += 1
				out_dict[locus][int(i)] = all_feature_dict[locus][feature_range]
				intragenic_dict[locus].append(int(i))
			elif snp_in_range != 0: # next snp not in the same feature
				break # jump to next feature range

# intergenic SNP
intergenic_dict = {}
for locus in snp_dict:
	intergenic_dict[locus] = [x for x in snp_dict[locus] if x not in intragenic_dict[locus]]
try:
	for locus in intergenic_dict:
		for i in intergenic_dict[locus]:
			out_dict[locus][int(i)] = ['intergenic']
except:
	pass

# write output
count_snp, count_cds, count_pseudo, count_RNAs, count_intergenic = 0, 0, 0, 0, 0

out_file_h = open(out_file, 'w')
for i in snp_list:
	count_snp += 1
	out_file_h.write(str(i[0]) + '\t' + str(i[1]) + '\t' + '\t'.join(out_dict[i[0]][int(i[1])]) + '\n')
	feature_type = out_dict[i[0]][int(i[1])][0]
	if feature_type == 'CDS':
		count_cds += 1
	elif feature_type == 'pseudo':
		count_pseudo += 1
	elif feature_type == 'RNAs':
		count_RNAs += 1
	elif feature_type == 'intergenic':
		count_intergenic += 1
	else:
		print('WTF is this:', str(i[0]), str(i[1]), '\t'.join(out_dict[i[0]][int(i[1])]))
out_file_h.close()

# print report
print('count_SNP = %i\ncount_in_cds = %i\ncount_in_pseudo = %i\ncount_in_RNAs = %i\ncount_in_intergenic = %i' % (count_snp, count_cds, count_pseudo, count_RNAs, count_intergenic))
# end of script


'''
