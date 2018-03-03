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
		vcf_reader = vcf.Reader(open(str(arg), 'r'))
	elif opt == "--cds_info":
		cds_info = str(arg)
	elif opt == "--nt_fasta":
		ntfasta_dict = SeqIO.to_dict(SeqIO.parse(str(arg), 'fasta'))
	elif opt == "--pseudo_info":
		pseudo_info = str(arg)
	elif opt == "--RNA_info":
		RNA_info = str(arg)
	elif opt == "--out_file":
		out_file = str(arg)
	else:
		assert False, "unhandled option"

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
		self.mut_type = []

	def add_feature_info(self, feature_type, locus_tag = '', name = '', product = ''): # CDS / pseudo / RNA / intergenic
		self.feature_type = feature_type
		self.locus_tag = locus_tag
		self.name = name
		self.product = product

	def add_vcf_type(self, vcf_type): # SNP / indel
		self.vcf_type = vcf_type

	def add_mut_type(self, synonymous = True): # if SNP, synonymous / non-synonymous
		if synonymous:
			self.mut_type.append('synonymous')
		else:
			self.mut_type.append('non-synonymous')

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

### read input ###
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
	RNA_ranges = get_range_from_info(RNA_info, 'RNA', 1, 3, 7)
	for locus in RNA_ranges:
		all_feature_dict[locus].update(RNA_ranges[locus])

# store characterized vcf record in list
out_vcf_record = []

# vcf characterization
for vcf_record_i in vcf_reader: # for each vcf record
	locus = str(vcf_record_i.CHROM)
	vcf_record_o = VcfOut(locus, vcf_record_i.POS, vcf_record_i.REF, vcf_record_i.ALT)
	intragenic_switch = False # init switch to False
	for feature_range in all_feature_dict[locus]: # for each range in the vcf locating locus 
		feature_record = all_feature_dict[locus][feature_range] # dict[locus][range(feature_start, feature_end)] = [feature_type, locus_tag, gene_name, product]
		#snp_in_range = 0
		if (vcf_record_i.POS in feature_range) and (not intragenic_switch):
			intragenic_switch = True # switch on
			#snp_in_range += 1			
			vcf_record_o.add_feature_info(feature_record[0], feature_record[1], feature_record[2], feature_record[3])
			if len(vcf_record_i.REF) == len(vcf_record_i.ALT[0]):
				#vcf_record_o.add_vcf_type('SNP')
				if feature_record[0] == 'CDS':
					ori_seq = str(ntfasta_dict[feature_record[1]])
					position = mutate_pos(feature_range, vcf_record_i.POS)
					for alt in vcf_record_i.ALT:
						new_seq = SNP_mutate(ori_seq, position, alt)
						vcf_record_o.add_mut_type(check_mutation(11, ori_seq, new_seq))
			elif len(vcf_record_i.REF) > len(vcf_record_i.ALT[0]):
				vcf_record_o.add_vcf_type('deletion')
			else:
				vcf_record_o.add_vcf_type('insertion')
		elif (vcf_record_i.POS in feature_range) and (intragenic_switch): # vcf in more than one feature
			print('WARNING: ' + str(vcf_record_o.CHROM) + ': ' + str(vcf_record_o.POS) + ' locates in more than one feature!!' )
	if not intragenic_switch: # this vcf is not in any gene => it's intergenic
		vcf_record_o.add_feature_info('intergenic')
	out_vcf_record.append(vcf_record_o)

# write output
count_vcf, count_cds, count_pseudo, count_RNA, count_intergenic = 0, 0, 0, 0, 0

out_file_h = open(out_file, 'w')
for i in out_vcf_record:
	count_vcf += 1
	vcf_record_list = [i.CHROM, i.POS, i.REF, i.ALT, i.vcf_type, i.feature_type, i.mut_type, i.locus_tag, i.name, i.product]
	out_file_h.write('\t'.join(vcf_record_list + '\n'))
	if i.feature_type == 'CDS':
		count_cds += 1
	elif i.feature_type == 'pseudo':
		count_pseudo += 1
	elif i.feature_type == 'RNA':
		count_RNA += 1
	elif i.feature_type == 'intergenic':
		count_intergenic += 1
	else:
		print('WTF is this:', str('\t'.join(vcf_record_list + '\n')))
out_file_h.close()

# print report
print('count_vcf = %i\ncount_in_cds = %i\ncount_in_pseudo = %i\ncount_in_RNAs = %i\ncount_in_intergenic = %i' % (count_vcf, count_cds, count_pseudo, count_RNA, count_intergenic))
# end of script
