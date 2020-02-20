#!/usr/bin/python3

# SNP_position.3.py

# Shu-Ting Cho <vivianlily6@hotmail.com>
# report SNP position, in inter/intragenic region
# v1 2018/01/11
# v2 2018/02/26: in_list => vcf; include genome(locus) in the record
# v3 2018/03/05: import vcf,biopython module; define nonsynonymous SNP
# v4 2020/02/19: in_list to read lastz diff output

# Usage: python3 /mnt/c/Users/vvn/pyscript/SNP_position.3.py --in_vcf=/mnt/c/Users/vvn/project/phyto29/phyto29.03/bwa.AS280.vcf --cds_info=/mnt/c/Users/vvn/project/phyto29/phyto29.03/PLY_v1.cds.03.2.info.ko.cog.desc.merg --nt_fasta=/mnt/c/Users/vvn/project/phyto29/phyto29.03/v2.nt.fasta --transl_table=11 --pseudo_info=/mnt/c/Users/vvn/project/phyto29/phyto29.03/PLY_v1.pseudo.02.info.ko.cog.desc.merg --RNA_info=/mnt/c/Users/vvn/project/phyto29/phyto29.03/PLY_v1.RNA.info --out_file=/mnt/c/Users/vvn/project/phyto29/phyto29.03/SNP.PLY.out

import sys, getopt, re, vcf

# for dna translation
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

# for reading fasta file
from Bio import SeqIO

# get options
opts, args = getopt.getopt(sys.argv[1:], '', longopts=[
	'in_vcf=',
	'in_list=',
	'index_CHROM=',
	'index_POS=',
	'index_REF=',
	'index_ALT=',
	'transl_table=',
	'cds_info=',
	'nt_fasta=',
	'pseudo_info=',
	'RNA_info=',
	'out_file='])

# get variables from opts
for opt, arg in opts:
	if opt == "--in_vcf":
		try:
			vcf_reader = vcf.Reader(open(str(arg), 'r'))
		except IOError:
			print("Could not read file:", str(arg))
	elif opt == "--in_list":
		diff = str(arg)
	elif opt == "--index_CHROM":
		try:
			index_CHROM = int(arg)
		except ValueError:
			print("Oops!  That was no valid number for index_CHROM.  Try again...")
	elif opt == "--index_POS":
		try:
			index_POS = int(arg)
		except ValueError:
			print("Oops!  That was no valid number for index_POS.  Try again...")
	elif opt == "--index_REF":
		try:
			index_REF = int(arg)
		except ValueError:
			print("Oops!  That was no valid number for index_REF.  Try again...")
	elif opt == "--index_ALT":
		try:
			index_ALT = int(arg)
		except ValueError:
			print("Oops!  That was no valid number for index_ALT.  Try again...")
	elif opt == "--cds_info":
		cds_info = str(arg)
	elif opt == "--nt_fasta":
		ntfasta_dict = SeqIO.index(str(arg), 'fasta')
	elif opt == "--transl_table":
		try:
			transl_table = int(arg)
		except ValueError:
			print("Oops!  That was no valid number for transl_table.  Try again...")
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
			try:
				out_dict[words[key_index]]
			except:
				out_dict[words[key_index]] = {} # set each locus as a dict
			# key is feature range, value is feature info
			# dict[locus][range(feature_start, feature_end)] = [feature_type, locus_tag, gene_name, product, strand]
			out_dict[words[key_index]][range(int(words[start_index]), int(words[start_index+1])+1)] = [feature_type, words[0], words[name_index], words[name_index+1], words[name_index-1]]
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
def mutate_pos(feature_range, strand, SNP_position):
	if strand == 1: # forward strand
		position = int(SNP_position) - feature_range[0] + 1
	else: # reverse strand
		position = feature_range[-1] - int(SNP_position) + 1
	return position

# get mut_coding_dna
def SNP_mutate(nt_seq, position, alt, strand): # position start from 1
	if not strand == 1:
		complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
		alt = complement[alt]
	new_seq = nt_seq[:int(position)-1] + alt + nt_seq[int(position):]
	return new_seq

# check non-synonymous mutation
def check_mutation(transl_table, nt1, nt2):
	coding_dna_1 = Seq(nt1, IUPAC.unambiguous_dna)
	coding_dna_2 = Seq(nt2, IUPAC.unambiguous_dna)
	if str(coding_dna_1.translate(table = transl_table)) == str(coding_dna_2.translate(table = transl_table)):
		return True
	else:
		return False


### read input ###
# get list of feature ranges
all_feature_dict = {}
if cds_info:
	try:
		cds_ranges = get_range_from_info(cds_info, 'CDS', 1, 3, 7)
		all_feature_dict = cds_ranges.copy()
	except IOError:
		print("Could not read file: ", cds_info)

if pseudo_info:
	try:
		pseudo_ranges = get_range_from_info(pseudo_info, 'pseudo', 1, 3, 7)
		for locus in pseudo_ranges:
			all_feature_dict[locus].update(pseudo_ranges[locus])
	except IOError:
		print("Could not read file: ", pseudo_info)

if RNA_info:
	try:
		RNA_ranges = get_range_from_info(RNA_info, 'RNA', 1, 3, 7)
		for locus in RNA_ranges:
			all_feature_dict[locus].update(RNA_ranges[locus])
	except IOError:
		print("Could not read file: ", RNA_info)



# check differences list input, there should be only one input
# read diff as vcf records
in_list_num = 0
try:
	vcf_reader
	in_list_num += 1
except:
	pass
try:
	diff
	in_list_num += 2
	try:
		with open(diff, 'r') as in_file_h:
			vcf_reader = []
			for line in in_file_h:
				words = line.strip('\n').split('\t')
				locus = str(words[index_CHROM])
				vcf_reader.append(VcfOut(locus, words[index_POS], words[index_REF], [words[index_ALT]]))
	except IOError:
		print("Could not read file: ", diff)
except:
	pass
if in_list_num not in [1,2]:
	print("!!ERROR: One differences list required, --in_vcf or --in_list")
	exit()


# store characterized vcf record in list
out_vcf_record = []

# vcf characterization
for vcf_record_i in vcf_reader: # for each vcf record
	locus = str(vcf_record_i.CHROM).split('.')[0]
	vcf_record_o = VcfOut(locus, vcf_record_i.POS, vcf_record_i.REF, vcf_record_i.ALT)
	if '-' in vcf_record_i.REF:
		vcf_record_i.REF = ''
	elif '-' in vcf_record_i.ALT:
		vcf_record_i.ALT = ['']
	if len(vcf_record_i.REF) == len(vcf_record_i.ALT[0]):
		vcf_record_o.add_vcf_type('SNP')
	elif len(vcf_record_i.REF) > len(vcf_record_i.ALT[0]):
		vcf_record_o.add_vcf_type('deletion')
	else:
		vcf_record_o.add_vcf_type('insertion')
	intragenic_switch = False # init switch to False
	for feature_range in all_feature_dict[locus]: # for each range in the vcf locating locus 
		feature_record = all_feature_dict[locus][feature_range] # dict[locus][range(feature_start, feature_end)] = [feature_type, locus_tag, gene_name, product, strand]
		if (vcf_record_i.POS in feature_range) and (not intragenic_switch):
			intragenic_switch = True # switch on
			vcf_record_o.add_feature_info(feature_record[0], feature_record[1], feature_record[2], feature_record[3])
			if feature_record[0] == 'CDS' and vcf_record_o.vcf_type == 'SNP':
				ori_seq = str(ntfasta_dict[feature_record[1]].seq)
				strand = int(feature_record[4])
				position = mutate_pos(feature_range, strand, vcf_record_i.POS)
				for alt in vcf_record_i.ALT:
					alt = str(alt)
					new_seq = SNP_mutate(ori_seq, position, alt, strand)
					vcf_record_o.add_mut_type(check_mutation(transl_table, ori_seq, new_seq))
		elif (vcf_record_i.POS in feature_range) and (intragenic_switch): # vcf in more than one feature
			print('WARNING: ' + str(vcf_record_o.CHROM) + ': ' + str(vcf_record_o.POS) + ' locates in more than one feature!!' )
	if not intragenic_switch: # this vcf is not in any gene => it's intergenic
		vcf_record_o.add_feature_info('intergenic')
	out_vcf_record.append(vcf_record_o)


# write output
count_vcf, count_cds, count_pseudo, count_RNA, count_intergenic, count_nonsynonymous = 0, 0, 0, 0, 0, 0

out_file_h = open(out_file, 'w')
for i in out_vcf_record:
	count_vcf += 1
	vcf_record_list = [i.CHROM, i.POS, i.REF, ','.join([str(j) for j in i.ALT]), i.vcf_type, i.feature_type, ','.join([str(j) for j in i.mut_type]), i.locus_tag, i.name, i.product]
	out_file_h.write('\t'.join([str(j) for j in vcf_record_list]) + '\n')
	if i.feature_type == 'CDS':
		count_cds += 1
		if i.vcf_type == 'SNP' and 'non-synonymous' in i.mut_type:
			count_nonsynonymous += 1
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
print('vcf = %i\ncds = %i\nnon-synonymous = %i\npseudo = %i\nRNAs = %i\nintergenic = %i' % (count_vcf, count_cds, count_nonsynonymous, count_pseudo, count_RNA, count_intergenic))
# end of script
