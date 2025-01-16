#!/usr/bin/env python
# -*- coding: ASCII -*-

"""
:Author: Christian Gross
:Contact: c.gross@tudelft.nl
:Date: 06-08.2018

This files goes through the ancestral sequence 
and compares it with the genome 
and the known population variants to identfy derived variants.

I make use only of high quality ancestor and assembly
Cases		:1		2		3		4		5
Ancestor	:C		C		C		C		C
Reference	:C		A>0.9	A		A>0.9	A
Alternative	:A>0.9	G		G>0.9	C		A


:Edited by: Seyan Hu
:Date: 19-10-2022
:Usage: python derived_var_gen.py <chr num> <path to ancestor seq> <path to genome> <path to frequency files>

:Example:
python derived_var_gen.py 19 ./ ../../../generate_ancestral_seq/data/genome/ ../../output/dir_freq_f/

Note:
Since this script only performs the generation of derived variants for one chromosome, 
a wrapper is needed to iterate over all chromosomes. 

:Modification: Julia Beets, 03-08-2023
"""


# Import dependencies
import os,sys
import io
from optparse import OptionParser
import time
import pysam
from pysam import VariantFile


# Assigns variables to input
chromosome_num = sys.argv[1]
ancestor_seq_path = sys.argv[2]
genome_path = sys.argv[3]
freq_path = sys.argv[4]


# Check if the nt are written in upper or lower case, 
# depending on this, the script writes it in one of the two output files.
def write_line(chrom, pos, reference, an, output_low, output_high):
	if (reference in 'ACGTacgt') and (an in 'ACGTacgt'):
		if an.islower():
			output_low.write("%s\t%s\t.\t%s\t%s\t.\t.\t.\n" %(chrom, pos + 1, reference, an))
		else:
			output_high.write("%s\t%s\t.\t%s\t%s\t.\t.\t.\n" %(chrom, pos + 1, reference, an))


# OptionParser for the chr number and start position. 
parser = OptionParser()
parser.add_option("-c", "--chromosome", dest= "chrom", help= "investigated Chromosome", default= chromosome_num)
parser.add_option("-s", "--start", dest= "start", help= "start position of the region", default= "0")
(options, args) = parser.parse_args()


# Checks whenever a input was given for OptionParser.
if len(options.chrom) == 0:
	sys.exit('No Chromosome specified, Program closed prematurely.')
if len(options.start) == 0:
	sys.exit('No Startposition specified, Program closed prematurely.')


# Define input files for the Ancestral sequence and the reference. #HARDCODED!!
ancestor_fasta = pysam.Fastafile(ancestor_seq_path + "Ancestor_117_chr%s.fa" %options.chrom)
print(ancestor_seq_path + "Ancestor_117_chr%s.fa" %options.chrom)
ref_fasta = pysam.Fastafile(genome_path + "ref%s.fa"%options.chrom)
print(genome_path + "ref%s.fa"%options.chrom)

# Open population frequency files
# Note: frequency files are not vcf files
input_snps = open(freq_path + "freq_%s.out"%options.chrom, 'r')
print(input_snps)

#HERE CREATE FIX FOR FINDING CHROMOSOME NUMBER, NOW NOT WORKING!!!
# Check if the ancestor seq is of same length as the reference
if ancestor_fasta.nreferences != 1:
	sys.exit('There are more than 1 record in the fasta file.')
else:
	ancestor_record = options.chrom
	print(ancestor_record)
	
	if '%s'%options.chrom in ancestor_record:
		print(ancestor_fasta.references[0])
		anc_length = ancestor_fasta.get_reference_length(ancestor_fasta.references[0])
		print("test "+ str(anc_length))
	elif 'chr%s'%options.chrom in ancestor_record:
		anc_length = ancestor_fasta.get_reference_length(ancestor_fasta.references[0])
		print(anc_length)
	else:
		sys.exit('The requested chromosome cannot be found in the record of the ancestor fasta\n%s\n%s'%('.%s:'%options.chrom, ancestor_record))
if anc_length != ref_fasta.get_reference_length(ref_fasta.references[0]):
	sys.exit('Ancestor fasta and ref fasta have not the same length, Chromosome: %s'%('%s'%options.chrom))
print('Ancestor sequence control: Chr' + str(chromosome_num) +', sequence is good')



# Create output files for upper and lower case variants.
output_high = open("derived_var_chr_%s_case_upper.vcf" % options.chrom, 'w')
output_low = open("derived_var_chr_%s_case_lower.vcf" % options.chrom, 'w')
output_low.write("##fileformat=VCFv4.0\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
output_high.write("##fileformat=VCFv4.0\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")


# Split the ancestor and reference sequences in chunks of 350000.
# Done so that the script reads this in chunks. 
index_outer = (int(options.start), int(options.start) + 350000)
snp_line = input_snps.readline()
snp_line = input_snps.readline()

snp_pos = int(snp_line.strip().split('\t')[1])

while index_outer[0] <= anc_length:
	print(index_outer[0], index_outer[1])
	# Fetch the reference and ancestral sequence.
	ref_seq = ref_fasta.fetch(str(ref_fasta.references[0]),index_outer[0], index_outer[1])
	anc_seq = ancestor_fasta.fetch(ancestor_record, index_outer[0], index_outer[1])
	# Iterate over the ancestral sequence position.
	for i, anc in enumerate(anc_seq):
		
		# Enumerate sets i always to 0, therefore i needs to be summed with the starting position of index_outer
		# Checking if current position in the sequences is equal to the snp position in the frequency files
		if i + index_outer[0] == (snp_pos-1):
			snp_allele = snp_line.strip().split('\t')[5].split(':')[0]
			
			# Checking if Case 1 derived snp is given
			if (anc == ref_seq[i]) and (anc != snp_allele) and (anc != '-') and (anc != '.'):
				# Case 1 write out
				print("case1")
				write_line(options.chrom, i + index_outer[0], ref_seq[i], snp_allele, output_low, output_high)
				# 'Try' is used, to avoid error when the end of snp file is reached but not yet the end of the sequence. 
				try: 
					snp_line = input_snps.readline()
					snp_pos = int(snp_line.strip().split('\t')[1])
				except:
					pass
			
			# Checking if case 3 (or 2) derived snp is given
			elif (anc != ref_seq[i]) and (anc != snp_allele) and (anc != '-'):
				print("case2,3")
				# Case 3 write out
				write_line(options.chrom, i + index_outer[0], ref_seq[i], snp_allele, output_low, output_high)
				try: 
					snp_line = input_snps.readline()
					snp_pos = int(snp_line.strip().split('\t')[1])
				except: 
					pass
			
			# Only take next line in snp file
			else: 
				print("here") 
				try:
					snp_line = input_snps.readline()
					snp_pos = int(snp_line.strip().split('\t')[1])
				except: 
					pass
					
		elif anc != '-' :
			print("anc is ref if not followed by ref is not")
			# Checking if case 5 (or 3) derived variant
			if (anc != ref_seq[i]):
				print("ref is not anc:"+ref_seq[i]+anc)
				# Case 5 write out
				write_line(options.chrom, i + index_outer[0], ref_seq[i], anc, output_low, output_high)
				
	index_outer = (index_outer[1], index_outer[1] + 350000)

output_high.close()
output_low.close()


