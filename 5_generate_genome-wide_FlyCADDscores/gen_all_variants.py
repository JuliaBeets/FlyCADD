'''
:Author: Seyan Hu
:Date: 10-1-2023
:Usage: python <script.py> -i <Path to genome files>

:Modified: Julia Beets
:Date: 23-08-2023

This script generates all possible variants for every position in the fasta files of the reference genome and outputs a CSV file for each chromosome. 
Note: Reference genome filenames should have this format: ref19.fa

'''

# Import dependencies.
import sys,os
from optparse import OptionParser


# OptionParser for input. 
parser = OptionParser()
parser.add_option("-i", "--input", dest="input", help="Path to genome files", default= "../data/")

(options, args) = parser.parse_args()


# Sorted list of files in directory.
input_list = []
for fn in os.listdir(options.input):
	if '.fa' in fn:
		input_list.append(fn)
input_list = sorted(input_list)


# Iterate over sorted list.
print('Iterating over input files')
for fn in input_list:
	print('Working on: ' + fn)
	
	# Determine chr number.
	chr_num = fn.split('.')[0].replace('ref', '')
	
	# Create output file.
	output = open(str(chr_num) + "_all_variants.vcf", 'w')
	output.write("##fileformat=VCFv4.0\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
	
	# Open chromosome file.
	open_f = open(options.input + fn, 'r')
	
	# Create one single continues sequence.
	seq = ''
	for line in open_f:
		if not line.startswith('>'):
			seq_line = line.replace('\n', '')
			seq += seq_line

	# Iterate over the nucleotides in 'seq'.
	pos = 1
	for nt in seq:
		# Write all possible variants to output.
		if nt in 'Aa':
			output.write("%s\t%s\t.\t%s\t%s\t.\t.\t.\n"%(chr_num, pos, 'A', 'T'))
			output.write("%s\t%s\t.\t%s\t%s\t.\t.\t.\n"%(chr_num, pos, 'A', 'C'))
			output.write("%s\t%s\t.\t%s\t%s\t.\t.\t.\n"%(chr_num, pos, 'A', 'G'))
		elif nt in 'tT':
			output.write("%s\t%s\t.\t%s\t%s\t.\t.\t.\n"%(chr_num, pos, 'T', 'A'))
			output.write("%s\t%s\t.\t%s\t%s\t.\t.\t.\n"%(chr_num, pos, 'T', 'C'))
			output.write("%s\t%s\t.\t%s\t%s\t.\t.\t.\n"%(chr_num, pos, 'T', 'G'))
		elif nt in 'Cc':
			output.write("%s\t%s\t.\t%s\t%s\t.\t.\t.\n"%(chr_num, pos, 'C', 'T'))
			output.write("%s\t%s\t.\t%s\t%s\t.\t.\t.\n"%(chr_num, pos, 'C', 'A'))
			output.write("%s\t%s\t.\t%s\t%s\t.\t.\t.\n"%(chr_num, pos, 'C', 'G'))
		elif nt in 'Gg':
			output.write("%s\t%s\t.\t%s\t%s\t.\t.\t.\n"%(chr_num, pos, 'G', 'T'))
			output.write("%s\t%s\t.\t%s\t%s\t.\t.\t.\n"%(chr_num, pos, 'G', 'C'))
			output.write("%s\t%s\t.\t%s\t%s\t.\t.\t.\n"%(chr_num, pos, 'G', 'A'))
		pos += 1