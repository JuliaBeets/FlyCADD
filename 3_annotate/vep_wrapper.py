#!/usr/bin/env python
# -*- coding: ASCII -*-

"""
:Author: Seyan Hu
:Date: 16-11-2022
:Usage: python <script.py> -i <Path to input (both derived and simulated variants)> -s <Species of interest>

Performes VEP on the derived and simulated variants. 
These variants needs to be in the same directory, split by chromosome and the script iterates over the files and performs VEP.

:Example: 
python vep_wrapper.py -i ../data/variants/ -s mus_musculus

:Modified: Julia Beets
:Date: 14-08-2023
"""

# Import dependencies.
import math
import sys,os
from optparse import OptionParser


# OptionParser for input. 
parser = OptionParser()
parser.add_option("-i", "--input", dest="input", help="Name and path of input files", default= "../data/variants/")
parser.add_option("-s", "--sp", dest="sp", help="Species of interest", default= "mus_musculus")
(options, args) = parser.parse_args()


# Sorted list of files in directory
file_list = os.listdir(options.input)
file_list = sorted(file_list)


# Iterate over sorted list.
print('Performing VEP!')
for fn in file_list:
	print('Working on: ' + fn) 
	
	# Check if input is a derived variant. 
	if 'derived' in fn:
		
		# Determine chr number. 
		chr_num_d = fn.split('chr_')[-1]
		chr_num_d = chr_num_d.replace('_case_upper.vcf', '')
		
		# Perform commandline.
		os.system('vep --input_file '+ options.input + fn +' --quiet --cache --offline --buffer 1000 --no_stats --species '+ options.sp +' --format vcf --regulatory --per_gene --ccds --domains --numbers --canonical --total_length --force_overwrite --output_file temp.vcf')
		os.system('''cat temp.vcf | awk 'BEGIN{ FS="\t"; OFS="\t"; }{ if ($1 ~ /^#/) { if ($1 ~ /^#Up/) { sub("#","",$1); print "#Chrom","Start","End",$0 } else { print } } else { split($2,a,":"); split(a[2],b,"-"); if (length(b) == 2) { print a[1],b[1],b[2],$0 } else { print a[1],b[1],b[1],$0 } }}' >> anc117/annotated/derived_chr'''+ chr_num_d +'''_VEP-annotated.vcf''')
		
	# Check if input is a simulated variant.
	elif 'simulated' in fn:
		
		# Determine chr number. 
		chr_num_s = fn.split('chr')[1]
		chr_num_s = chr_num_s.replace('.vcf', '')

		# Perform commandline.
		os.system('vep --input_file '+ options.input + fn +' --quiet --cache --offline --buffer 1000 --no_stats --species '+ options.sp +' --format vcf --regulatory --per_gene --ccds --domains --numbers --canonical --total_length --force_overwrite --output_file temp.vcf')
		os.system('''cat temp.vcf | awk 'BEGIN{ FS="\t"; OFS="\t"; }{ if ($1 ~ /^#/) { if ($1 ~ /^#Up/) { sub("#","",$1); print "#Chrom","Start","End",$0 } else { print } } else { split($2,a,":"); split(a[2],b,"-"); if (length(b) == 2) { print a[1],b[1],b[2],$0 } else { print a[1],b[1],b[1],$0 } }}' >> anc117/annotated/simulated_chr'''+ chr_num_s +'''_VEP-annotated.vcf''')