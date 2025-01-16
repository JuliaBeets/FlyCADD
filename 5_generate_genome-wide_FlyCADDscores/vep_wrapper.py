#!/usr/bin/env python
# -*- coding: ASCII -*-

"""
:Author: Seyan Hu
:Date: 16-11-2022
:Usage: python <script.py> -i <Path to input (all possible variants)> -s <Species of interest>

Performes VEP on all variants. 
These variants needs to be in the same directory and the script iterates over the files and performs VEP.

:Example: 
python vep_wrapper.py -i ../data/variants/ -s mus_musculus

:Modified: Julia Beets
:Date: 23-08-2023

"""

# Import dependencies.
import math
import sys,os
from optparse import OptionParser


# OptionParser for input. 
parser = OptionParser()
parser.add_option("-i", "--input", dest="input", help="Name and path of input files", default= "output/dir_all_variants/")
parser.add_option("-s", "--sp", dest="sp", help="Species of interest", default= "mus_musculus")
(options, args) = parser.parse_args()


# Sorted list of files in directory
file_list = os.listdir(options.input)
file_list = sorted(file_list)


# Iterate over sorted list.
print('Performing VEP!')
for fn in file_list:
	print('Working on: ' + fn) 
		
	# Determine chr number. 
	chr_num_d = fn.split('_')[0]
		
	# Perform commandline.
	os.system('vep --input_file '+ options.input + fn +' --cache --offline --buffer 1000 --no_stats --species '+ options.sp +' --format vcf --regulatory --per_gene --ccds --domains --numbers --canonical --total_length --force_overwrite --output_file temp.vcf --quiet')
	os.system('''cat temp.vcf | awk 'BEGIN{ FS="\t"; OFS="\t"; }{ if ($1 ~ /^#/) { if ($1 ~ /^#Up/) { sub("#","",$1); print "#Chrom","Start","End",$0 } else { print } } else { split($2,a,":"); split(a[2],b,"-"); if (length(b) == 2) { print a[1],b[1],b[2],$0 } else { print a[1],b[1],b[1],$0 } }}' >> '''+ chr_num_d +'''_VEP-annotated.vcf''')