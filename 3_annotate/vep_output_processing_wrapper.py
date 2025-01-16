#!/usr/bin/env python
# -*- coding: ASCII -*-

"""
:Author: Seyan Hu
:Date: 16-11-2022
:Usage: python <script.py> -v <VEP output path> -r <Reference chromosome path> -s <Original vcf file before VEP annotation path> -e <exome file>

Wrapper for processing the VEP output on the derived and simulated variants. 

:Example:
python vep_output_processing_wrapper.py -s output/dir_bgzip_tabix/ -r ../generate_ancestral_seq/data/genome/ -v output/dir_vep_output/ -g data/grantham_matrix/grantham_matrix_formatted_correct.tsv

:Modified: Julia Beets
:Date: 11-08-2023 

"""

# Import dependencies.
import sys,os
from optparse import OptionParser


# OptionParser for input. 
parser = OptionParser()
parser.add_option("-v","--vep", dest="vep", help="VEP (variant effect predictor) output path",default="")
parser.add_option("-r","--ref", dest="ref", help="Path to reference chr",default='')
parser.add_option("-s","--vcf_source", dest="vcf", help="Path to bgzipped & tabix index vcf source file (vcf file of the generated variants before VEP annotation)",default="") 

parser.add_option("-g","--grantham", dest="grantham", help="Path to Grantham score annotation file",default='')

(options, args) = parser.parse_args()


# Sorted list of files in directory
VEP_list = []
for file_n in os.listdir(options.vep):
	if file_n.endswith('vcf'):
		VEP_list.append(file_n)
VEP_list = sorted(VEP_list)


# Iterate over sorted list.
print('Start processing VEP!')
for fn in VEP_list:
	print('Working on: ' + fn)
	
	# Check if input is a derived variant. 
	if 'derived' in fn:
		
		# Determine chr number. 
		chr_num_d = fn.split('chr')[-1]
		chr_num_d = chr_num_d.split('_')[0]
		
		# Perform commandline. Hardcoded script name!
		os.system('python VEP-processing.py -v '+options.vep+fn+' -r '+options.ref+'ref'+chr_num_d+'.fa'+' -s '+options.vcf+'derived_var_chr_'+chr_num_d+'_case_upper.vcf.gz'+' -g '+options.grantham)
		continue
	# Check if input is a simulated variant. 
	if 'simulated' in fn:
		print("Sim file", flush=True)
		# Determine chr number. 
		chr_num_s = fn.split('chr')[-1]
		chr_num_s = chr_num_s.split("_")[0]
		
		# Perform commandline.
		os.system('python VEP-processing.py -v '+options.vep+fn+' -r '+options.ref+'ref'+chr_num_s+'.fa'+' -s '+options.vcf+'simulated_chr'+chr_num_s+'.vcf.gz'+' -g '+options.grantham)