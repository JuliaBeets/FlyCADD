#!/usr/bin/env python
# -*- coding: ASCII -*-

"""
:Author: Seyan Hu
:Date: 16-11-2022
:Usage: python <script.py> -v <VEP output path> -r <Reference chromosome path> -g <grantham file>

Wrapper for processing the VEP output on the derived and simulated variants. 

:Example:
python ver_output_processing_wrapper.py -r ../generate_ancestral_seq/data/genome/ -v output/dir_vep_output/ -g data/grantham_matrix/grantham_matrix_formatted_correct.tsv

"""

# Import dependencies.
import sys,os
from optparse import OptionParser


# OptionParser for input. 
parser = OptionParser()
parser.add_option("-v","--vep", dest="vep", help="VEP (variant effect predictor) output path",default="")
parser.add_option("-r","--ref", dest="ref", help="Path to reference chr",default='')

parser.add_option("-g","--grantham", dest="grantham", help="Path to Grantham score annotation file",default='')

(options, args) = parser.parse_args()

print("wrapper started", flush=True) 

# Sorted list of files in directory
VEP_list = []
for file_n in os.listdir(options.vep):
	if file_n.endswith('vcf'):
		VEP_list.append(file_n)
VEP_list = sorted(VEP_list)


# Iterate over sorted list.
print('Start processing VEP!')
for fn in VEP_list:
	print('Working on: ' + fn, flush=True)
		
	# Determine chr number. 
	chr_num_d = fn.replace('_VEP-annotated.vcf', '')

	# Perform commandline.
	os.system('python VEP-processing.py -v '+options.vep+fn+' -r '+options.ref+'ref'+chr_num_d+'.fa'+' -g '+options.grantham)