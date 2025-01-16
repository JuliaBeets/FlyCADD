#!/usr/bin/env python
# -*- coding: ASCII -*-

"""
:Author: Seyan Hu
:Date: 16-11-2022
:Usage: python <script.py> -i <path to merged annotation files>

Wrapper for handling encodings and missing values. 
It takes the merged tsv files as input and outputs csv files. 

:Modified: Julia Beets
:Date: 15-08-2023

"""

# Import dependencies.
import sys,os
from optparse import OptionParser


# OptionParser for input. 
parser = OptionParser()
parser.add_option("-i", "--input", dest="input", help="Path to merged annotation files",default="./")

(options, args) = parser.parse_args()


# Sorted list of files in directory with only annotations of simulated variants. 
simulated_files_list = []
for fn in os.listdir(options.input):
	if 'simulated' in fn:
		simulated_files_list.append(fn)
simulated_files_list = sorted(simulated_files_list)


# Iterate over sorted list.
print('Performing encoding and imputation!')
for fn in simulated_files_list:
	print('Working on: ' + fn) 
	
	# Determine chr number.
	chr_num = fn.split('_')[0]
	
	# Perform encoding and imputation script on simulated variant files first. Hardcoded script name!
	print('python data_encoding_mv_handling.py -i ' + options.input + fn)
	os.system('python data_encoding_mv_handling.py -i ' + options.input + fn)
	
	# Perform encoding and imputation on derived variant files. Hardcoded script name!
	print('python data_encoding_mv_handling.py -i ' + options.input + chr_num + '_derived_merged_features.tsv -d')
	os.system('python data_encoding_mv_handling.py -i ' + options.input + chr_num + '_derived_merged_features.tsv -d')