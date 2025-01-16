#!/usr/bin/env python
# -*- coding: ASCII -*-

"""
:Author: Seyan Hu
:Date: 16-11-2022
:Usage: python <script.py> -i <path to merged annotation files>

Wrapper for handling encodings and missing values. 
It takes the merged tsv files as input and outputs csv files. 

"""

# Import dependencies.
import sys,os
from optparse import OptionParser


# OptionParser for input. 
parser = OptionParser()
parser.add_option("-i", "--input", dest="input", help="Path to merged annotation files",default="./")

(options, args) = parser.parse_args()


# Sorted list of files in directory. 
files_list = []
for fn in os.listdir(options.input):
	if '_chunk' in fn:
		files_list.append(fn)
files_list = sorted(files_list)


# Iterate over sorted list.
print('Performing encoding and imputation!')
for fn in files_list:
	print('Working on: ' + fn) 
	
	# Determine chr number.
	chr_num = fn.split('_')[0]
	
	# Perform encoding and imputation script. 
	os.system('python data_encoding_mv_handling.py -i ' + options.input + fn + ' -p /step3/output/dir_encoding_mv/')


