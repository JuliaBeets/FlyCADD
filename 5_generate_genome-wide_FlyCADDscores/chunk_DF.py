#!/usr/bin/env python
# -*- coding: ASCII -*-

"""
:Author: Seyan Hu
:Date: 31-1-2022
:Usage: python <script.py> -i <path to merged annotation files>

Split merged DF into chunks. 

"""

# Import dependencies.
import sys,os
from optparse import OptionParser
import pandas
import numpy as np


# OptionParser for input. 
parser = OptionParser()
parser.add_option("-i", "--input", dest="input", help="Path to merged annotation files",default="./test_runs/")
parser.add_option("-l", "--lines", dest="lines", help="The number of lines in each chunk",default="1000000")

(options, args) = parser.parse_args()


# Sorted list of files in directory. 
files_list = []
for fn in os.listdir(options.input):
	if 'merged_features.tsv' in fn:
		files_list.append(fn)
files_list = sorted(files_list)

# Iterate over sorted list.
print('Seperate DF into chunks')
for fn in files_list:
	print('Working on: ' + fn) 
	
	# Determine chr number.
	chr_num = fn.split('_')[0]

	# Perform splitting
	os.system('split -l'+str(options.lines)+' --verbose '+options.input+fn+' '+chr_num+'_chunk.')
	
	# Make list of output. 
	outp_list = []
	for chunk_f in os.listdir('./'):
		if not '_chunk.aa' in chunk_f:
			if '_chunk.' in chunk_f:
				outp_list.append(chunk_f)
	outp_list = sorted(outp_list)

	# Check first line in chunk.aa
	with open(chr_num+'_chunk.aa') as f:
		first_line = f.readline().strip('\n')

	# Insert first line: 
	for chunk_file in outp_list:
		outp_fn = chunk_file.replace('chunk', 'chunk_ln')
		os.system("sed '1 i "+ first_line +"' "   +chunk_file+' > '+outp_fn)