#!/usr/bin/env python
# -*- coding: ASCII -*-

"""
:Author: Seyan Hu
:Date: 5-12-2022
:Usage: python <script.py> 

These fasta files containing repeats should be downloaded from the UCSC Genome Browser database for the species of interest. 
And should be put in the 'data/repeats/' directory. 
UCSC should have masked fasta files (repeats are in lower case) per chromosome for the species of interest. 
These files should be decompressed. 
The script creates per chromosome a output file containing a list of the position of its repeats. 


:Example:
python get_position_repeats.py -r ./

"""

# Import dependencies.
import sys, os
from os import listdir
from optparse import OptionParser
import json


# OptionParser for input. 
parser = OptionParser()
parser.add_option("-r","--rep", dest="rep", help="Path to repeats",default="../../../data/repeats/") # ../output/dir_toWig/

(options, args) = parser.parse_args()


# Define input files.
fn_list = []
for fn in os.listdir(options.rep):
	if fn.endswith('.fa'):
		fn_list.append(fn)
print(fn_list)


# Iterate over the file list.
for fn in fn_list:
	print('Working on: ' + fn)

	
	# Open file.
	open_fn = open(options.rep + fn, 'r')
	
	# Iterate over each line and convert the fasta sequence to one string. 
	print('Creating one sequence string for: ' + fn)
	seq_string = ''
	for line in open_fn:
		if not line.startswith('>'):
			seq_line = line.replace('\n', '')
			seq_string += seq_line
		else:
			chr_num = line.split(" ")[0].strip('\n').split(">chr")[1]

	# Iterate over position in sequence and add positions of repeats to list.
	print('Creating list of repeat positions for: ' + fn)
	repeat_pos = [] 
	for pos, nt in enumerate(seq_string):
		pos += 1
		if nt.islower():
			print(pos, nt)
			repeat_pos.append(pos)
	
	# Create output with json. 
	with open("list_repeats_chr" + str(chr_num) + '.txt', "w") as fp:
		json.dump(repeat_pos, fp)