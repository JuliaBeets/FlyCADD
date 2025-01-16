#!/usr/bin/env python
# -*- coding: ASCII -*-

"""
:Author: Seyan Hu
:Date: 5-12-2022
:Usage: python <script.py> -f <Path to formatted file containing phastCons/phyloP scores> -c <List of considered chromosomes>

The script splits the formatted file into subfiles per chromosome. 

:Example:
python split_pC_pP_scores.py -f ./phastCons_scores.txt -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,X

"""

# Import dependencies.
import sys, os
from os import listdir
from optparse import OptionParser


# OptionParser for input. 
parser = OptionParser()
parser.add_option("-f","--file", dest="file", help="Path to formatted file containing phastCons/phyloP scores",default="./phastCons_scores.txt")
parser.add_option("-c","--chr", dest="chr", help="List of considered chromosomes",default="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,X")

(options, args) = parser.parse_args()


# Iterate over the chromosomes.  
for chr_num in options.chr.split(','):
	
	# Create output.
	if 'phastCons' in options.file:
		outp = open(chr_num + '_phastCons.tsv', 'w')
	elif 'phyloP' in options.file:
		outp = open(chr_num + '_phyloP.tsv', 'w')
	outp.write('#chr_num\tstart_position\tscore\n')
	
	# Open input.
	open_inp = open(options.file, 'r')
	
	# Iterate over lines in opened file. 
	for line in open_inp:
		split_line = line.split('\t')
		chr_number = split_line[0]
		if chr_number == chr_num:
			outp.write(line)