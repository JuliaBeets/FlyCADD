#!/usr/bin/env python
'''
:Author: Seyan Hu
:Date: 19-9-2022
:Usage: python wrapper_derived_gen.py <list of chromosomes> <path to ancestor seq> <path to genome> <path to frequency files>

:Example:
python scripts/wrapper_derived_gen.py '1','2','3','4','5' ../generate_ancestral_seq/output/dir_generated_ancestor_seq/ ../generate_ancestral_seq/data/genome/ output/dir_freq_f/

This script is a wrapper for derived_var_gen.py. 
Since that script can only take one input file, 
this script is used as a wrapper to loop through multiple input files.

It is used to run this command line, looping through multiple chromosomes: 
python scripts/derived_var_gen.py <one chr num> <path to ancestor seq> <path to genome> <path to frequency files>
'''

# Import dependencies
import sys, os
from optparse import OptionParser


# Assign input to variables
chromosome_list = sys.argv[1]
ancestor_seq_path = sys.argv[2]
genome_path = sys.argv[3]
freq_path = sys.argv[4]

# OptionParser for the directories of the input 
parser = OptionParser()
parser.add_option("-a", "--ancestor_path", dest="an_p", help="Path to ancestor folder", default= ancestor_seq_path)
parser.add_option("-g", "--genome_path", dest="g_p", help="Path to genome folder", default= genome_path)
parser.add_option("-f", "--frequency_path", dest="f_p", help="Path to frequency files folder", default= freq_path)

(options, args) = parser.parse_args()

# Creates a list for the chromosomes
chr_list = chromosome_list.split(',')


# Loops through each chromosome and performs the commandline
for chr_number in chr_list:
	os.system('python derived_var_gen.py ' + chr_number + ' ' + ancestor_seq_path + ' ' + genome_path + ' ' + freq_path)
