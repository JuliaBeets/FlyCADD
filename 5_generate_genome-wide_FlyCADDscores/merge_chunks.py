'''
:Author: Seyan Hu
:Date: 6-2-2023
:Usage: python <script.py> -i <Path to processed chunks> 

The script checks for the available chromosomes in the given path. 
Iterates over the chromosome numbers and merger per chromosome their chunks. 

:Modified: Julia Beets
:Date: 9-5-2024

'''

# Import dependencies.
import sys,os
from optparse import OptionParser


# OptionParser for input. 
parser = OptionParser()
parser.add_option("-i", "--input", dest="input", help="Path to processed chunks", default= "./")

(options, args) = parser.parse_args()


# Sorted list of chromosomes in directory.
input_chr_list = []
for fn in os.listdir(options.input):
	if 'processed' in fn:
		fn = fn.split('_')[0]
		input_chr_list.append(fn)
		
input_chr_list = list(set(input_chr_list))
input_chr_list = sorted(input_chr_list)


# Iterate over sorted list.
print('Iterating over input files')
for chr_num in input_chr_list:
	
	# Merge chunks. 
	# Performs this: cat * > merged-file
	os.system('cat ' + options.input+chr_num + '* > ' + chr_num + '_merged_chunks.csv')