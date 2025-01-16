#!/usr/bin/env python

"""
:Author: Seyan Hu
:Date: 29-9-2022
:Extension and modification: Julia Höglund
:Date: 20-3-2023
:Usage: python maftools_sorter.py -p <maf files from chr sorting path> -s <species> -o <new path to new files>

Sorts alignment blocks with respect to coordinates of species of interest.

"""

# Import dependencies
import sys, os
from optparse import OptionParser
import gzip
from sh import gunzip

parser = OptionParser()

parser.add_option("-p", "--path", dest="path", help="Path to input folder (i.e. output from previous step)", default= './sorted')
parser.add_option("-s", "--species", dest="species", help="Species of interest", default= '')
parser.add_option("-o", "--ordered", dest="sorted", help="directory within which the new output will be saved", default="sortedMSA")
parser.add_option("-c", "--clean", dest="clean", help="remove intermediate files after computation [yes/no]", default="no")

if len(sys.argv)==1:
   	parser.print_help()
   	parser.exit()
(options, args) = parser.parse_args()

# Checking if the path ends with '/'
if (not options.path.endswith('/')):
	options.path = options.path+'/'

## if path not exist makedir sorted else not do it
if not os.path.isdir(options.sorted):
    os.system('mkdir ' + str(options.sorted))

# Create list for input
file_list = []
for fn in os.listdir(options.path):
	if fn.startswith(options.species): 
		file_list.append(fn)

file_list = sorted(file_list)

# Loop through input list
for inp in file_list:
	inp_f = options.path + inp
	print('Working on: ' + inp)

	# Check chr number.
	chrom = inp.split('.')[0]
	chr_num = chrom.replace('chr', '')

	# Determine sequence name.
	seq_ident = options.species +'.'+ chr_num

	# Perform mafSorter
	os.system('mafSorter --maf ' + inp_f + ' --seq ' + seq_ident + ' > ' + options.sorted + '/mS_' + str(chrom) + '.maf')

if options.clean=='yes':
	os.system('rm -r ' + str(options.path))