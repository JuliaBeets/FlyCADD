#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Original Author: Christian Gross
Contact: cgross@tudelft.nl
Date: 18-11-2017

Modification:
Julia Beets (j.beets@vu.nl)
15-7-2023

"""

# Import dependencies
import sys, os, argparse
from optparse import OptionParser, OptionGroup
import gzip
from Bio import Phylo
from Bio import AlignIO
from io import StringIO
from sh import gunzip

# OptionParser for the input directory
parser = OptionParser()
parser.add_option("-p", "--path", dest="path", help="path to folder with alignment files", default="./")
parser.add_option("-a", "--ancestor", dest="ancestor", help="sequence which should be marked as ancestror", default="")
parser.add_option("-i", "--identifier", dest="identifier", help="identifier to get the correct files; (eg. chr1[.maf])",
                  default="")
parser.add_option("-f", "--full", dest="name",
                  help="full scientific name of species of interest, separated with an underscore; (eg. homo_sapiens)",
                  default="")

if len(sys.argv) == 1:
    parser.print_help()
    parser.exit()
(options, args) = parser.parse_args()

# CHECK input
if not os.path.exists(options.path):
    sys.stderr.write("Invalid path to given msa directory.\n")
    sys.exit()
if (not options.path.endswith('/')):
    options.path = options.path + '/'

if options.identifier == '':
    sys.stderr.write("Found no files with given identifier (prefix).\n")
    sys.exit()
if options.name == '':
    sys.stderr.write("No full species name given!\n")
    sys.exit()

# Reformat input to fit MAF format
name_sp1 = 's\t' + options.name + '.'
name_ancestor = 's\t' + options.ancestor + '.'

# creates a list of files to be analysed
file_list = os.listdir(options.path)
processed_file_list = []
for file in file_list:
    if (options.identifier in file) and (not 'README' in file):
        processed_file_list.append(file)

if not os.path.isdir('marked'):
    os.system('mkdir marked')

for file in processed_file_list:
    print('Marking ancestor in file: {}'.format(file))

    outfile = open('marked/marked_' + file.replace('.gz', ''), "w")

    if file.endswith(".gz"):
        gunzip(options.path + file)
    alignment_file = open(options.path + file.replace('.gz', ''), "r")
    for lines in alignment_file:
        if lines.startswith(name_sp1):
            # Write out given species
            outfile.write(lines)
        elif lines.startswith(name_ancestor):
            newline = lines.strip().split('.')[1]
            outfile.write('s\tAncestor_' + options.ancestor + '.' + newline + '\n')
        else:
            outfile.write(lines)
    os.system('gzip marked/marked_' + file.replace('.gz', ''))
    outfile.close()
    os.system('gzip ' + str(options.path + file.replace('.gz', '')))
    alignment_file.close()

if not os.path.isdir('output'):
    os.system('mkdir output')
