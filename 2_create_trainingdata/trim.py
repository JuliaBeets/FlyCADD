#!/usr/bin/env python
'''
:Author: Julia Beets
:Date: 09-08-2023
:Usage: python <script.py> -s <path to simulated variants> -d <path to derived variants>

Trims the vcf file of simulated variants to the same amount of derived variants, whilst keeping the same proportions of variants per chromosome. 

'''

import sys, os, random
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-s", "--simu", dest="simu", help="Simulated variants", default=None)
parser.add_option("-d", "--derived", dest="derived", help="Derived variants", default=None)

(options, args) = parser.parse_args()

simu_file_list = []
derived_file_list = []
for fn in os.listdir(options.simu):
    if fn.startswith('snps_sim_variants_') and fn.endswith('_filtered.vcf'):
        simu_file_list.append(fn)
for fn in os.listdir(options.derived):
    if fn.startswith('derived_var_') and fn.endswith('_upper.vcf'):
        derived_file_list.append(fn)

simu_file_list = sorted(simu_file_list)
derived_file_list = sorted(derived_file_list)

# Count the total number of derived variants across all files
total_derived_variants = 0
for derived_fn in derived_file_list:
    with open(options.derived + derived_fn, 'r') as f:
        for line in f:
            if not line.startswith("#"):
                total_derived_variants += 1
print("Total derived variants: " + str(total_derived_variants))

# Shuffle the simulated variants and take as many as the total derived variants
for simu_fn in simu_file_list:
	output = open('trimmed_' + simu_fn, 'w')
	with open(options.simu + simu_fn, 'r') as f:
		all_lines = f.readlines()
		headers = [line for line in all_lines if line.startswith("#")]
		lines = [line for line in all_lines if not line.startswith("#")]
		if len(lines) < total_derived_variants:	
			print(f"Error: The simulated file {simu_fn} has fewer variants ({len(lines)}) than the total derived variants ({total_derived_variants}). Exiting.")
			sys.exit(1)

		random.shuffle(lines)

		for header in headers:
			print("header: "+ header)
			output.write(header)

        # Take the same number of lines as total_derived_variants
		for line in lines[:total_derived_variants]:
			output.write(line)

		output.close()