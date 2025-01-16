#!/usr/bin/env python
# -*- coding: ASCII -*-

"""
:Author: Seyan Hu
:Date: 1-12-2022


Adds the y (label) to the processed VEP files and adds for every variant the correct phastCons, phyloP scores and repeat position and DNAshapeR annotation for Roll, MGW, . 
This script can be expanded for more annotations (features). 

:Example:
python merge_all_annotations.py -v ./ -g ./ -p ./ -r ./

:Modified: Julia Beets
:Date: 15-08-2023
:Usage: python <script.py> -v <Path to the processed VEP output> -p <Path to phastCons/phyloP output> -r <Path to repeats output> -g <Path to in-gene positions>
"""

# Import dependencies.
import sys, os
from os import listdir
from os.path import isfile, join
from optparse import OptionParser
import pandas as pd
import json


# OptionParser for input. 
parser = OptionParser()
parser.add_option("-v","--vep", dest="vep", help="Path to processed VEP output",default="./")
parser.add_option("-p","--pCpP", dest="pCpP", help="Path to phastCons/phyloP output",default="./")
parser.add_option("-r","--re", dest="re", help="Path to repeats output",default="./")
parser.add_option("-g","--gene", dest="gene", help="Path to in-gene output",default="./")
parser.add_option("-d","--dnashaper", dest="dnashaper", help="Path to dnashaper output",default="./")
(options, args) = parser.parse_args()


# Creates list of VEP files. 
VEP_list = []
for fn in os.listdir(options.vep):
	if fn.endswith('.vcf'): # and fn.startswith('derived'):
		VEP_list.append(fn)
VEP_list = sorted(VEP_list)


# Create list of phastCons/phyloP, repeats and in-gene files. 
pC_list = []
pP_list = []
for fn in os.listdir(options.pCpP):
	if fn.endswith('phastCons.tsv'):
		pC_list.append(fn)
	elif fn.endswith('phyloP.tsv'):
		pP_list.append(fn)

roll_list = []
EP_list = []
HelT_list = []
MGW_list = []
ProT_list = []
for fn in os.listdir(options.dnashaper):
	if fn.startswith('Roll'):
		roll_list.append(fn)
	elif fn.startswith('EP'):
		EP_list.append(fn)
	elif fn.startswith('HelT'):
		HelT_list.append(fn)
	elif fn.startswith('MGW'):
		MGW_list.append(fn)
	elif fn.startswith('ProT'):
		ProT_list.append(fn)
		
repeats_list = []
for fn in os.listdir(options.re):
	if fn.startswith('list_repeats'):
		repeats_list.append(fn)

gene_list = []
for fn in os.listdir(options.gene):
	print("hello" + fn)
	if fn.startswith('positions'):
		print("yes")
		gene_list.append(fn)


# Iterate over derived VEP list (over the chromosomes). 
for filen in VEP_list:
	
	# Open in pandas and add labels for y (derived = 1, simulated = 0).
	df_vep = pd.read_csv(options.vep + filen, sep='\t', lineterminator='\n')
	
	# Determine chr number.
	chr_num = filen.split('_')[0]
	print('Working on: Chr' + chr_num)	
	
	for fn in pC_list:
		num = fn.replace('_phastCons.tsv', '')
		if chr_num == num:
			pC_f = fn
			break
	for fn in pP_list:
		num = fn.replace('_phyloP.tsv', '')
		if chr_num == num:
			pP_f = fn
			break

	for fn in roll_list:
		num = fn.split("_")[-1].replace(".txt", "")
		if chr_num == num:
			roll_f = fn
			break

	for fn in EP_list:
		num = fn.split("_")[-1].replace(".txt", "")
		if chr_num == num:
			EP_f = fn
			break

	for fn in HelT_list:
		num = fn.split("_")[-1].replace(".txt", "")
		if chr_num == num:
			HelT_f = fn
			break

	for fn in MGW_list:
		num = fn.split("_")[-1].replace(".txt", "")
		if chr_num == num:
			MGW_f = fn
			break

	for fn in ProT_list:
		num = fn.split("_")[-1].replace(".txt", "")
		if chr_num == num:
			ProT_f = fn
			break
	 
	for fn in repeats_list:
		num = fn.split('chr')[1].replace('.txt', '')
		if chr_num == num:
			repeats_f = fn
			break
 
	for fn in gene_list:
		num = fn.split('positions_')[1].replace('.txt', '')
		if chr_num == num:
			gene_f = fn
			break



	# Open files. 
	print('Opening files as pandas DF. ')
	df_pC = pd.read_csv(options.pCpP + pC_f, sep='\t', lineterminator='\n')
	df_pP = pd.read_csv(options.pCpP + pP_f, sep='\t', lineterminator='\n')

	df_roll = pd.read_csv(options.dnashaper + roll_f, sep='\t', lineterminator='\n')
	df_EP = pd.read_csv(options.dnashaper + EP_f, sep='\t', lineterminator='\n')
	df_HelT = pd.read_csv(options.dnashaper + HelT_f, sep='\t', lineterminator='\n')	
	df_MGW = pd.read_csv(options.dnashaper + MGW_f, sep='\t', lineterminator='\n')
	df_ProT = pd.read_csv(options.dnashaper + ProT_f, sep='\t', lineterminator='\n')

	open_repeats = open(options.re + repeats_f, 'r')
	repeats_json = json.load(open_repeats)
	
	open_gene = open(options.gene + gene_f, 'r')
	gene_json = json.load(open_gene)
	

	# Duplicate positions needs to be removed from the DFs.
	print('Remove duplicate positions from DF.')

	df_pC.drop_duplicates(subset='start_position', keep="last", inplace=True)
	df_pP.drop_duplicates(subset='start_position', keep="last", inplace=True)
	
	
	## Merge annotations.
	# Merge DataFrames based on the positions of the variants.
	print('Merging DataFrames')


	# Merge PhastCons scores and retain rows of the left DF while merging.
	df_pC = df_pC.drop("#chr_num", axis='columns')
	df_pC.rename(columns = {'start_position':'Pos'}, inplace = True)
	df_merged = pd.merge(df_vep, df_pC, on='Pos', how='left')
	df_merged.rename(columns = {'score':'PhastCons'}, inplace = True)
	
	# Merge PhyloP scores and retain rows of the left DF while merging.
	df_pP = df_pP.drop("#chr_num", axis='columns')
	df_pP.rename(columns = {'start_position':'Pos'}, inplace = True)
	df_merged = pd.merge(df_merged, df_pP, on='Pos', how='left')
	df_merged.rename(columns = {'score':'PhyloP'}, inplace = True)
	print("pcpp merged")

	df_roll = df_roll.drop("Chromosome", axis='columns')
	df_roll.rename(columns={'Position': 'Pos', }, inplace=True)
	df_merged = pd.merge(df_merged, df_roll, on='Pos', how='left')
	print("roll merged")

	df_EP = df_EP.drop("Chromosome", axis='columns')
	df_EP.rename(columns={'Position': 'Pos'}, inplace=True)
	df_merged = pd.merge(df_merged, df_EP, on='Pos', how='left')
	print("EP merged")

	df_MGW = df_MGW.drop("Chromosome", axis='columns')
	df_MGW.rename(columns={'Position': 'Pos'}, inplace=True)
	df_merged = pd.merge(df_merged, df_MGW, on='Pos', how='left')
	print("mgw merged")

	df_HelT = df_HelT.drop("Chromosome", axis='columns')
	df_HelT.rename(columns={'Position': 'Pos'}, inplace=True)
	df_merged = pd.merge(df_merged, df_HelT, on='Pos', how='left')
	print("helt merged")

	df_ProT = df_ProT.drop("Chromosome", axis='columns')
	df_ProT.rename(columns={'Position': 'Pos'}, inplace=True)
	df_merged = pd.merge(df_merged, df_ProT, on='Pos', how='left')		
	print("prot merged")

	# Add Repeats data and replace all NaN values with '-'.
	df_merged['Repeats'] = df_merged['Pos'].isin(repeats_json).astype(int)
	df_merged = df_merged.fillna('-')

	# Add in-gene data and replace all NaN values with '-'.
	df_merged['ingene'] = df_merged['Pos'].isin(gene_json).astype(int)
	df_merged = df_merged.fillna('-')
	
	
	# Write new dataframe to output.
	df_merged.to_csv(str(chr_num) + '_merged_features.tsv', index=None, sep='\t')