#!/usr/bin/env python
# -*- coding: ASCII -*-

"""
:Author: Julia Beets
:Date: 26-10-2023


Adds additional annotations to merged_annotations: repeattype, TFBS, miRNA, ReMap, CRM, 124Conservation, BG3 and S2 chromatin states. 
This script can be expanded for more annotations (features). 

:Example:
python add_annotations.py -o <Path to merged annotations> -a <Path to annotation to add> -b <Path to second annotation to add> -c <Path to third annotation to add> -d <Path to fourth annotation to add> etc

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
parser.add_option("-o","--annotated", dest="annotated", help="Path to annotated variants",default="./")

parser.add_option("-a","--a", dest="a", help="Path to annotation a",default="./")
parser.add_option("-b","--b", dest="b", help="Path to annotation b",default="./")
parser.add_option("-c","--c", dest="c", help="Path to annotation c",default="./")
parser.add_option("-d","--d", dest="d", help="Path to annotation d",default="./")
parser.add_option("-e","--e", dest="e", help="Path to annotation e",default="./")
parser.add_option("-f","--f", dest="f", help="Path to annotation f",default="./")
parser.add_option("-g","--g", dest="g", help="Path to annotation g",default="./")
parser.add_option("--hh", dest="hh", help="Path to annotation h",default="./")
parser.add_option("-i","--i", dest="i", help="Path to annotation i",default="./")
parser.add_option("-j","--j", dest="j", help="Path to annotation j",default="./")
#parser.add_option("-k","--k", dest="k", help="Path to annotation k",default="./")
#parser.add_option("-l","--l", dest="l", help="Path to annotation l",default="./")
#parser.add_option("-m","--m", dest="m", help="Path to annotation m",default="./")


(options, args) = parser.parse_args()

# Creates list merged annotation files 
annotated_list = []
for fn in os.listdir(options.annotated):
	if fn.endswith('.tsv'):
		annotated_list.append(fn)
annotated_list = sorted(annotated_list)

#Create list of annotation files
repeattype_list = []
for fn in os.listdir(options.a):
	if fn.startswith('repeattype'):
		repeattype_list.append(fn)

TFBS_list = []
for fn in os.listdir(options.b):
	if fn.startswith('TFBS'):
		TFBS_list.append(fn)

miRNA_list = []
for fn in os.listdir(options.c):
	if fn.startswith('miRNApos'):
		miRNA_list.append(fn)

remap_list = []
for fn in os.listdir(options.d):
	if fn.startswith('remap'):
		remap_list.append(fn)

crm_list = []
for fn in os.listdir(options.e):
	if fn.startswith('CRM'):
		crm_list.append(fn)

pC_list = []
pP_list = []
for fn in os.listdir(options.f):
	if fn.endswith('124phastCons.tsv'):
		pC_list.append(fn)
	elif fn.endswith('124phyloP.tsv'):
		pP_list.append(fn)

BG3_list = []
for fn in os.listdir(options.g):
	if fn.startswith('ChromStateBG3'):
		BG3_list.append(fn)

S2_list = []
for fn in os.listdir(options.hh):
	if fn.startswith('ChromStateS2'):
		S2_list.append(fn)

gerpRS_list = []
gerpN_list = []
for fn in os.listdir(options.i):
	if fn.endswith('GERPRS_all.tsv'):
		gerpRS_list.append(fn)
	elif fn.endswith('GERPN_all.tsv'):
		gerpN_list.append(fn)

phylop_list = []
for fn in os.listdir(options.j):
	if fn.endswith('phylop_all.tsv'):
		phylop_list.append(fn)

# Iterate over derived annotated list (over the chromosomes). 
for filen in annotated_list:
	
	# Open in pandas and add labels for y (derived = 1, simulated = 0).
	df_annotated = pd.read_csv(options.annotated + filen, sep='\t', lineterminator='\n')
	df_annotated = df_annotated.drop(['motifECount', 'motifEHIPos', 'motifEScoreChng', 'SIFTcat', 'SIFTval'], axis='columns')


	# Determine chr number. 
	chr_num = filen.split('_')[0]
	print('Working on: Chr' + chr_num)
	
	## Open Annotations as a pd dataframe. 

	# Define annotation file names. 

	for fn in repeattype_list:
		num = fn.split('repeattype_')[1].replace('.tsv', '')
		if chr_num == num:
			repeattype_f = fn
			break

	for fn in TFBS_list:
		num = fn.split('TFBSS_')[1].replace('.txt', '')
		if chr_num == num:
			TFBS_f = fn
			break

	for fn in miRNA_list:
		num = fn.split('miRNApos_')[1].replace('.txt', '')
		if chr_num == num:
			miRNA_f = fn
			break

	for fn in remap_list:
		num = fn.split('remap_chr')[1].replace('.txt', '')
		if chr_num == num:
			remap_f = fn
			break

	for fn in crm_list:
		num = fn.split('CRMpos_')[1].replace('.txt', '')
		if chr_num == num:
			crm_f = fn
			break

	for fn in pC_list:
		num = fn.replace('_124phastCons.tsv', '')
		if chr_num == num:
			pC_f = fn
			break
	for fn in pP_list:
		num = fn.replace('_124phyloP.tsv', '')
		if chr_num == num:
			pP_f = fn
			break

	for fn in BG3_list:
		num = fn.split('ChromStateBG3_')[1].replace('.txt', '')
		if chr_num == num:
			BG3_f = fn
			break

	for fn in S2_list:
		num = fn.split('ChromStateS2_')[1].replace('.txt', '')
		if chr_num == num:
			S2_f = fn
			break

	for fn in gerpRS_list:
		num = fn.replace('_GERPRS_all.tsv', '')
		if chr_num == num:
			gerpRS_f = fn
			break

	for fn in gerpN_list:
		num = fn.replace('_GERPN_all.tsv', '')
		if chr_num == num:
			gerpN_f = fn
			break

	for fn in phylop_list:
		num = fn.replace('_phylop_all.tsv', '')
		if chr_num == num:
			phylop_f = fn
			break


	# Open files. 
	print('Opening files as pandas DF. ')

	df_repeattype = pd.read_csv(options.a + repeattype_f, sep='\t', lineterminator='\n')

	open_TFBS = open(options.b + TFBS_f, 'r')
	TFBS_json = json.load(open_TFBS)

	open_miRNA = open(options.c + miRNA_f, 'r')
	miRNA_json = json.load(open_miRNA)	
	
	df_remap = pd.read_csv(options.d + remap_f, sep='\s+', lineterminator='\n')

	open_crm = open(options.e + crm_f, 'r')
	crm_json = json.load(open_crm)

	df_pC = pd.read_csv(options.f + pC_f, sep='\t', lineterminator='\n')
	df_pC.drop_duplicates(subset='start_position', keep="last", inplace=True)

	df_pP = pd.read_csv(options.f + pP_f, sep='\t', lineterminator='\n')
	df_pP.drop_duplicates(subset='start_position', keep="last", inplace=True)

	df_BG3 = pd.read_csv(options.g + BG3_f, sep='\t', lineterminator='\n')

	df_S2 = pd.read_csv(options.hh + S2_f, sep='\t', lineterminator='\n')

	df_gerpRS = pd.read_csv(options.i + gerpRS_f, sep='\t', lineterminator='\n')

	df_gerpN = pd.read_csv(options.i + gerpN_f, sep='\t', lineterminator='\n')

	df_phylop = pd.read_csv(options.j + phylop_f, sep='\t', lineterminator='\n')

	## Merge annotations.
	# Merge DataFrames based on the positions of the variants.
	print('Merging DataFrames')

	df_repeattype = df_repeattype.drop("Chr", axis='columns')
	df_merged = pd.merge(df_annotated, df_repeattype, on='Pos', how='left')
	df_merged['RepeatType'].fillna('No_repeat', inplace=True)		

	df_merged['TFBS'] = df_merged['Pos'].isin(TFBS_json).astype(int)

	df_merged['miRNA'] = df_merged['Pos'].isin(miRNA_json).astype(int)

	df_merged = pd.merge(df_merged, df_remap, on='Pos', how='left')

	df_merged['CRM'] = df_merged['Pos'].isin(crm_json).astype(int)

	df_pC = df_pC.drop("#chr_num", axis='columns')
	df_pC.rename(columns = {'start_position':'Pos'}, inplace = True)
	df_merged = pd.merge(df_merged, df_pC, on='Pos', how='left')
	df_merged.rename(columns = {'score':'PhastCons124'}, inplace = True)

	df_pP = df_pP.drop("#chr_num", axis='columns')
	df_pP.rename(columns = {'start_position':'Pos'}, inplace = True)
	df_merged = pd.merge(df_merged, df_pP, on='Pos', how='left')
	df_merged.rename(columns = {'score':'PhyloP124'}, inplace = True)

	df_merged = pd.merge(df_merged, df_BG3, on='Pos', how='left')
	df_merged.rename(columns={'Cat': 'BG3_state'}, inplace = True)

	df_merged = pd.merge(df_merged, df_S2, on='Pos', how='left')
	df_merged.rename(columns={'Cat': 'S2_state'}, inplace = True)

	df_gerpRS = df_gerpRS.drop("Chrom", axis='columns')
	df_merged = pd.merge(df_merged, df_gerpRS, on='Pos', how='left')
	
	df_gerpN = df_gerpN.drop("Chrom", axis='columns')
	df_merged = pd.merge(df_merged, df_gerpN, on='Pos', how='left')

	df_phylop = df_phylop.drop("Chrom", axis='columns')
	df_merged = pd.merge(df_merged, df_phylop, on='Pos', how='left')	

	# Write new dataframe to output.
	df_merged.to_csv(str(chr_num) + '_merged_features.tsv', index=None, sep='\t')