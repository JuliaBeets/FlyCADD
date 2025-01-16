#!/usr/bin/env python
# -*- coding: ASCII -*-

"""
:Author: Christian Gross
:Contact: c.gross@tudelft.nl
:Date: 15-08-18

This script takes fully annotated variant files. 
The script encodes and imputes data of features. 

Instead of reading in and manipulating both files at the same time it will be performed one at the time, 
but the simulated data first because mean imputation is used from the simulated data.
It check automatically if imputation file is present in current directory, if not and if derived is specified then an error is thrown.

:Edited by: Seyan Hu
:Date: 13-12-2022
:Usage: python data_encoding_mv_handling.py -i <Merged annotation file> -p <Path to the imputation file>

:Modified: Julia Beets
:Date: 02-01-2024

"""

# Import dependencies
print('Import dependencies')
import sys, os
import numpy as np
import math
from optparse import OptionParser
import pandas
from pandas import HDFStore,DataFrame
import ast
import time


# OptionParser for input. 
print('Load OptionParser')
parser = OptionParser()
parser.add_option("-i", "--input", dest="input", help="annotated infile containing simulated variants",default="./19_test_inp_features.tsv")
parser.add_option("-d","--derived", dest="derived", help="Flag, if given it will relabel Ref and Alt. (def OFF) (can be ignored for encoding/imputation when generating genome-wide CADD scores)",default=False,action="store_true")
parser.add_option("-p", "--impdict", dest="impdict", help="Path to imputation dict file")

(options, args) = parser.parse_args()
	
	
# Function for checking if classes are properly encoded. 
def class_encoded_check(classlabel, selection, data):
	for i in selection: 
		try:
			data[classlabel+'_'+i]
		except KeyError:
			data[classlabel+'_'+i] = np.zeros(data.shape[0],dtype=float)
	return data
	
	
# Determine chr number.
print('Chr number:')
chr_num = options.input.split('/')[-1].split('_')[0]


# Open DF in pandas. 
print('Open DF')
myData = pandas.read_csv(options.input, sep='\t', na_values=['-'], dtype={"#Chrom":object, "Pos":int, "Ref":object , "Alt":object , "isTv":bool, "Consequence":object , "GC":float,  "CpG":float,"Domain":object,"oAA":object,"nAA":object,"Grantham":float,"cDNApos":float,"relcDNApos":float,"CDSpos":float, "relCDSpos":float,"protPos":float,"relprotPos":float, "PhastCons":float,"PhyloP":float, "Roll":float, "EP":float,"MGW":float,"HelT":float,"ProT":float,"Repeats":int, "ingene":int, "RepeatType":object, "TFBS": int, "miRNA":int, "ReMap": float, "CRM": int, "PhastCons124": float, "PhyloP124": float, "BG3_state": float, "S2_state":float, "GERPRS_all":float, "GERPN_all":float, "phylop_all": float})

# Set motifEHIPos to 1 if True and else to 0. 
# Special treatment NaN has to be set to 0.0, everything else 1.0
print("Replace NaN values with '0.0'")
#myData['motifEHIPos'] = myData['motifEHIPos'].map({True: 1.0, False: 0.0, np.NaN: 0.0})


# Handles missing values by fillin it with 0.0. 
# Log transform Gerp, but before that I also have to handle the missing values.
#! Scores were not log transformed due to it causing -inf for these scores. 
np.seterr(divide = 'ignore') 



# Other categorical features which need a different category for undefined features
myData['Domain'] = myData['Domain'].map({ 'lcompl':'lcompl', 'ndomain':'ndomain','ncoils':'ncoils','tmhmm':'tmhmm','sigp':'sigp',np.nan:'UD'})
myData['oAA'] = myData['oAA'].map({np.nan:'UD', 'P':'P', 'V':'V', 'G':'G', 'I':'I', 'F':'F', 'E':'E', 'N':'N', 'T':'T', 'Y':'Y', 'L':'L', 'C':'C', 'S':'S','A':'A', 'K':'K', 'H':'H', 'R':'R', 'D':'D', 'Q':'Q', '*':'*', 'M':'M', 'W':'W','U':'U'})
myData['nAA'] = myData['nAA'].map({np.nan:'UD', 'P':'P', 'V':'V', 'G':'G', 'I':'I', 'F':'F', 'E':'E', 'N':'N', 'T':'T', 'Y':'Y', 'L':'L', 'C':'C', 'S':'S','A':'A', 'K':'K', 'H':'H', 'R':'R', 'D':'D', 'Q':'Q', '*':'*', 'M':'M', 'W':'W','U':'U'})
myData['RepeatType'] = myData['RepeatType'].map({ 'DNA':'DNA', 'LINE':'LINE','Low_complexity':'LCR','LTR':'LTR','No_repeat':'No_repeat',np.nan:'UD','RC':'RC','Satellite':'Satellite','Simple_repeat':'Simple_repeat','rRNA':'rRNA'})

# 'isTv' first has to be translated to floats and then missing values imputed by inserting 0.5.
print('Imputation of isTv')
myData['isTv'] = myData['isTv'].astype(float).fillna(0.5)


# Create indicator columns before imputation.
print('Handle imputation for missing values')
for label in ['cDNApos', 'CDSpos', 'protPos', 'Grantham']:
	myData['IND_'+label] = myData[label].isnull().astype(float)


# (Generate &) Load the mean imputations from the simulated variants. 
# According to the hCADD supplementary data the mean imputed values are generated based on the means from the simulated data only, 
# not from the whole data set.
print('Load imputation dict')
selection = ['GC','CpG','PhastCons', 'PhyloP', 'Roll', 'EP', 'MGW', 'HelT', 'ProT', 'PhastCons124', 'PhyloP124', 'phylop_all', 'ReMap', "GERPRS_all", 'GERPN_all']
means = {}
if options.derived==False:
	try:
		with open(options.impdict + chr_num + '_imputation_dict.txt') as f: 
			means = ast.literal_eval(f.read())
	except IOError:
		sys.exit("The imputation_dict.txt file which contains the mean imputed values from the simulated data is missing, the file should be present in the path: '../generate_annotations/output/dir_complete_dataset/'.")
elif options.derived==True:
	try:
		with open(options.impdict + chr_num + '_imputation_dict.txt') as f:
			means = ast.literal_eval(f.read())
	except IOError:
		sys.exit("The imputation_dict.txt file which contains the mean imputed values from the simulated data is missing. This information is necessary vtherefore the program stops prematurely")

print('Perform imputation of missing values')
for label in selection:
	try:
		myData[label] = myData[label].fillna(means[label])
	except KeyError:
		sys.exit("imputation_dict.txt was loaded correctly but does not contain all Keys. %s is missing" %label)


# The mean imputation has finished, derived variants can be arranged in accordence to the specifications mentioned in the supplementary data 
# and both data sets can be joined into one (later). 


# Fil in missing values for these features. 
selection = ['cDNApos', 'relcDNApos', 'CDSpos', 'relCDSpos', 'protPos', 'relprotPos','Grantham', "BG3_state", "S2_state", "GERPRS_all", "GERPN_all"]
for label in selection:
	myData[label] = myData[label].fillna(0.0)


# Output DF with imputation applied. 
chunk_num = options.input.split('.')[1]


# Encoding via pandas.get_dummies for these features. 
print('Handle one-hot encoding')
for label in ['Ref','Alt','Consequence','Domain','oAA','nAA', 'RepeatType']:
	myData = pandas.get_dummies(myData,prefix=label,columns=[label])


# Checking if all possible classes are encoded even if they are not available in this particular data set.
print('Check encoding')
myData = class_encoded_check('Ref',['A','C','G','T'],myData)
myData = class_encoded_check('Alt',['A','C','G','T'],myData)
myData = class_encoded_check('oAA',['*', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'UD'],myData)
myData = class_encoded_check('nAA',['*', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'UD'],myData)
myData = class_encoded_check('Consequence',['U3', 'U5', 'DN', 'IG', 'I', 'NC', 'NS', 'R', 'CS', 'S', 'SG', 'SL', 'SN', 'UP','O'],myData)
myData = class_encoded_check('Domain',['ncoils', 'tmhmm','sigp','lcompl','UD','ndomain'],myData)
myData = class_encoded_check('RepeatType',['DNA', 'LINE','LCR','LTR','No_repeat','UD','RC','Satellite','Simple_repeat','rRNA'],myData)


# Add data to encoded ref, alt nucleotides. 
# Label all nucleotides with .1 therefore it is not necessary anymore to check for duplicated column names.
# Exclude nucleotide substitutions in which the reference and alternative are the same because these do not exist.
# Label all substitutions as "Nuc.1_Nuc.1", then it is not needed to save nucleotide substitutions 
# and compare them with aminoacid substitutions to avoid duplicates.
# Ref x Alt
# A.1_C.1,A.1_G.1,A.1_T.1,C.1_A.1,C.1_G.1,C.1_T.1,G.1_A.1,G.1_C.1,G.1_T.1,T.1_A.1,T.1_C.1,T.1_G.1
print('Encoding for nt substitutions')
start1 = time.time()
selection = ['A','C','G','T']
Alt1 = pandas.DataFrame()
for i in selection:
	Ref = myData['Ref_'+i].values
	for j in selection:
		if i==j:
			continue
		else:
			Alt = myData['Alt_'+j].values
			Alt1[i+'.1_'+j+'.1'] = pandas.Series([o*a for o,a in zip(Ref,Alt)],index=myData.index)

myData = pandas.concat([myData,Alt1],axis=1)
del Alt1


# Add data to encoded AA substitutions. 
# oAA x nAA
# ['*_*', '*_A', '*_C', '*_D', '*_E', '*_F', '*_G', '*_H', '*_I', '*_K', '*_L', '*_M', '*_N', '*_P', '*_Q', '*_R', '*_S', '*_T', '*_V', '*_W', '*_Y', 'A_*', 'A_A', 'A_C', 'A_D', 'A_E', 'A_F', 'A_G', 'A_H', 'A_I', 'A_K', 'A_L', 'A_M', 'A_N', 'A_P', 'A_Q', 'A_R', 'A_S', 'A_T', 'A_V', 'A_W', 'A_Y', 'C_*', 'C_A', 'C_C', 'C_D', 'C_E', 'C_F', 'C_G', 'C_H', 'C_I', 'C_K', 'C_L', 'C_M', 'C_N', 'C_P', 'C_Q', 'C_R', 'C_S', 'C_T', 'C_V', 'C_W', 'C_Y', 'D_*', 'D_A', 'D_C', 'D_D', 'D_E', 'D_F', 'D_G', 'D_H', 'D_I', 'D_K', 'D_L', 'D_M', 'D_N', 'D_P', 'D_Q', 'D_R', 'D_S', 'D_T', 'D_V', 'D_W', 'D_Y', 'E_*', 'E_A', 'E_C', 'E_D', 'E_E', 'E_F', 'E_G', 'E_H', 'E_I', 'E_K', 'E_L', 'E_M', 'E_N', 'E_P', 'E_Q', 'E_R', 'E_S', 'E_T', 'E_V', 'E_W', 'E_Y', 'F_*', 'F_A', 'F_C', 'F_D', 'F_E', 'F_F', 'F_G', 'F_H', 'F_I', 'F_K', 'F_L', 'F_M', 'F_N', 'F_P', 'F_Q', 'F_R', 'F_S', 'F_T', 'F_V', 'F_W', 'F_Y', 'G_*', 'G_A', 'G_C', 'G_D', 'G_E', 'G_F', 'G_G', 'G_H', 'G_I', 'G_K', 'G_L', 'G_M', 'G_N', 'G_P', 'G_Q', 'G_R', 'G_S', 'G_T', 'G_V', 'G_W', 'G_Y', 'H_*', 'H_A', 'H_C', 'H_D', 'H_E', 'H_F', 'H_G', 'H_H', 'H_I', 'H_K', 'H_L', 'H_M', 'H_N', 'H_P', 'H_Q', 'H_R', 'H_S', 'H_T', 'H_V', 'H_W', 'H_Y', 'I_*', 'I_A', 'I_C', 'I_D', 'I_E', 'I_F', 'I_G', 'I_H', 'I_I', 'I_K', 'I_L', 'I_M', 'I_N', 'I_P', 'I_Q', 'I_R', 'I_S', 'I_T', 'I_V', 'I_W', 'I_Y', 'K_*', 'K_A', 'K_C', 'K_D', 'K_E', 'K_F', 'K_G', 'K_H', 'K_I', 'K_K', 'K_L', 'K_M', 'K_N', 'K_P', 'K_Q', 'K_R', 'K_S', 'K_T', 'K_V', 'K_W', 'K_Y', 'L_*', 'L_A', 'L_C', 'L_D', 'L_E', 'L_F', 'L_G', 'L_H', 'L_I', 'L_K', 'L_L', 'L_M', 'L_N', 'L_P', 'L_Q', 'L_R', 'L_S', 'L_T', 'L_V', 'L_W', 'L_Y', 'M_*', 'M_A', 'M_C', 'M_D', 'M_E', 'M_F', 'M_G', 'M_H', 'M_I', 'M_K', 'M_L', 'M_M', 'M_N', 'M_P', 'M_Q', 'M_R', 'M_S', 'M_T', 'M_V', 'M_W', 'M_Y', 'N_*', 'N_A', 'N_C', 'N_D', 'N_E', 'N_F', 'N_G', 'N_H', 'N_I', 'N_K', 'N_L', 'N_M', 'N_N', 'N_P', 'N_Q', 'N_R', 'N_S', 'N_T', 'N_V', 'N_W', 'N_Y', 'P_*', 'P_A', 'P_C', 'P_D', 'P_E', 'P_F', 'P_G', 'P_H', 'P_I', 'P_K', 'P_L', 'P_M', 'P_N', 'P_P', 'P_Q', 'P_R', 'P_S', 'P_T', 'P_V', 'P_W', 'P_Y', 'Q_*', 'Q_A', 'Q_C', 'Q_D', 'Q_E', 'Q_F', 'Q_G', 'Q_H', 'Q_I', 'Q_K', 'Q_L', 'Q_M', 'Q_N', 'Q_P', 'Q_Q', 'Q_R', 'Q_S', 'Q_T', 'Q_V', 'Q_W', 'Q_Y', 'R_*', 'R_A', 'R_C', 'R_D', 'R_E', 'R_F', 'R_G', 'R_H', 'R_I', 'R_K', 'R_L', 'R_M', 'R_N', 'R_P', 'R_Q', 'R_R', 'R_S', 'R_T', 'R_V', 'R_W', 'R_Y', 'S_*', 'S_A', 'S_C', 'S_D', 'S_E', 'S_F', 'S_G', 'S_H', 'S_I', 'S_K', 'S_L', 'S_M', 'S_N', 'S_P', 'S_Q', 'S_R', 'S_S', 'S_T', 'S_V', 'S_W', 'S_Y', 'T_*', 'T_A', 'T_C', 'T_D', 'T_E', 'T_F', 'T_G', 'T_H', 'T_I', 'T_K', 'T_L', 'T_M', 'T_N', 'T_P', 'T_Q', 'T_R', 'T_S', 'T_T', 'T_V', 'T_W', 'T_Y', 'V_*', 'V_A', 'V_C', 'V_D', 'V_E', 'V_F', 'V_G', 'V_H', 'V_I', 'V_K', 'V_L', 'V_M', 'V_N', 'V_P', 'V_Q', 'V_R', 'V_S', 'V_T', 'V_V', 'V_W', 'V_Y', 'W_*', 'W_A', 'W_C', 'W_D', 'W_E', 'W_F', 'W_G', 'W_H', 'W_I', 'W_K', 'W_L', 'W_M', 'W_N', 'W_P', 'W_Q', 'W_R', 'W_S', 'W_T', 'W_V', 'W_W', 'W_Y', 'Y_*', 'Y_A', 'Y_C', 'Y_D', 'Y_E', 'Y_F', 'Y_G', 'Y_H', 'Y_I', 'Y_K', 'Y_L', 'Y_M', 'Y_N', 'Y_P', 'Y_Q', 'Y_R', 'Y_S', 'Y_T', 'Y_V', 'Y_W', 'Y_Y']
print('Encoding for AA substitutions')
selection = ['*', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
Alt1 = pandas.DataFrame()
start2 = time.time()
for i in selection:
	Ref = myData['oAA_'+i].values
	for j in selection:
		Alt = myData['nAA_'+j].values
		Alt1[i+'_'+j] = pandas.Series([o*a for o,a in zip(Ref,Alt)],index=myData.index)

myData = pandas.concat([myData,Alt1],axis=1)
del Alt1
end2 = time.time()


# Add data to encoded consequences. 
# consequence and set D
print('Encoding for consequences')
selection = ['U3', 'U5', 'CS','DN','I','IG','NC','NS', 'R','S','SG','SL','SN','UP']
selection1 = ['cDNApos', 'CDSpos', 'PhastCons', 'PhyloP', 'protPos', 'relcDNApos', 'relCDSpos', 'relprotPos', 'PhastCons124', 'PhyloP124', 'phylop_all', 'ReMap', "GERPRS_all", 'GERPN_all']
Alt1 = pandas.DataFrame()

start3 = time.time()
for i in selection:
	for j in selection1:
		Ref = myData['Consequence_'+i].values
		Alt = myData[j].values
		Alt1[i+'_'+j] = pandas.Series([o*a for o,a in zip(Ref,Alt)],index=myData.index)

myData = pandas.concat([myData,Alt1],axis=1)
del Alt1


# Output DF to csv file. 
print('Write output')
chunk_num = options.input.split('.')[1]
myData.to_csv(chr_num + '_' + chunk_num + '_enc_mv.csv', index=False, na_rep='-')