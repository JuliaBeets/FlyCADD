'''
:Author: Seyan Hu
:Date: 23-1-2023
:Usage: python <script.py> -i <Path to prediction output>


Iterates over the chunks and processes the Dataframe by removing the probability for class 1 (benign). 
Since only class 0 is needed.

'''

# Import dependencies.
import sys,os
from optparse import OptionParser
import csv
import pandas as pd


# OptionParser for input. 
parser = OptionParser()
parser.add_option("-i", "--input", dest="input", help="Path to prediction output", default= "./")

(options, args) = parser.parse_args()


# Sorted list of files in directory.
input_list = []
for fn in os.listdir(options.input):
	if '.csv' in fn:
		input_list.append(fn)
input_list = sorted(input_list)


# Iterate over sorted list.
print('Iterating over input files')
for fn in input_list:

	# Determine chr number and chunk.
	chr_num = fn.split('_')[0]
	chunk_num = fn.split('_')[1]
	
	opened_df = pd.read_csv(options.input+fn, header=None)
	
	
	# Convert first column of df to list.
	series = list(opened_df.iloc[:,0])
	series = [item.strip("array('d', [").strip("])") for item in series]

	
	# Define proba of class 0 (deleterious). 
	class_0_list = []
	for proba in series:
		class_0 = proba.split(', ')[0]

		# Convert string to float. 
		if 'e-' in class_0:
			num = class_0.split('e-')[0]
			power = class_0.split('e-')[1] # *10^-
			class_0 = float(num) * 10**-float(power)
		else:
			class_0 = float(class_0)
		class_0_list.append(class_0)
	del series

	
	# List to DataFrame and merge on index. 	
	opened_df.pop(opened_df.columns[0])
	df_class_0 = pd.DataFrame (class_0_list, columns = ['proba_class_0'])
	joined_df = pd.merge(opened_df, df_class_0, left_index=True, right_index=True)
	joined_df.to_csv(chr_num + '_' + chunk_num + '_processed_pred.csv', sep=',', encoding='utf-8', index=False, header=False)