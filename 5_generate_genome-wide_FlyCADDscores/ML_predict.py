'''
:Author: Seyan Hu
:Date: 23-1-2023
:Usage: python <script.py> -i <Path to scaled Dataframes> -m <Path to saved models>

:Modified: Julia Beets, 30-05-2024

Iterates over the scaled DataFrames.
For each DataFrame its previously trained model is loaded and the variants posterior likelihoods are predicted. 
Hardcoded model name!
'''

# Import dependencies.
import turicreate as tc
import sys,os
from optparse import OptionParser
import csv
import pandas as pd


# OptionParser for input. 
parser = OptionParser()
parser.add_option("-i", "--input", dest="input", help="Path to scaled Dataframes", default= "./")
parser.add_option("-m", "--model", dest="model", help="Path to saved models", default= "./")

(options, args) = parser.parse_args()


# Sorted list of files in directory.
input_list = []
for fn in os.listdir(options.input):
	if '_scaled.csv' in fn:
		input_list.append(fn)
input_list = sorted(input_list)
print(str(input_list))

# Iterate over sorted list.
print('Iterating over input files')
for fn in input_list:
	
	
	## Determine chr number and chunk.
	chr_num = fn.split('_')[0]
	chunk_num = fn.split('_')[1]
	sub_chunk = list(str(chunk_num))[0]
	
	
	## Load the data.
	print('Loading data for chr: ' + str(chr_num))
	loaded_data =  tc.SFrame(options.input + fn)
	loaded_data = loaded_data.remove_column('#Chrom')
	
	## Load model.
	loaded_model = tc.load_model(options.model + 'WG_cross_model_file')
	
	
	## Predict posterior likelihoods (left: proba for class 0 simulated; right proba for class 1 derived).
	predictions = loaded_model.predict(loaded_data, output_type = 'probability_vector')
	predictions =  tc.SFrame(predictions)


	## Get extract correct columns from sframe of loaded input.
	df_pred = predictions.to_dataframe()
	
	
	## Merge predictions and their corresponding positions, ref nt and alt nt. 
	if chunk_num == 'aa':
		df_opened = pd.read_csv('output/dir_chunked_df/dir_a/'+chr_num+'_chunk.aa', sep='\t', header=0)
		df_opened = df_opened[['Pos', 'Ref', 'Alt']]
	else:
		df_opened = pd.read_csv('output/dir_chunked_df/dir_'+sub_chunk+'/'+chr_num+'_chunk_ln.'+chunk_num, sep='\t', header=0)
		df_opened = df_opened[['Pos', 'Ref', 'Alt']]
	
	joined_df = pd.merge(df_pred, df_opened, left_index=True, right_index=True)

	
	## Save the predictions into csv format.
	joined_df.to_csv(chr_num + '_' + chunk_num + '_predictions.csv', sep=',', encoding='utf-8', index=False, header=False)

