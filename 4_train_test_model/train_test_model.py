'''
:Author: Seyan Hu
:Date: 10-1-2023
:Usage: python <script.py> -i <Path to scaled Dataframes>


Iterates over the scaled DataFrames.
For each DataFrame it splits the data into a training and testing set. 
The data is then fitted to the parameters of a logistic classifer, with the L2 norm penalization.
Afterwards the fitted model is used to predict the test set, 
followed by an evaluation of the model on this test set.

'''

# Import dependencies.
import turicreate as tc
import sys,os
from optparse import OptionParser
import csv


# OptionParser for input. 
parser = OptionParser()
parser.add_option("-i", "--input", dest="input", help="Path to scaled Dataframes", default= "./")

(options, args) = parser.parse_args()


# Sorted list of files in directory.
input_list = []
for fn in os.listdir(options.input):
	if 'concat_scaled.csv' in fn:
		input_list.append(fn)
input_list = sorted(input_list)


# Iterate over sorted list.
print('Iterating over input files')
for fn in input_list:
	
	
	# Load the data.
	print('Loading data')
	loaded_data = tc.SFrame(options.input + fn)
	

	# Split data.
	print('Create train and test set.')
	train_data, test_data = loaded_data.random_split(0.90)


	# Create model.
	print('Creating model.')
	model = tc.logistic_classifier.create(train_data, feature_rescaling=True, target='#Derived(1)_Simulated(0)', solver='newton',l1_penalty=0.0, l2_penalty=1.0 , max_iterations=100, verbose=True)

	# Get coefficients of features. 
	top_coefs = model.summary
	print(top_coefs)
	coefs = model.coefficients
	feature_num = int(test_data.num_columns())
	coefs.export_csv('model_all_coefficients.csv', delimiter=',', line_terminator='\n', header=True)
	print(coefs)
	

	# Predict on test set. 
	print('Predict on test set.')
	predictions = model.classify(test_data)
	print(predictions)
	
	# Evaluate model on test set. 
	print('Evaluate model.')
	results = model.evaluate(test_data)
	print(results['accuracy'])
	print(results)
	with open('model_evaluation.csv', 'w') as f:
		for key in results.keys():
			f.write("%s,%s\n"%(key, results[key]))
	
	
	# Save model (saves to a folder).
	model.save('model_file')
