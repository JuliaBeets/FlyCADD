'''
:Author: Julia Beets
:Date: 27-05-2024
:Usage: python <script.py> -i <Path to data>

The script scales every feature in the data with their stdev.

It first merges the simulated and derived datasets, afterwards it scales the features.
Empty features (that only contains 0) are removed from the DF.
Eventually it outputs a merged scaled csv file and a dict file that contains the features as key and its stdev (used for scaling) as value.

'''

# Import dependencies.
import sys,os
from optparse import OptionParser
import json
import pandas as pd


# OptionParser for input. 
parser = OptionParser()
parser.add_option("-i", "--input", dest="input", help="Path to merged Dataframes", default= "./")

(options, args) = parser.parse_args()

# Function for scaling the features with their stdev.
def stdev_norm(df, column):
    col = df[column]
    mean = col.mean()
    stdev = col.std()
    df[column] = (col - mean) / stdev  # Standard scaler
    return mean, stdev

# Sorted list of files in directory.
simulated_files_list = []
derived_files_list = []
for fn in os.listdir(options.input):
    print(fn)
    if 'simulated' in fn:
        simulated_files_list.append(fn)
    if 'derived' in fn:
        derived_files_list.append(fn)
simulated_files_list = sorted(simulated_files_list)

# Iterate over sorted list.
print('Perform scaling!')
for i in range(len(simulated_files_list)):
    simulated_fn = simulated_files_list[i]
    derived_fn = derived_files_list[i]
    
    # Open DataFrames.
    print('Opening csv files.')
    open_simulated = pd.read_csv(os.path.join(options.input, simulated_fn), sep=',')
    open_derived = pd.read_csv(os.path.join(options.input, derived_fn), sep=',')
    
    # Merge DF.
    print('Merging DFs.')
    frames = [open_simulated, open_derived]
    concat_df = pd.concat(frames)
    
    ## Perform scaling on DF and output scaled DF.
    print('Perform scaling.')
    ftr_stdev_dict = dict()
    ftr_mean_dict = dict()  # Dictionary to store mean values

    # Scale DF and append data to dict (key: feature, value: used stdev).
    # The label, chromosome number, position are not scaled. 
    features = list(concat_df.columns)[3:]

    for ftr in features:
        mean, stdev = stdev_norm(concat_df, ftr)
        ftr_stdev_dict[ftr] = stdev
        ftr_mean_dict[ftr] = mean

    # Drop empty features.
    print('Dropping empty features.')
    nan_value = float("NaN")
    concat_df.replace("", nan_value, inplace=True)
    concat_df.dropna(how='all', axis=1, inplace=True)

    columns_to_drop = ['Pos', '#Chrom']
    concat_df.drop(columns=columns_to_drop, inplace=True)
    
    # Output scaled DF and dict.
    print('Generate output.')
    concat_df.to_csv('concat_scaled.csv', index=False)
    
    with open('stdev_dict.txt', 'w') as dict_file:
        dict_file.write(json.dumps(ftr_stdev_dict))

    with open('means_trainingset.txt', 'w') as mean_file:
        mean_file.write(json.dumps(ftr_mean_dict))
