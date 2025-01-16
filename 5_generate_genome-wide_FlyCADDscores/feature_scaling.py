'''
:Author: Julia Beets
:Date: 28-05-2024
:Usage: python <script.py> -i <Path to merged Dataframes> -d <Path to the directory containing means and stdev files>

'''

# Import dependencies.
import sys, os
from optparse import OptionParser
import json
import pandas

parser = OptionParser()
parser.add_option("-i", "--input", dest="input", help="Path to merged Dataframes", default="./")
parser.add_option("-d", "--dict_p", dest="dict_p", help="Path to the directory containing means and stdev files", default="../train_test_model/output/dir_scaled_df/")

(options, args) = parser.parse_args()

def scale_feature(df, column, col_mean, col_stdev):
    col = df[column]
    df[column] = (col - col_mean) / col_stdev

# Load the means and stdev dicts.
with open(os.path.join(options.dict_p, 'means_trainingset.txt'), 'r') as means_file:
    means_dict = json.load(means_file)

with open(os.path.join(options.dict_p, 'stdev_dict.txt'), 'r') as stdev_file:
    stdev_dict = json.load(stdev_file)

files_list = [fn for fn in os.listdir(options.input) if '_enc_mv.csv' in fn]
files_list = sorted(files_list)

print('Perform scaling!')
for fn in files_list:
    print('Working on: ' + fn)

    print('Determine chr number and chunk number.')
    chr_num, chunk_num = fn.split('_')[:2]

    print('Opening csv file.')
    open_df = pandas.read_csv(os.path.join(options.input, fn), sep=',')

    # Perform scaling on DataFrame and output scaled DataFrame.
    print('Perform scaling.')
    features = list(open_df.columns)[2:]
    for ftr in features:
        if ftr in means_dict and ftr in stdev_dict:
            scale_feature(open_df, ftr, means_dict[ftr], stdev_dict[ftr])

    print('Dropping empty features.')
    nan_value = float("NaN")
    open_df.replace("", nan_value, inplace=True)
    open_df.dropna(how='all', axis=1, inplace=True)

    print('Generate output.')
    output_file = os.path.join(chr_num + '_' + chunk_num + '_scaled.csv')
    open_df.to_csv(output_file, index=False)

print('Scaling complete!')
