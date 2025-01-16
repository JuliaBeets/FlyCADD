"""
:Author: Julia Beets
:Date: 29-05-2024
:Usage: python <script.py> -i <Path to scaled Dataframes>

Iterates over the scaled DataFrames.
For each DataFrame, it splits the data into a training and testing set.
The data is then fitted to the parameters of a logistic classifier, with the L2 norm penalization.
Afterward, the fitted model is used to predict the test set,
followed by an evaluation of the model on this test set.
Cross-validation is applied to determine the optimal L2 penalty.
"""

# Import dependencies.
import turicreate as tc
import sys, os
from optparse import OptionParser
import csv
import numpy as np

# OptionParser for input.
parser = OptionParser()
parser.add_option("-i", "--input", dest="input", help="Path to scaled Dataframes", default="./")

(options, args) = parser.parse_args()

# Sorted list of files in directory.
input_list = []
for fn in os.listdir(options.input):
    if 'concat_scaled.csv' in fn:
        input_list.append(fn)
input_list = sorted(input_list)

# Function for cross-validation.
def cross_validate(data, target, l2_penalty_values, folds=5):
    results = []
    for l2_penalty in l2_penalty_values:
        accuracy = []
        for i in range(folds):
            train_data, valid_data = data.random_split(0.8)
            model = tc.logistic_classifier.create(train_data, target=target, l2_penalty=l2_penalty, validation_set=None, verbose=False)
            evaluation = model.evaluate(valid_data)
            accuracy.append(evaluation['accuracy'])
        mean_accuracy = np.mean(accuracy)
        results.append((l2_penalty, mean_accuracy))
    return results

# Iterate over sorted list.
print('Iterating over input files')
for fn in input_list:
    # Load the data.
    print('Loading data')
    loaded_data = tc.SFrame(options.input + fn)

    # Split data.
    print('Create train and test set.')
    train_data, test_data = loaded_data.random_split(0.90)

    # Define target column and L2 penalty values for grid search.
    target_column = '#Derived(1)_Simulated(0)'
    l2_penalty_values = [0.01, 0.1, 1.0, 10.0, 100.0]

    # Perform cross-validation to find the optimal L2 penalty.
    print('Performing cross-validation.')
    cv_results = cross_validate(train_data, target_column, l2_penalty_values)
    best_l2_penalty = max(cv_results, key=lambda item: item[1])[0]
    print(f'Best L2 penalty: {best_l2_penalty}')

    # Create model with the optimal L2 penalty.
    print('Creating model with the optimal L2 penalty.')
    model = tc.logistic_classifier.create(train_data, feature_rescaling=True, target=target_column, solver='newton', l1_penalty=0.0, l2_penalty=best_l2_penalty, max_iterations=100, verbose=True)

    # Get coefficients of features.
    top_coefs = model.summary
    print(top_coefs)
    coefs = model.coefficients
    coefs.export_csv(f'model_{fn}_all_coefficients.csv', delimiter=',', line_terminator='\n', header=True)
    print(coefs)

    # Predict on test set.
    print('Predicting on test set.')
    predictions = model.classify(test_data)
    print(predictions)

    # Evaluate model on test set.
    print('Evaluating model.')
    results = model.evaluate(test_data)
    print(results['accuracy'])
    print(results)
    with open(f'model_{fn}_evaluation.csv', 'w') as f:
        for key in results.keys():
            f.write("%s,%s\n" % (key, results[key]))

    # Save model (saves to a folder).
    model.save(f'model_{fn}_file')
