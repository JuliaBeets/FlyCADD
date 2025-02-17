# Training and testing of the machine learning model
The fourth step in generating the FlyCADD is training and testing of the machine learning (ML) model. The ML model is a logistic regression classifier of TuriCreate, trained on 90 % of the training dataset and tested with the remaining 10 % of the training (performance evaluation). The model is trained to predict the posterior probability of a variant being deleterious.

### Feature scaling
Processed (annotated, imputed and hot-encoded) derived and simulated variants per chromosome from the previous step are first merged into one file per chromosome and then scaled according to their standard deviation. In addition, empty features are removed from the dataframe. <br />
Usage: `python feature_scaling.py -i <path to directory containing processed variants>`

### Training and testing
Cross-validation was apply to determine the L2 penalization for training. <br />
Usage: `python <script.py> -i <Path to scaled Dataframes>`

The sets of variants are split into a training dataset (90%) and a testing dataset (10%). This data is fitted to the parameters of a logistic classifier (TuriCreate) with the L2 norm penalization. The fitted model is tested by predicting the class (derived/neutral (1) or simulated/deleterious (0)) and the performance of the model on this test dataset is evaluated. <br/>
Usage: `python train_test_model.py -i <path to directory with scaled datasets per chromosome from the previous step>`