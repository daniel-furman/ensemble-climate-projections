#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 11:26:05 2020

@author: danielfurman
"""

# Now that we have trained our ML classifiers, we are ready to deploy the .pkl
# files we saved from PyCaret's training and examine the model predictions
# on the validation set and examine the models' performance. We first print 
# validation set accuracy, as well as the F statistic and the 2x2 confusion
# matrix. Finally, we visualize the AUC statistic with a ROC curve for each
# model. We also examine the results of a blended model constructed from the 
# aforementioned five most predictive learners.

from sklearn.metrics import confusion_matrix
from matplotlib import pyplot as plt
from matplotlib import style
from sklearn.metrics import roc_curve, auc, f1_score
from sklearn.neural_network import MLPClassifier
from pycaret.classification import *
import warnings
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

# create a dictionary of ML models, name -> (line format, classifier)
# models deployed from ML_sdms_train.py

CLASS_MAP = {
    '5-model Blend': ('-.', load_model('xant_blended')[23]),
    'Random Forest': ('-', load_model('xant_rf')[23]),
    'Extra Trees': ('--', load_model('xant_etrees')[23]),
    'XGBoost': ('-.', load_model('xant_xgb')[23]),
    'LGBoost Machine': ('-.', load_model('xant_lgbm')[23]),
    'Catboost': ('-.', load_model('xant_cboost')[23])
    }

# load training (80%) and test (20%) sets

env_data = pd.read_csv(
    '/Users/danielfurman/Data_science_code/xantusia-codebase/data/envtrain_xv.csv')
training_class = env_data['pa']
env_data = env_data.drop(['Unnamed: 0'], axis=1)
training_data = env_data.drop(['pa'], axis=1) # drop class

env_data_test = pd.read_csv(
    '/Users/danielfurman/Data_science_code/xantusia-codebase/data/envtest_xv.csv')
validation_class = env_data_test['pa']
env_data_test = env_data_test.drop(['Unnamed: 0'], axis=1)
validation_data  = env_data_test.drop(['pa'], axis=1)

# perform validation set analyses, iterate over dictionary :

f_score = np.zeros(len(CLASS_MAP)) # initialize F score array
col_names = [] # col_names for F score
i = 0  # iterator
style.use('default')
plt.rcParams["figure.figsize"] = (6, 4)

for name, (line_fmt, model) in CLASS_MAP.items():
    result = model.fit(training_data, training_class)
    # array w one col per label
    preds = model.predict_proba(validation_data)
    pred = pd.Series(preds[:, 1])
    fpr, tpr, thresholds = roc_curve(validation_class, pred)
    auc_score = auc(fpr, tpr) 
    label = '%s: auc=%.4f' % (name, auc_score)
    plt.plot(fpr, tpr, line_fmt,
             linewidth=2, label=label)
    # compute confusion matrix and F1 stat
    predicted_class_type = model.predict(validation_data)
    print('\n\nFraction correct validation ' + name + ' :',
          np.sum(predicted_class_type == validation_class)
          / len(validation_class))
    cnf_matrix_test = confusion_matrix(validation_class, predicted_class_type)
    print(cnf_matrix_test)
    print('The F1 validation score is : ',
          f1_score(validation_class, predicted_class_type))
    f_score[i] = f1_score(validation_class, predicted_class_type)
    col_names.append(name)
    i = + (i + 1)

# annotate AUC Plot     
plt.legend(loc="lower right", shadow=True)
plt.title('Comparing Classifiers: Validation AUC')
# plt.plot([0, 1], [0, 1], 'k-', alpha = 0.3) #x=y line. Visual aid
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.savefig('auc.png', dpi=400)

# create pandas dataframe with F statistic scores
f_score = pd.DataFrame(data=f_score.reshape(-1, len(f_score)),
                       columns=col_names)
f_score = f_score.rename(index={0: "F-statistic :"})
f_score = f_score.sort_values(by="F-statistic :", axis=1,
                              ascending=False)
