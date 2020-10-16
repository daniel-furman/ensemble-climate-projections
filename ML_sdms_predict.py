#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 11:26:05 2020

@author: danielfurman
"""

# Now that we have trained and tuned our ML models, we are ready to predict
# the validation set and examine the models' performance. We first print the
# validation set accuracy, as well as the F statistic and the 2x2 confusion
# matrix. Finally, we visualize the AUC statistic with a ROC curve for each
# model. We also examine the results of a blended model constructed from the 
# afformentioned five most predictive learners.

import numpy as np
import pandas as pd

from sklearn.metrics import confusion_matrix
from matplotlib import pyplot as plt
from matplotlib import style
from sklearn.metrics import roc_curve, auc, f1_score
from sklearn.neural_network import MLPClassifier
from pycaret.classification import *
import sklearn.ensemble as ensemble
from lightgbm import LGBMClassifier
from xgboost import XGBClassifier

env_data = pd.read_csv(
    '/Volumes/HardDrive/xvigilis-data-main/envtrain_corr.csv')
class_type = env_data['pa']
env_data.head()

env_data1 = pd.read_csv(
    '/Volumes/HardDrive/xvigilis-data-main/testbackg_corr.csv')
env_data2 = pd.read_csv(
    '/Volumes/HardDrive/xvigilis-data-main/testpres_corr.csv')

env_data1.insert(0, 'pa', 0)
env_data2.insert(0, 'pa', 1)

env_data_test = pd.concat([env_data1, env_data2])

class_type_test = env_data_test['pa']
env_data = env_data.drop(['Unnamed: 0'], axis=1)
env_data = env_data.drop(['pa'], axis=1)

env_data_test = env_data_test.drop(['Unnamed: 0'], axis=1)
env_data_test = env_data_test.drop(['pa'], axis=1)


#print('head validation features (length = 1282):\n\n', env_data_test.head())
#validation set, 20 percent of the data
#print('\nhead training features (length = 5125):\n\n', env_data.head())
#train set, 80 percent of the data

import warnings
warnings.filterwarnings("ignore")

# create a dictionary of ML models, name -> (line format, classifier)
# params set to those determined in ML_sdms_train.py
# by printing load_model('pycaret-model')


CLASS_MAP = {
'Random Forest':('-', ensemble.RandomForestClassifier(bootstrap=False, ccp_alpha=0.0,
                                        class_weight=None, criterion='entropy',
                                        max_depth=60, max_features='log2',
                                        max_leaf_nodes=None, max_samples=None,
                                        min_impurity_decrease=0.0,
                                        min_impurity_split=None,
                                        min_samples_leaf=1, min_samples_split=9,
                                        min_weight_fraction_leaf=0.0,
                                        n_estimators=40, n_jobs=-1,
                                        oob_score=False, random_state=8143,
                                        verbose=0, warm_start=False)),

'XGBoost': ('-.', XGBClassifier(base_score=0.5, booster='gbtree', colsample_bylevel=1,
              colsample_bynode=1, colsample_bytree=1, gamma=0, gpu_id=-1,
              importance_type='gain', interaction_constraints='',
              learning_rate=0.300000012, max_delta_step=0, max_depth=6,
              min_child_weight=1, monotone_constraints='()',
              n_estimators=100, n_jobs=-1, num_parallel_tree=1,
              objective='binary:logistic', random_state=6289, reg_alpha=0,
              reg_lambda=1, scale_pos_weight=1, subsample=1,
              tree_method='exact', validate_parameters=1, verbosity=0)),

'Extra Trees':('--', ensemble.ExtraTreesClassifier(bootstrap=False, ccp_alpha=0.0,
                                      class_weight=None, criterion='gini',
                                      max_depth=None, max_features='auto',
                                      max_leaf_nodes=None, max_samples=None,
                                      min_impurity_decrease=0.0,
                                      min_impurity_split=None,
                                      min_samples_leaf=1, min_samples_split=2,
                                      min_weight_fraction_leaf=0.0,
                                      n_estimators=100, n_jobs=-1,
                                      oob_score=False, random_state=8143,
                                      verbose=0, warm_start=False)),

'LGBoost Machine':('-.', LGBMClassifier(boosting_type='gbdt',
                                                     class_weight=None,
                                colsample_bytree=1.0, importance_type='split',
                                learning_rate=0.1, max_depth=-1,
                                min_child_samples=20, min_child_weight=0.001,
                                min_split_gain=0.0, n_estimators=100, n_jobs=-1,
                                num_leaves=31, objective=None,
                                random_state=8143, reg_alpha=0.0,
                                reg_lambda=0.0, silent=True, subsample=1.0,
                                subsample_for_bin=200000, subsample_freq=0)),

'MLP neural-net':('--', MLPClassifier(solver='adam'))

    }


training_data = env_data
training_class = class_type

validation_data = env_data_test
validation_class = class_type_test

f_score = np.zeros(len(CLASS_MAP)+1)
i = 0 #iterator

# iterate over dictionary :
style.use('default')
plt.rcParams["figure.figsize"] = (6,4)


for name, (line_fmt, model) in CLASS_MAP.items():

    result = model.fit(training_data, training_class)
    # array w one col per label}|
    preds = model.predict_proba(validation_data)
    
    pred = pd.Series(preds[:,1])
    fpr, tpr, thresholds = roc_curve(validation_class, pred)
    auc_score = auc(fpr, tpr)
    label='%s: auc=%.4f' % (name, auc_score)
    plt.plot(fpr, tpr, line_fmt,
        linewidth=2, label=label)


    #compute confusion matrix and F1 stat    
    predicted_class_type = model.predict(validation_data)
    print('\n\nFraction correct validation ' + name +' :' ,
          np.sum(predicted_class_type == validation_class)
              / len(validation_class))
    cnf_matrix_test = confusion_matrix(validation_class, predicted_class_type)
    print(cnf_matrix_test)
    print('The F1 validation score is : ', 
          f1_score(validation_class, predicted_class_type))
    f_score[i] = f1_score(validation_class, predicted_class_type)
    i =+ (i + 1)

# annotate AUC Plot     
plt.legend(loc="lower right", shadow = True) 
plt.title('Comparing Classifiers: Validation AUC')
plt.plot([0, 1], [0, 1], 'k-', alpha = 0.3) #x=y line. Visual aid 
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate') 
plt.savefig('images-xant/auc.png', dpi = 400)


# finally, print model validation statistics for blended model
env_data1 = pd.read_csv(
    '/Volumes/HardDrive/xvigilis-data-main/testbackg_corr.csv')
env_data2 = pd.read_csv(
    '/Volumes/HardDrive/xvigilis-data-main/testpres_corr.csv')
env_data1.insert(0, 'pa', 0)
env_data2.insert(0, 'pa', 1)
env_data_test = pd.concat([env_data1, env_data2])
env_data_test = env_data_test.drop(['Unnamed: 0'], axis=1)
print('\n')
blender = load_model('xant_blended')
predictions = predict_model(blender, data = env_data_test)
pred_blend = predictions['Label'].astype(str).astype(int)

print('\n\nFraction correct validation ' + 'Blended model' +' :' , 
      np.sum(pred_blend == validation_class)/len(validation_class))
cnf_matrix_test = confusion_matrix(validation_class, pred_blend)
print(cnf_matrix_test)

print('The F1 validation score is : ', f1_score(
    validation_class, pred_blend))

f_score[5] = f1_score(validation_class, pred_blend)
columns = ['RForest', 'XGBoost', 'Extra Trees', 'LGBM', 'MLP-net', 'Blended']
f_score = pd.DataFrame(data = f_score.reshape(-1, len(f_score)),
                       columns=columns)
f_score = f_score.rename(index={0: "F-statistic :"})
f_score = f_score.sort_values(by = "F-statistic :", axis = 1,
                              ascending = False)






