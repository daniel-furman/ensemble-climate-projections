#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 11:26:03 2020

@author: danielfurman

"""

# Binary classification with ten BioClim features. Five are a function of 
# precipitation and five of temperature, all decorrelated below a 0.5
# correlation threshold. We use PyCaret to train and tune (10-fold cv) our
# models from a  train set that contains 80%  of the total data. We find that
# random forest, xgboost, lgbm, catboost, and extra trees algorithms perform
# best. The last step was to create blended model from these five.

import pandas as pd
from pycaret.utils import version
version()
from pycaret.classification import *

data = pd.read_csv("/Users/danielfurman/Data_science_code/xantusia-data-main/xant-pycaret.csv")
data = data.sample(frac=1)
data = data.drop(['Unnamed: 0','Unnamed: 0.1'], axis = 1)
exp_clf = setup(data, target = 'pa');
compare_models();


# model statistics are the mean of 10-fold cross validation:

etrees = create_model('et') #Accuracy = .9609, F1 = .9226, AUC = .9897

xgboost = create_model('xgboost') #Accuracy = .9612, F1 = .9234, AUC = .9885

catboost = create_model('catboost') #Accuracy = .9609, F1 = .9233, AUC = .9892

rf = create_model('rf') #Accuracy = .9626, F1 = .9258, AUC = .9904

lgbm = create_model('lightgbm') #Accuracy = .9593, F1 = .9202, AUC = .9895

# tuned_models were only included when performance is improved

#print(etrees)
finalize_model(etrees)
save_model(etrees, 'xant_etrees')

#print(tuned_xgboost)
finalize_model(xgboost)
save_model(tuned_xgboost, 'xant_xgb')

#print(catboost)
finalize_model(catboost)
save_model(catboost, 'xant_cboost')

#print(tuned_rf)
finalize_model(rf)
save_model(tuned_rf, 'xant_rf')

#print(lgbm)
finalize_model(lgbm)
save_model(lgbm, 'xant_lgbm')

blender_specific = blend_models(estimator_list = [rf, lgbm],
                method = 'soft')

# print(blender_specific)
finalize_model(blender_specific)
save_model(blender_specific, 'xant_blended')


print('PyCaret training ended \n\n')






