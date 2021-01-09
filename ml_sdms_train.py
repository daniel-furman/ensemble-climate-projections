#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 11:26:03 2020

@author: danielfurman

"""

# Binary classification with ten BioClim features. Five are a function of
# precipitation and five of temperature, all de-correlated below a 0.5
# correlation threshold. We use PyCaret to train and tune (10-fold cv) our
# models from a  train set that contains 80%  of the total data. We find that
# random forest, xgboost, lgbm, catboost, and extra trees algorithms perform
# best. The last step was to create blended model from these five.

from pycaret.classification import setup
from pycaret.classification import create_model
from pycaret.classification import finalize_model
from pycaret.classification import save_model
from pycaret.classification import blend_models
from pycaret.classification import compare_models
from pandas import read_csv

#data = read_csv('data/envtrain_xv.csv')
data = read_csv('data_2.0/envtrain_xv.csv')

#data = data.drop(['Unnamed: 0'], axis=1)
exp_clf = setup(data, target='pa', log_experiment = True,
                experiment_name = 'xv-21', session_id = 110,
                numeric_features = ['bclim14'])

# create models
etrees = create_model('et')
xgboost = create_model('xgboost')
catboost = create_model('catboost')
rf = create_model('rf')
lgbm = create_model('lightgbm')
log = create_model('lr')

# save models as .pkl files
finalize_model(etrees)
save_model(etrees, 'classifier_models(pkl)/xant_etrees')

finalize_model(xgboost)
save_model(xgboost, 'classifier_models(pkl)/xant_xgb')

finalize_model(catboost)
save_model(catboost, 'classifier_models(pkl)/xant_cboost')

finalize_model(rf)
save_model(rf, 'classifier_models(pkl)/xant_rf')

finalize_model(lgbm)
save_model(lgbm, 'classifier_models(pkl)/xant_lgbm')

finalize_model(log)
save_model(log, 'classifier_models(pkl)/xant_log')

blender_specific = blend_models(estimator_list=[
    etrees, lgbm, rf], method='soft')
#blender_specific = blend_models(estimator_list=[
    #etrees, lgbm, catboost], method='soft')


finalize_model(blender_specific)
save_model(blender_specific, 'classifier_models(pkl)/xant_blended')

print('PyCaret training ended \n\n')

compare_models()  # print ordered 10-fold cv scores
