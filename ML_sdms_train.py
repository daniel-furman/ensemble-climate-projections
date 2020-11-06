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

data = read_csv('data/envtrain_xv.csv')
data = data.drop(['Unnamed: 0'], axis=1)
exp_clf = setup(data, target='pa')

# create models
etrees = create_model('et')
xgboost = create_model('xgboost')
catboost = create_model('catboost')
rf = create_model('rf')
lgbm = create_model('lightgbm')
log = create_model('lr')

# save models as .pkl files
finalize_model(etrees)
save_model(etrees, 'xant_etrees')

finalize_model(xgboost)
save_model(xgboost, 'xant_xgb')

finalize_model(catboost)
save_model(catboost, 'xant_cboost')

finalize_model(rf)
save_model(rf, 'xant_rf')

finalize_model(lgbm)
save_model(lgbm, 'xant_lgbm')

blender_specific = blend_models(estimator_list=[
    rf, etrees, xgboost, lgbm, catboost], method='soft')

finalize_model(log)
save_model(log, 'xant_log')

finalize_model(blender_specific)
save_model(blender_specific, 'xant_blended')

print('PyCaret training ended \n\n')

compare_models()  # print ordered 10-fold cv scores
