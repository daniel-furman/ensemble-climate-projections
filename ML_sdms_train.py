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
exp_clf = setup(data, target = 'pa')

# create models
etrees = create_model('et')
xgboost = create_model('xgboost')
catboost = create_model('catboost') 
rf = create_model('rf') 
lgbm = create_model('lightgbm')

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

blender_specific = blend_models(estimator_list = [rf, lgbm],
                method = 'soft')
finalize_model(blender_specific)
save_model(blender_specific, 'xant_blended')

print('PyCaret training ended \n\n')

compare_models() # print ordered 10-fold cv scores
