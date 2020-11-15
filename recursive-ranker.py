#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 11:26:05 2020

@author: danielfurman
"""

import pandas as pd
import numpy as np
from scipy.stats import spearmanr
import sys


def recursive_ranker(
        covariance,
        feature_importance,
        threshold,
        raw_data):

    '''
    
    This function recursively removes the most correlated modeling variables,
    such that the final set is below a Spearman's threshold, using the rank of
    feature importance scores to decide between variables to drop. 

    covariance: Pandas object containing the covariance matrix, with
        correlations between modeling variables, by definition containing
        ones along the diagonal. Variable names should only be above the
        entries.

    feature_importance: Pandas object containing a model's feature importance
        scores in the first row, with the same order of variables as the
        covariance matrix. Feature importance is generally defined as
        techniques that assign a score to input features based on how useful
        they are at predicting a target variable during classification.

    threshold: A correlation value for which features are filtered below,
        Thresholds between 0.5 - 0.7 are commonly used (e.g. Dormann et al.,
        2013, doi: 10.1111/j.1600-0587.2012.07348.x).

    raw_data: The database from which the original covariance matrix was
        created from.

    Warning:
    --------
    * The Pandas dataframes should have the same variables order.
    * Make sure dependencies are installed: pandas, np, scipy.stats.spearmanr
    
    Example:
    --------
    * See https://nbviewer.jupyter.org/github/daniel-furman/ensemble-climate-
    projections/blob/main/Comparing_MLs.ipynb, Appendix 1

    '''
    
    # initial transformations
    covariancenp = pd.DataFrame.to_numpy(covariance)
    covariancenp = np.triu(covariancenp)    
    covariance = pd.DataFrame(covariancenp, columns=list(covariance))
    for i in np.arange(0,len(covariance)):
        covariance.rename(index={i: list(covariance)[i]}, inplace=True)
    covariance = covariance.abs()
    covar = covariance.copy(deep=True)
    for i in np.arange(0, len(covar)):
        for p in np.arange(0, len(covar)):
            if covar.iloc[i, p] < threshold:
                covar.iloc[i, p] = np.NaN
        covar.iloc[i, i] = np.NaN
    covariance_bool = covar.isna()
    # check order of variables
    if list(covariance) != list(feature_importance):
        sys.exit('Variable names need to be consistent')

    # stopping case
    if covariance_bool.all(axis=None):
        fin = list(covariance)
        print('\nfinal set of variables: ', fin)
        print('\nCovariance matrix (r < ', str(threshold), '):\n')
        print(covariance)

    # recursion call
    else:
        for i in np.arange(0, len(covariance)):
            for p in np.arange(0, len(covariance)):
                if covariance.iloc[i, p] < threshold:
                    covariance.iloc[i, p] = np.NaN
            covariance.iloc[i, i] = np.NaN
        maximum_corr = np.max(np.max(covariance))
        for i in np.arange(0, len(covariance)):
            for p in np.arange(0, len(covariance)):
                if covariance.iloc[i, p] == maximum_corr:
                    colname = list(covariance)[p]
                    rowname = list(covariance)[i]
        min_imp = np.min([feature_importance[colname],
                          feature_importance[rowname]])
        if feature_importance[colname].loc['importance'] == min_imp:
            raw_data.drop([colname], axis=1, inplace=True)
            feature_importance.drop([colname], axis=1, inplace=True)
            print('Comparing', rowname, 'or', colname, ' | Dropping', colname)
        else:
            raw_data.drop([rowname], axis=1, inplace=True)
            feature_importance.drop([rowname], axis=1, inplace=True)
            print('Comparing', rowname, 'or', colname, ' | Dropping', rowname)
        covariance = pd.DataFrame(spearmanr(raw_data).correlation,
                                  columns=list(feature_importance))
        # recursion call
        recursive_ranker(covariance, feature_importance,
                         threshold, raw_data)
