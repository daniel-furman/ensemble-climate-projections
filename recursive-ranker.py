#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 11:26:05 2020

@author: danielfurman
"""

import pandas as pd
import numpy as np
from scipy.stats import spearmanr


def recursive_ranker(
        covariance,
        feature_importance,
        covariance_bool,
        threshold,
        raw_data):

    '''
    
    This function recursively removes the most correlated modeling variables,
    such that the final set is below a Spearman's threshold, using the rank of
    feature importance scores.

    covariance: Pandas object containing the covariance matrix, with
        correlations between modeling variables, by definition containing
        ones along the diagonal. Variable names should be above the entries

    feature_importance: Pandas object containing model feature importance
        scores in the first row. Feature importance is generally defined as
        techniques that assign a score to input features based on how useful
        they are at predicting a target variable during classification.

    covariance_bool: Initialize by passing the covariance dataframe to
        covariance.isna().

    threshold: A correlation value for which features are filtered below,
        Thresholds between 0.5 - 0.7 are commonly used (e.g. Dormann et al.,
        2013, doi: 10.1111/j.1600-0587.2012.07348.x).

    raw_data: The database from which your original covariance matrix was
        created from.

    Warning
    --------
    * The Pandas dataframes should all have variable names in the same order.
    * Make sure dependencies are installed: pandas, np, scipy.stats.spearmanr
    
    Example
    --------
    * See https://nbviewer.jupyter.org/github/daniel-furman/ensemble-climate-
    projections/blob/main/Comparing_MLs.ipynb, Appendix 1

    '''

    if covariance_bool.all(axis=None):
        print(covariance, '\n')
        fin = list(covariance)
        print('final set of variables: ', fin)

    else:
        for i in np.arange(0, len(covariance)):
            for p in np.arange(0, len(covariance)):
                if covariance.iloc[i, p] < threshold or covariance.iloc[
                        i, p] == 1:
                    covariance.iloc[i, p] = np.NaN
        maximum_corr = np.max(np.max(covariance))

        for i in np.arange(0, len(covariance)):
            for p in np.arange(0, len(covariance)):
                if covariance.iloc[i, p] == maximum_corr:
                    colname = list(covariance)[p]
                    rowname = list(covariance)[i]
        min_imp = np.min([feature_importance.iloc[0, p],
                          feature_importance.iloc[0, p]])
        print('Dropping', rowname, 'or', colname, '\n')

        if feature_importance[colname].loc['importance'] == min_imp:
            raw_data.drop([colname], axis=1, inplace=True)
            feature_importance.drop([colname], axis=1, inplace=True)
        else:
            raw_data.drop([rowname], axis=1, inplace=True)
            feature_importance.drop([rowname], axis=1, inplace=True)
        covariance = pd.DataFrame(spearmanr(raw_data).correlation,
                                  columns=list(feature_importance))
        covariancenp = pd.DataFrame.to_numpy(covariance)
        covariancenp = np.triu(covariancenp)
        covariance = pd.DataFrame(covariancenp, columns=list(covariance))

        for i in np.arange(0, len(covariance)):
            covariance.rename(index={i: list(covariance)[i]}, inplace=True)
        covariance = covariance.abs()
        covar = covariance.copy(deep=True)

        for i in np.arange(0, len(covar)):
            for p in np.arange(0, len(covar)):
                if covar.iloc[i, p] < threshold or covar.iloc[i, p] == 1:
                    covar.iloc[i, p] = np.NaN
        covariance_bool = covar.isna()

        recursive_ranker(covariance, feature_importance, covariance_bool,
                         threshold, raw_data)
