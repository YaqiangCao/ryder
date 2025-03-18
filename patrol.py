#!/usr/bin/env python
#--coding:utf-8--
"""
PATROL algorithm for cross-sample normalized ChIP-seq/DNase-seq/ATAC-seq differential peaks detection. 
"""

import numpy as np
import pandas as pd
from scipy import linalg
from scipy.stats import chi2

def mahalanobis(mat):
    """
    Caculate the mahalanobis distance.

    according to: https://www.statology.org/mahalanobis-distance-python/

    @param mat: np.array, row is each item and column is each variable 

    @return dis: np.array,Mahalanobis distance
    @return ps: np.array, chi-square test p-values
    """
    #covariance matrix 
    cov = np.cov(mat, rowvar=False)
    #inverse covariance matrix 
    invCov = np.linalg.inv(cov)
    #center 
    center = np.mean(mat, axis=0)
    #mahalanobis distance
    mu = mat - center
    dis = np.dot( np.dot(mu,invCov), mu.T ).diagonal()
    #Chi-square test p-values for detecting outliers 
    ps = 1 - chi2.cdf(dis, mat.shape[1]-1)
    return dis, ps


def twoPassesMDTest(data, pcut=0.01):
    """
    Perform MD test with two passes. First pass use all data, second pass using cov and center from data without outliers. 
    @param data: pd.DataFrame, row is item and column is sample.
    @param pcut: Chi-Square p-value cutoffs
    """
    #first pass test
    #mahalanobis distance and pvaues
    dis, ps = mahalanobis(data.values)
    dis = pd.Series(dis, index=data.index)
    ps = pd.Series(ps, index=data.index)
    inds = ps[ps < pcut].index

    #second pass test, only using data without outliers to caculate cov matrix and center
    ndata = data.drop(inds)
    #covariance matrix
    cov = np.cov(ndata.values, rowvar=False)
    #inverse covariance matrix
    invCov = np.linalg.inv(cov)
    #center
    center = np.mean(ndata.values, axis=0)

    #mahalanobis distance for all data
    mu = data.values - center
    dis = np.dot(np.dot(mu, invCov), mu.T).diagonal()
    #Chi-square test p-values for detecting outliers from all data
    ps = 1 - chi2.cdf(dis, data.shape[1] - 1)

    dis = pd.Series(dis, index=data.index)
    ps = pd.Series(ps, index=data.index)
    return dis, ps


