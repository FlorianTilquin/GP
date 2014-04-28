#! /usr/bin/env python 
# -*- encoding: utf-8 -*-
# Filename: Toep.py

# Create a Toeplitz Covariance array from data Y,T,C

import numpy as np
from scipy.spatial.distance import *
from scipy.special import *
def eucl(X1,X2):
    X1,X2 = np.array(X1),np.array(X2)
    if len(np.shape(X1)) == 1:
        X1 = np.reshape(X1,(len(X1),1))
    if X1.shape[0] < X1.shape[1]:
        X1 = X1.T
    if len(np.shape(X2)) == 1:
        X2 = np.reshape(X2,(len(X2),1))
    if X2.shape[0] < X2.shape[1]:
        X2 = X2.T
    D = cdist(X1,X2,'sqeuclidean')
    return D

#########################
#Some Kernels definition#
#########################


#Square-exponential kernel (SE)
def SE(X1, X2, theta, white_noise = False):
    '''
    Squared exponential covariance function (SE)
    '''
    X1,X2 = np.array(X1),np.array(X2)
    s = np.shape(X1)
    n = s[0]
    D = eucl(X1,X2)
    K = theta[0]**2 * np.exp(- D / (2.*(theta[1]**2))) # calculate covariance array
    if white_noise == True: # add white noise
        K += (np.eye(n) * (theta[2]**2))
    K = np.array(K)
    return K

#Rational-Quadratic kernel (RQ)
def RQ(X1, X2, theta, white_noise = False):
    ''' 
    Rational quadratic covariance function (RQ)
    '''
    X1, X2 = np.array(X1), np.array(X2)
    D2     = eucl(X1,X2) # calculate squared Euclidean distance
    K      = (1 + D2 /(2*theta[0]*(theta[1]**2)))**(-theta[0]) #calculate covariance array
    if white_noise == True: # add white noise
        K += (np.identity(X1[:,0].size) * (theta[2]**2))
    K = np.array(K)
    return K

#periodic kernel (PER)
def PER(X1, X2, theta, white_noise = False):
    ''' 
    Periodic covariance function (PER)
    '''
    X1, X2 = np.array(X1), np.array(X2)
    s  = np.shape(X1)
    n  = s[0]
    D2 = cdist(X1, X2,'euclidean') # calculate squared Euclidean distance
    K  = theta[0]*np.exp((np.sin(-D2*np.pi/theta[1])**2)/(2*theta[2]**2)) #calculate covariance array
    if white_noise == True: # add white noise
        K += (np.eye(n) * (theta[3]**2))
    K = np.array(K)
    return K
#Matern kernel (Matern)
def Matern(X1, X2, m, theta, white_noise = False):
    '''
    Matern of order m covariance function
    '''
    X1,X2 = np.array(X1),np.array(X2)
    if len(np.shape(X1)) == 1:
        X1 = np.reshape(X1,(len(X1),1))
    if X1.shape[0] < X1.shape[1]:
        X1 = X1.T
    if len(np.shape(X2)) == 1:
        X2 = np.reshape(X2,(len(X2),1))
    if X2.shape[0] < X2.shape[1]:
        X2 = X2.T
    D = cdist(X1,X2,'euclidean')
    K = (theta[0]**2)/(2**(m-1)*gamma(m))*((np.sqrt(2*m)/theta[1]*D)**m)*kv(m,np.sqrt(2*m)/theta[1]*D)
    np.fill_diagonal(K,theta[0]**2)
    if white_noise:
        K += (np.identity(X1[:,0].size)*(theta[2]**2))
    K = np.array(K)
    return K
