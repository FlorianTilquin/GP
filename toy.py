#! /usr/bin/env python
# -*- encoding: utf-8 -*-
from numpy import *
from numpy import random as rnd

from scipy.linalg import *
import scipy.optimize as opt
from Toep import *
from KF import *
import matplotlib.pyplot as plt

n = 500
u = range(0,n)

# Input Parameters
l = 10
sigma1 = 2
sigma2 = 0.1

# Output Parameters
phi = 0.2
sigma = sigma1
Z = sigma2

## Data points
# Mean Vector and Covariance matrix
mean = zeros(n)
cov = SE(u,u,[sigma1,l,sigma2],1)

# Drawns from the covariance matrix
X = rnd.multivariate_normal(mean, cov, 1)
X = reshape(X,(n,1))
#plt.figure()
#plt.plot(X,'+b')

## Prediction inputs test
def NLL(x):
    try:
        NLL = KF(X,x[0],x[1],x[2],1)
    except:
        NLL = inf
    return NLL
#x = opt.fminbound(NLL,array([0.01,1,0.05]),array([0.5,5,0.5]))
#print x
N = []
for a in arange(0.01,0.1,0.01):
    for b in arange(1,5,0.1):
        for c in arange(0.001,0.1,0.001):
            N.append(NLL([a,b,c]))
print min(N)
#print NLL([0.01,4,0.085])
#print argmin(N),min(N)
#Xp,a,b = KF(X,x,sigma,Z)
#plt.plot(Xp,'-k')
#plt.show()
