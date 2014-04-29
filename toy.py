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
l = 20
sigma1 = 2
sigma2 = 0.5

# Output Parameters
phi = 1/(l*sqrt(2.))
sigma = sigma1
Z = sigma2

## Data points
# Mean Vector and Covariance matrix
mean = zeros(n)
cov = SE(u,u,[sigma1,l,sigma2],1)

# Drawns from the covariance matrix
X = rnd.multivariate_normal(mean, cov, 1)
X = reshape(X,(n,1))
plt.figure()
plt.plot(X,'.k')

## Prediction inputs test
Xp,a,b = KF(X,phi,sigma1,sigma2)


plt.plot(Xp,'-b')
plt.show()
