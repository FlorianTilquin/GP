#! /usr/bin/env python
# -*- encoding: utf-8 -*-
# Filename:BAS.py

#Sandbox
from time import *
import matplotlib.pyplot as plt
from numpy import *
import scipy.optimize as opt
from Toep import *
from KF import *
from Kepler_data import *
from Kepler_GP import *

n = 100
p = 10
u = range(0,n)
A = 10   # f(t) = A*cos(w*t)
w = 2*pi/n
sigma = 0.1

# True distribution 
X = range(0,n)
X = reshape(array(X),(n,1))
Y = A * cos(w*X)

# Noised distribution
err = sigma*random.randn(n)
err = reshape(err,(n,1))
Yerr = Y + A*err

# Prediction inputs test
Xpre = arange(0,n,0.1)
Xpre = reshape(Xpre,(10*n,1))

# Kalman Filter
def LLH(theta):
    Ypre, ELL = KF(Yerr,theta[0],theta[1],theta[2])
    return ELL
L = 22
phi = 1./(sqrt(2)*L)
Sigma = 0.0015
Z = 0.001

Kepler_GP(1295333,4,7,phi,Sigma,Z,6,1)

#print ELL1,ELL2

#mini = opt.fmin(LLH,[0.01,0.1,0.25])
#print mini
#Ypre,ELL = KF(Yerr,mini[0],mini[1],mini[2])
#Ypre,ELL = KF(Yerr,0.3,20,2)

# Plotting the model
#plt.figure()
#plt.plot(Yerr,'+b')
#plt.plot(Ypre,'-k')
#plt.show()
