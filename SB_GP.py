#! /usr/bin/env python
# -*- encoding: utf-8 -*-

# Filename: Kepler_GP
# GP fitting with state-based model approach to Kepler Data

from matplotlib import pyplot as plt
from lti_disc import *
from materngp2ss import *

def SB_GP(Y , phi, sigma, Z , p = 7, alpha = 1, plot = 0, plottitle = 0):  
    
    # Sampling
    u = range(0,len(Y),alpha)

    # Settings for the method
    H = eye(1,p+1)
    I = eye(p+1)

    # Kalman Filter Matrices
    [F , L , tilde , q , p0] = materngp2ss(phi , sigma , p)
    [A , Q] = lti_disc(F , L , q)

    # Initialisation 
    Xtt0 = zeros((p+1,1))
    Ptt0 = p0
    Xpre = []
    Ppre = p0
    ELL  = log(2*pi)*0.5
    for k in u:
        Xtt0 = dot(A,Xtt0)
        Ptt0 = dot(A,dot(Ptt0,A.T)) + Q
        ELL += log(Ptt0[0,0]+Z)/2 + ((Y[k]-Xtt0[0])**2)/(2*(Ptt0[0,0]+Z))
        Kpre = Ptt0[:,0]/(Ptt0[0,0]+Z)
        Kpre = reshape(Kpre,(p+1,1))
        Xtt0 = Xtt0 + Kpre*(Y[k]-Xtt0[0])
        Ptt0 = dot(I-hstack((Kpre,zeros((p+1,p)))),Ptt0)
        if plot:
            Xpre = hstack((Xpre , Xtt0[0]))
            Ppre = hstack((Ppre , Ptt0))
            #Ptt0 = dot(I-dot(Kpre,H),Ptt0)
            #Xtt0 = Xtt0 + dot(Kpre,Y[k]-dot(H,Xtt0))
            #Kpre = dot(Ptt0,H.T)/(dot(dot(H,Ptt0),H.T)+Z)
    if plot == 0:
        return ELL
    if plot:
        V = Ppre[0,range(0,Ppre.shape[1],p+1)]
        Var = sqrt(V[1:])
        plt.figure()
        plt.plot(Y , '+b')
        plt.plot(Xpre , '-k')
        plt.plot(Xpre + Var, '--b')
        plt.plot(Xpre - Var, '--b')
        if plottitle != 0:
            plt.savefig(plottitle+'.png')
        else:
            plt.show()
