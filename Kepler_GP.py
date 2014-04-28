#! /usr/bin/env python
# -*- encoding: utf-8 -*-

# Filename: Kepler_GP
# GP fitting with state-based model approach to Kepler Data

from Kepler_data import *
from matplotlib import pyplot as plt
from lti_disc import *
from materngp2ss import *
from KF import *

def Kepler_GP(Kid , Beg , End , phi , sigma, Z , p = 6, plot = 0, plottitle = 0):  
    
    # Data recuperation
    T,Y,C = reg(Kid, Beg, End)
    Y = Y - 1
    n = len(Y)
    
    # Kalman Filter application
    Xpre, Ppre, ELL = KF(Y,phi,sigma,Z,p)

    if plot == 0:
        return ELL
    if plot:
        print ELL     # For information
        Var = sqrt(abs(Ppre))
        plt.figure()
        plt.plot(T , Y , '.r')
        plt.plot(T , Xpre , '-k')
        #plt.plot(T , Xpre + Var, '--b')
        #plt.plot(T , Xpre - Var, '--b')
        if plottitle != 0:
            plt.savefig(plottitle+'.png')
        else:
            plt.show()
