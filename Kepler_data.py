#! /usr/bin/env python 
# -*- encoding: utf-8 -*-

# Filename: Kepler_data.py

# With inputs Kid of the star, begining quarter number beg, and finishing quarter (included) number end, returns flux, time and cadence for the star durind period of [Qbeg;Qend], and plot flux = f(time) if the option plot is added

import numpy as np
from matplotlib import pyplot as plt
from scipy.io.matlab import mio
from Kepler_rotate import *


# Give the position of the star in each quarter from the quarter Qbeg to the quarter Qend (included)
def read_Kepler(Kid,beg,end):
    # Find the position of the star in the first module of the quarter
    text = "../../Data/q"+str(beg)+"_objlist.txt"
 
    for line in open(text):
        if str(Kid) in line:
            a = int(float( line.split(  )[0] ))#mod
            b = int(float( line.split(  )[1] ))#out
            Pos = np.array([[a,b]])  #Pos = [mod , out] 
    for k in range(0,end-beg):
        Pos = np.vstack((Pos,Kepler_rotate(Pos[k][0],Pos[k][1])))
    return Pos

#Get the corresponding Kepler data and stack the raw data in vectors to be treated
def Kepler_data(Kid, beg, end, plot = 0):
    Pos  = read_Kepler(Kid, beg, end)
    Time = np.array([])
    Flux = np.array([])
    Cad  = np.array([])
    for k in range(0,end+1-beg):
        matflux = "../../Data/q"+str(beg+k)+"_mod"+str(Pos[k][0])+"_out"+str(Pos[k][1])+"_raw.mat" # Gets the corresponding data, provided they are in the right folder
        data = mio.loadmat(matflux)
        T    = data['time'].T
        Y    = data['flux_arr_pdc']
        K    = data['kid_arr']
        C    = data['cadence'].T
        n    = next(i for i in range(0,len(K)) if K[i] == Kid)
        Time = np.hstack((Time,T[0][2:len(T[0][:])]))
        Flux = np.hstack((Flux,Y[n][2:len(Y[n][:])]/np.median(Y[n][2:len(Y[n][:])])))
        Cad  = np.hstack((Cad, C[0][2:len(C[0][:])]))
    if plot:
        plt.plot(Time,Flux)
        plt.xlabel('Time')
        plt.ylabel('Flux')
        plt.title('Flux curve of star'+str(Kid)+'over the period [Q'+str(beg)+';Q'+str(end)+']')
        plt.show()
    return Time,Flux,Cad

def reg(Kid,Beg,End):
    T,Y,C = Kepler_data(Kid,Beg,End)
    l     = np.where(np.isnan(Y) == False) # Set of indexes for which both time and flux are finite, as NaN(time) C NaN(flux)
    C_reg = C[l]
    N     = C_reg.max()+1-C_reg.min()   # Number of shots that should have been taken
    C_reg = np.r_[C.min():C.max()+1]
    T_reg = T[l]
    delta = np.median(T[1:]-T[0:-1])  # Basicly, delta is the mean dt for which cadence * dt = laps time (where we have an output/input)
    T_reg = delta * C_reg +T[0]
    Y_reg = Y[l]
    T_reg = T_reg[l]
    C_reg = C_reg[l]
    return T_reg,Y_reg,C_reg

