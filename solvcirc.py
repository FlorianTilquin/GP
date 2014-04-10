#!/usr/bin/env python
# -*- encoding: utf-8 -*-
# Filename: solvcirc.py

from numpy import *
from numpy.fft import *
from math import *
def solvcirc(c,b):
    # c Column of the circulant matrix
    # b Vector
    c = array(c)
    b = array(b)

    if len(shape(b)) == 1:
        b = reshape(b,(len(b),1))
    if len(shape(c)) == 1:
        c = reshape(c,(len(c),1))
    if c.shape[0] > c.shape[1]:
        c = c.T
    if b.shape[0] > b.shape[1]:
        b = b.T
    n = c.shape[0]
    s = b.shape[0]

    if s != n:
        print 'Dimensions are inconsitent'
        return('pbl dim')
    else:
        # Construction de la premiere colone de la matrice circulante
        # de maniere a ce qu'elle soit de taille 2^i
        i = floor(log(n)/log(2))
        if (i == log(n)/log(2)):
            C = c
        else:
            C = hstack((c,zeros((1,2**(i+1)-n))))
            b = hstack((b,zeros((1,2**(i+1)-n))))
    # Calcul de la solution
    P = fft(b)/fft(C)
    Y = ifft(P).real
    Y = Y[0:n].T
    return Y
