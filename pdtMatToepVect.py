#!/usr/bin/env python
# -*- encoding: utf-8 -*-
# Filename: pdtMatToepVect.py

from numpy import *
from numpy.fft import *
from math import *
def pdtMatToepVect(c,r,b):
    # c,r column et r of the toeplitz matrib
    # b Vector
    c = array(c)
    r = array(r)
    b = array(b)

    if len(shape(b)) == 1:
        b = reshape(b,(len(b),1))
    if len(shape(r)) == 1:
        r = reshape(r,(len(r),1))
    if len(shape(c)) == 1:
        c = reshape(c,(len(c),1))
    if c.shape[0] < c.shape[1]:
        c = c.T
    if r.shape[0] < r.shape[1]:
        r = r.T
    if b.shape[0] < b.shape[1]:
        b = b.T
    n = c.shape[0]
    m = r.shape[0]
    s = b.shape[0]

    if n != m or m != s or s != n:
        print 'Must use a square matrix'
        return('pbl dim')
    else:
        # Construction de la premiere colone de la matrice circulante
        # de maniere a ce qu'elle soit de taille 2^i
        i = floor(log(2*n)/log(2))
        if (i == log(2*n)/log(2)):
            Ct1 = vstack((c,r[:0:-1]))
            p   = i
        else:
            Ct1 = vstack((c,zeros((2**(i+1)-2*n+1,1)),r[:0:-1]))
            p   = i+1
    # calcul des vp de Ct
    D = fft(Ct1.T)
    # Modification de b
    b = vstack((b,zeros((2**p-n,1))))
    # Calcul de la solution
    P = D*fft(b.T)
    Y = ifft(P).T.real
    Y = reshape(Y[0:n],(n,1))
    return Y
