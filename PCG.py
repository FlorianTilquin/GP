#! /usr/bin/env python
# -*- encoding: utf-8 -*-
# Filename: PCG.py

from numpy import *
from numpy.fft import *
from math import *
from pdtMatToepVect import pdtMatToepVect
from scipy.linalg import *
from solvcirc import *

def PCG(col,row,b,tol = 10**-3):
# solve Ax = b with A = toeplitz(col,row)    
# with preconditionnate conjuguate gradient with a circulante matrix precontionner (Strang preconditionner)
# Considerations on the type of vectors we are dealing with. All vectors should have a vertical shape (n,1) (the 1 is important)   
    row,col,b  = array(row),array(col),array(b)
    
    if len(shape(b)) == 1:
        b = reshape(b,(len(b),1))
    if b.shape[0] < b.shape[1]:
        b = b.T
    if len(shape(col)) == 1:
        col = reshape(col,(len(col),1))
    if col.shape[0] < col.shape[1]:
        col = col.T
    if len(shape(row)) == 1:
        row = reshape(row,(len(row),1))
    if row.shape[0] < row.shape[1]:
        row = row.T
#Thus, this is the true length of our vectors
    m = col.shape[0]  
    n = row.shape[0]
    s = b.shape[0]
    if m != n or m != s or n != s:
        print 'incompatible dimensions for col, row, b'
        return 0
        
# Proper algorithm    
    # Computation of the precondition
    tolb = tol*norm(b)
    M = int(n/2)
    N = n-1-M
    Pc = row[:N+1]
    Pr = col[M:0:-1]
    Pf = vstack((Pr,Pc))
    Pf = reshape(Pf,(n,1))
    # Set up for the method
    step = 0
    x = zeros((n,1))
    r = b # - A*x
    z = solvcirc(Pf,r)
    p = z
    rho1 = float(dot(r.T,z).real)
    for i in range(0,n):    
        rho  = rho1
        d    = pdtMatToepVect(col,row,p)
        if float(dot(p.T,d).real) == 0:
            return 0,step
        a    = rho/float(dot(p.T,d).real)
        x   += a*p
        r   -= a*d
        normr = norm(r)
        if normr < tolb:
            break
        z    = solvcirc(Pf,r)
        rho1 = float(dot(z.T,r).real)
        beta = rho1/rho
        p    = z + beta*p
        if norm(pdtMatToepVect(col,row,x)-b) < tolb:
            break
        step += 1    
    return x,step
