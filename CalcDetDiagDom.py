#! /usr/bin/env python
# -*- encoding: utf-8 -*-
# Filename: CalcFilterDiagDom.py

# Computes the determinant of a matrix in which most of the weight is concentrated in the near-diagonal coefficients
from numpy import *
from scipy.linalg import *
from Toep import SE

# Computes the  base changing matrix T and its inverse U, used later to compute A' = T*A*U, which is even more diagonaly centered
def ChgBase(n):
    T = zeros((n,n))
    U = zeros((n,n)) # U = T^-1
    if n % 2 == 0:
        p = n/2
        for k in range(0,2*p): #raws of T
            for l in range(0,p): #columns of T
                if k == 2*l+1 or k == 2*l:
                    T[k,l] = 1.
            for l in range(p,2*p):
                if k == 2 * (l - p):
                    T[k,l] = 1./2.
                if k == 2 * (l - p)+1:
                    T[k,l] = -1./2.
        for kk in range(0,2*p): #columns
            for ll in range(0,p): #raws
                if kk == 2*ll or kk == 2*ll+1:
                    U[ll,kk]  = 1./2.
            for ll in range(p,2*p):
                if kk == 2*(ll-p):
                    U[ll,kk] = 1.
                if kk == 2*(ll-p)+1:
                    U[ll,kk] = -1.
    return T,U                
    
    if n % 2 == 1:
        T,U = ChgBase(n+1)
        T = T[0:n,0:n] # This is one of the properties of the base changing matrix
        U = U[0:n,0:n]
        U[(n-1)/2,n-1] = 1
        return T,U

def Filter(A,minsize):
    m,n = shape(A)
    d   = 1
    k   = 0
    if m != n:
        print 'You must use a square matrix'
        return 0
        T,U = ChgBase( m ) #We compute the Haar base changing matrices
        A   = dot(dot(T,A),U) 
    if m > minsize:
        if m % 2 == 0: # We have to make a disjonction on the matrix size
            p = m/2
            B = A[0:p,0:p]            
            C = A[p:m,p:m] 
            D = Filter(B,minsize),Filter(C,minsize) # Apply recursively the algorithm to the sub-matrices of A we're interested in:UL and DR. 
            return D
        if m % 2 == 1: # Could use else: , just for clarity
            A1 = eye(m+1) # If the matrix size isn't even, we construct a bigger matrix by adding 1 row/col with only a 1 in the bottom right position, thus not changing the determinant value, by product
            A1[0:m,0:m] = A
            A = A1
            p = (m-1)/2
            B = A[0:p+1,0:p+1]            
            C = A[p+1:m+1,p+1:m+1] 
            D = Filter(B,minsize),Filter(C,minsize) # Apply recursively the algorithm to the sub-matrices of A we're interested in:UL and DR. 
    else:
        D = A[0:m,0:m]
    return D

# Computes the determinant of a diagonaly dominant matrix using the recursive algorithm above : computes each det of the minimum size matrix and multiplies them together as if the big matrix was diagonal by blocks.
def Det(A,p):
    D  = Filter(A,p)
    D  = array(D)
    dim = prod(D.shape[:-2])
    dim = (dim,D.shape[-2],D.shape[-1])
    D = reshape(D,dim)
    d  = 1
    for k in range(0,len(D)):
        d *= det(D[k])
    return d
