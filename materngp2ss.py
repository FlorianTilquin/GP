#! /usr/bin/env python
# -*- encoding:utf-8 -*-
#
# MATERNGP2SS State space representation of GP with Matern kernel
#
# Syntax:
#   [F,L,H,q] = materngp2ss(phi,sigma,p)
#
# In:
#     phi - Precision parameter (range)
#   sigma - Amount of variation (partial sill)
#       p - Integer part of the smoothness parameter nu = p + 1/2
# Out:
#   F - Feedback matrix
#   L - Noise gain
#   H - Measurement matrix
#   q - Spectral density of noise
#   P_inf - Stationary covariance
#
# Description:
#
#   Form a state space representation to Gaussian process with
#   Matern covariance function
#
#     C(t,t') = sigma**2 exp(-sqrt(2*(p+1)/2) tau/l) Gamma(p+1)/Gamma(2*p+1)
#             * sum_{i=0}**p [(p+i)!/(i!(p-i)!)(sqrt(8(p+1/2))tau/l)**(p-i)]
#
#   where tau = |t-t'| and l is the length-scale (phi = sqrt(2)/2*l).
#

# Copyright (C) 2010 Jouni Hartikainen, Simo Särkkä
#
# This software is distributed under the GNU General Public
# Licence (version 2 or later) please refer to the file
# Licence.txt, included with the software, for details.
from scipy import special
from scipy.linalg import *
from numpy import *

def care(a, b, q, r):
    try:
        g = inv(r)
    except LinAlgError:
        raise ValueError('Matrix R in the algebraic Riccati equation solver is ill-conditioned')
    g = dot(dot(b, g), b.conj().transpose())

    z11 = a
    z12 = -1.0*g
    z21 = -1.0*q
    z22 = -1.0*a.conj().transpose()

    z = vstack((hstack((z11, z12)), hstack((z21, z22))))

    # Note: we need to sort the upper left of s to have negative real parts,
    #       while the lower right is positive real components (Laub, p. 7)
    [s, u, sorted] = schur(z, sort='lhp')

    (m, n) = u.shape

    u11 = u[0:m//2, 0:n//2]
    u21 = u[m//2:m, 0:n//2]
    u11i = inv(u11)

    return dot(u21, u11i)


def pascal(n):
    triangle = ones((n,n))
    for i in range(1,n):
        for j in range(1,n):
            triangle[i,j] = triangle[i,j-1]+triangle[i-1,j]
    return triangle

def materngp2ss(phi,sigma,p):
    #scale = 0
    #if sigma < 0.001:
    #    sigma *= 1000
    #    scale = 1
    la = 2*sqrt(p+0.5)*phi
    c = (sigma**2)*2*pi**(0.5)*exp(special.gammaln(p+1)-special.gammaln(p+0.5))*la**(2*p+1)
    q = c

    ppoly = pascal(p+2)
    ppoly = diag(ppoly[::-1,:]).T
    ppoly = ppoly[1:]

    lav = la**array(range(1,p+2))
    F = diag(ones(p),1)
    lp = -lav*ppoly
    F[p,:] = lp[::-1]

    L = zeros((p+1,1))
    L[p] = 1
    H = zeros((1,1+p))
    H[0,0] = 1
    (n,p) = shape(F)
    P_inf = care(F.T,zeros((n,p)),dot(L,q*L.T),eye(p))
    #if scale:
    #    q /= 10**6
    #    P_inf /= 10**6
    return F,L,H,q,P_inf
