#! :usr/bin/env python
# -*- encoding:utf-8 -*-

#LTI_DISC  Discretize LTI ODE with Gaussian Noise
#
# Syntax:
#   [A,Q] = lti_disc(F,L,Qc,dt)
#
# In:
#   F  - NxN Feedback matrix
#   L  - NxL Noise effect matrix        (optional, default identity)
#   Qc - LxL Diagonal Spectral Density  (optional, default zeros)
#   dt - Time Step                      (optional, default 1)
#
# Out:
#   A - Transition matrix
#   Q - Discrete Process Covariance
#
# Description:
#   Discretize LTI ODE with Gaussian Noise. The original
#   ODE model is in form
#
#     dx/dt = F x + L w,  w ~ N(0,Qc)
#
#   Result of discretization is the model
#
#     x[k] = A x[k-1] + q, q ~ N(0,Q)
#
#   Which can be used for integrating the model
#   exactly over time steps, which are multiples
#   of dt.

# History:
#   11.01.2003  Covariance propagation by matrix fractions
#   20.11.2002  The first official version.
#
# Copyright (C) 2002, 2003 Simo Särkkä
#
# $Id: lti_disc.m 111 2007-09-04 12:09:23Z ssarkka $
#
# This software is distributed under the GNU General Public 
# Licence (version 2 or later) please refer to the file 
# Licence.txt, included with the software, for details.
from numpy import *
from scipy.linalg import expm
def lti_disc(F, L =[] , Q =[], dt = 1):
    if L == []:
        L = eye(shape(F)[0])
    if Q == []:
        Q = zeros(shape(F))
       
    # Closed form integration of transition matrix
    A = expm(F*dt)
    # Closed form integration of covariance by matrix fraction decomposition
    QQ  = dot(L,dot(Q,L.T))
    n   = F.shape[0]
    h1 = hstack((F,QQ))
    h2 = hstack((zeros((n,n)),-F.T))
    Phi = vstack((h1,h2))
    AB  = dot(expm(Phi*dt),vstack((zeros((n,n)),eye(n))))
    Q   = AB[:n,:]/AB[n:(2*n),:]
    return A,Q 
