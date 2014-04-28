#! /usr/bin/env python
# -*- encoding: utf-8 -*-

# Filename: Kepler_GP
# GP fitting with state-based model approach to Kepler Data

from lti_disc import *
from materngp2ss import *

def KF(Y, phi, sigma, Z,llh = 0, p = 7):
    n = len(Y)

    # Use the Särkkä routines to compute the transition matrices for the filter
    F, L, H, Qc, P_inf = materngp2ss(phi, sigma, p)
    G, Q = lti_disc(F, L,Qc)
    H = H.T  # H should be a line vector
    
    # Initialisation
    X = Y[0]*H.T
    P = P_inf
    Xpre = []
    Ppre = []
    NLL  =  0.5 * log(2*pi)
    
    for k in range(0,len(Y)):
        X = dot(G , X)
        P = dot(G , dot(P , G.T)) + Q
        NLL += 0.5*log(abs(P[0,0]+Z))+0.5*((Y[k]-X[0])**2)/abs(P[0,0]+Z) #### /¡!¡\ Ici le abs n'est pas naturel !! Il ya  un problème !!
        
        # Kalman gain
        K = reshape(P[:,0]/(P[0,0]+Z),(p+1,1))
            
        # Predictions
        Xp = X + K * ( Y[k]-X[0] ) 
        Pp = dot(eye(p+1)-dot(K , H) , P)
        Xpre.append( Xp[0] )
        Ppre.append( Pp )

        # Update
        X = Xp
        P = Pp
    # Reshaping: Predictions should be vectors, not lists    
    Xpre = reshape(Xpre , (n,1) )
    #Ppre = reshape(Ppre , (p+1 , (n+1)*(p+1) ) )
    #Ppre = reshape(Ppre[0,range(p+1, (n+1)*(p+1) , p+1)].T, (n,1) )
    if llh==1:
        return NLL

    return Xpre,Ppre,NLL
