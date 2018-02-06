#!/usr/bin/python
"""
Using optimal estimation for estimating the best physical parameter values 
based on measured brightness temperatures. 
Iterations start from a set of input parameters, also containing our priors

Usage:
 >>>   optimize(<Tb>, <inputpar>)
"""

import numpy as np
import amsrrtm2
import icemodel

def optimize(Tb,inputpar):
    Se = inputpar["Se"] # Covariance matrix for measured brightness temperatures
    Sp = inputpar["Sp"] # Covariance matrix for parameters 
    P0 = inputpar["P"] # Initial values (for aabentvandsatmosfaeren?), also our priors 
    lowerlimits = inputpar["lowerlimits"] # Lower limits for physical parameters
    upperlimits = inputpar["upperlimits"] # Upper limits for physical parameters
    
    # Compute Ts_amsr (effective temperature) at the following frequencies: [6,10,19,23,37,50,52,89] 
    Ts_amsr = icemodel.Ts_amsr(Tb) # Based om e.g. Lise's equation for snow depth
    
    # Initialize adjoint matrix, M:
    len_p = len(P0)
    len_tb = len(Tb)
    M = np.matrix(np.ones((len_tb,len_p)))
 
    # Start with initial set of parameters: 
    P = np.copy(P0)
    
    # Iterations:
    max_iter = 6 # Maximum number of iterations
    eps = -1000 # 0.5 # Convergence criterion, if negative (&large): Iterate until convergence
    cost = np.inf
    dcost = np.inf
    
    n=0
    while (n<=max_iter and dcost>eps):
        n=n+1

        # Apply forward model using current parameters (P):         
        # To be called with the following input: 
        # amsrrtm2.amsr(V, W, L, Ti, Ts, Ts_amsr, Ts_amsr, c_ice, ev_ice, eh_ice)
        #   We have: P = [Ts,V,W,L,c_ice,sd]
        #   Assume: Ti = Ts = P[0]
        #   Ts_amsr is given for following frequencies: [6,10,19,23,37,50,"52",89] 
        #   ev_ice and eh_ice are given at following frequencies: [6,10,18,18,37,37,37,89] (approximating those in Ts_amsr)
        Ta=amsrrtm2.amsr(P[1],P[2],P[3],P[0],P[0],Ts_amsr,Ts_amsr,P[4],icemodel.ev_ice(Tb,P[0],P[5]),icemodel.eh_ice(Tb,P[0],P[5]))
        # Output brightness temperatures are given for the following frequencies:
        #   [Tv6,Th6,Tv10,Th10,Tv18,Th18,Tv23,Th23,Tv36,Th36,Tv89,Th89]
      
        # Compute adjoint matrix (M): 
        for i in range(len_p):
            # Compute partial derivates (parameters are pertubed by 1% of current value):
            Pa = np.copy(P); Pa[i]=P[i]+0.01*P[i]
            #Pb = np.copy(P); Pb[i]=P[i]-0.01*P[i]
            Pb = np.copy(P); #Pb[i]=P[i]-0.01*P[i]
            dp = 0.01*P[i]

            # Partial differences:
            dTdp =(amsrrtm2.amsr(Pa[1],Pa[2],Pa[3],Pa[0],Pa[0],Ts_amsr,Ts_amsr,Pa[4],icemodel.ev_ice(Tb,Pa[0],Pa[5]),icemodel.eh_ice(Tb,Pa[0],Pa[5])) - 
                   amsrrtm2.amsr(Pb[1],Pb[2],Pb[3],Pb[0],Pb[0],Ts_amsr,Ts_amsr,Pb[4],icemodel.ev_ice(Tb,Pb[0],Pb[5]),icemodel.eh_ice(Tb,Pb[0],Pb[5])))/dp
                 
            M[:,i]=dTdp.reshape(len_tb,1)
        
        # Compute S:
        S = (Sp.I + M.T*Se.I*M).I
        
        # New estimate for P (eq. 5.9, page 85 in Rodgers)
        P1 = P0.reshape((len_p,1)) + S*M.T*Se.I*(Tb.reshape((len_tb,1))-Ta.reshape((len_tb,1)) + 
                        M*(P.reshape((len_p,1))-P0.reshape((len_p,1))))
        P1 = np.squeeze(np.asarray(P1))
        
        # Ensure parameters stay within physical limits:
        #dlim = [0.01, 0.01, 0.01, 0.001, 0.01, 0.01] # why not set as max value? - to avoid 0'es
        dlim = [0.1, 0.1, 0.1, 0.01, 0.1, 0.1] # why not set as max value? - to avoid 0'es
        for k in range(len_p):
            if (P1[k] < lowerlimits[k]): P1[k] = lowerlimits[k]+dlim[k]
            if (P1[k] > upperlimits[k]): P1[k] = upperlimits[k]
            
        # Compute new cost:
        cost1=np.sqrt(sum((Tb.reshape((len_tb,1))-Ta.reshape(len_tb,1))**2))
        dcost = cost-cost1

        # Check for convergence
        if n==max_iter or dcost<eps:
            P_final = np.copy(P1) # Use as final value
            S_diag = np.diag(S)
            cost=np.sqrt(sum((Tb.reshape((len_tb,1))-Ta.reshape((len_tb,1)))**2))
            
        else:
            # Set as new value of parameters:
            P = np.copy(P1)
            # New cost value:
            cost = np.copy(cost1)
            
    # Convert ice concentration back to original, physical, values:
    P_final[4]=P_final[4]-1.0
    return P_final, S_diag, cost