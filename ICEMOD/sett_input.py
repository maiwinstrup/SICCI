#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Set priors for physical parameters: mean, covariance, limits, and set covariance 
for measured brightness temperatures.

The priors also functions as first guess for the iterations.

Usage:
 >>>   setinputpar(<Tb>)
"""
import numpy as np
import comisoiceconc
from icemodel import liseequation

def setinputpar(Tb):
    # Covariance matrix for measured brightness temperatures (Tb): 
    # A diagonal nxn matrix, here n = 12 (number of temperature channels)
    #v = (np.ones([12]))*0.5
    v = (np.ones([12]))*(1**2)
    v[10:12]=2**2 #1 # 89 GHz channels are allowed higher variance 
    # Previous values ??: 
    # v = [0.0900, 0.1089, 0.2209, 0.2916, 0.2304, 0.2116, 0.2025, 0.1936, 0.2025, 0.1600, 0.1600, 0.1600]
    
    ##kovariansmatricen for de estimerede straalingstemperatur vaerdier S
    #Se=np.matrix([[0.2, 0.15, 0.1, 0.05, 0, 0, 0, 0, 0, 0, 0, 0],\
                  #[0.15, 0.2, 0.15, 0.1, 0.05, 0, 0, 0, 0, 0, 0, 0],\
                  #[0.1, 0.15, 0.2, 0.15, 0.1, 0.05, 0, 0, 0, 0, 0, 0],\
                  #[0.05, 0.1, 0.15, 0.2, 0.15, 0.1, 0.05, 0, 0, 0, 0, 0],\
                  #[0, 0.05, 0.1, 0.15, 0.2, 0.15, 0.1, 0.05, 0, 0, 0, 0],\
                  #[0, 0, 0.05, 0.1, 0.15, 0.2, 0.15, 0.1, 0.05, 0, 0, 0],\
                  #[0, 0, 0, 0.05, 0.1, 0.15, 0.2, 0.15, 0.1, 0.05, 0, 0],\
                  #[0, 0, 0, 0, 0.05, 0.1, 0.15, 0.2, 0.15, 0.1, 0.05, 0],\
                  #[0, 0, 0, 0, 0, 0.05, 0.1, 0.15, 0.2, 0.15, 0.1, 0.05],\
                  #[0, 0, 0, 0, 0, 0, 0.05, 0.1, 0.15, 0.2, 0.15, 0.1],\
                  #[0, 0, 0, 0, 0, 0, 0, 0.05, 0.1, 0.15, 0.2, 0.15],\
                  #[0, 0, 0, 0, 0, 0, 0, 0, 0.05, 0.1, 0.15, 0.2]])
    
    # Covariance matrix:
    Se = np.matrix(np.diag(v))
    
    # Input values (first guess and priors) for physical parameters in 
    # radiative balance equation: 
    # Using initial vaerdier for aabentvandsatmosfaeren (???)
    Ts=275.0 # Default surface temperature (IST or SST); assumed equal to the 2m air temperature [K] 
    #Ts=279.5 #sst, overfladetemperatur, effektiv temperatur
    V=3.6 # Default total column water vapour [kg/m2]
    W=5.0 # Default wind speed [m/s]
    L=0.08 # Default total column cloud liquid water content [kg/m2]
    # Sea ice concentration: 
    cf=comisoiceconc.bootstrap_f(Tb[4],Tb[8])
    # c_ice: 1.0 is open water and c_ice=2.0 is sea ice
    c_ice=cf+1.0   #OBS: openwater eller hva? vi mangler en isemissivitetsfunktion
    # Snow depth:
    sd = liseequation(Tb)  # Linear regresion for snow depth [m]
    
    # Create parameter vector (1xm)
    P0 = np.array([Ts,V,W,L,c_ice,sd]) 

    # Set covariance matrix for physical parameters (P):
    # A diagonal mxm matrix (m=6).        
    w = [15, 5, 5, 0.3, 0.4, 0.1] # Ordered as above: Ts,V,W,L,c_ice,sd
    Sp = np.matrix(np.diag(w))
    
    ##kovariansmatricen for de estimerede fysiske vardier P
    #Sp=np.matrix([[0.5, 0.1, 0.1, 0.0],\
                  #[0.1, 0.5, 0.0, 0.0],\
                  #[0.0, 0.0, 15.0, 0.1],\
                  #[0.0, 0.1, 0.1, 0.5]])

    
    # Lower and upper limits for physical parameters
    lowerlimits=np.array([200.0, 0.0, 0.0, 0.0, 0.9, 0.0])
    upperlimits=np.array([280.0,40.0,50.0, 1.0, 2.1, 1.0])
    ##oevre og nedre graenser til de fysiske parametre
    #L=np.array([270.0,320.0,0.0,40.0,0.0,50.0,0.0,1.0])
    ##oevre og nedre graenser til de fysiske parametre
    #L=np.array([0.0,0.55,0.0,1.0,190.0,273.15,0.0,1.0])
    
    # Create dictionary with input parameters
    inputpar = {'Sp':Sp, 'Se':Se, 'P':P0, 'lowerlimits':lowerlimits, 'upperlimits':upperlimits}
    
    return inputpar