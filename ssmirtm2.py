#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""---------------------------------------------------------------------------
wentz: a well calibrated ocean algorithm for ssm/i. JGR 102(C4), 8703-8718, 1997
computes the brightness temperature of the ocean surface at 19 22 and 37 GHz
Tb=f(V,W,L,Ts,Ti,c_ice)
V: columnar water vapor [mm]
W: windspeed over water [m/s]
L: columnar cloud liquid water [mm]
Ts: sea surface temperature [K]
Ti: ice surface temperature [K]
c_ice: ice concentration [0-1]
-------------------------------------------------------------------------------"""
import sys, os, string, math
import numpy as np

def wentz(V,W,L,Ts,Ti,c_ice,e_ice):
#def wentz():

    #V=0.0
    #W=0.0
    #L=0.0
    #Ts=275.0
    #Ti=275.0
    #c_ice=0.0
    #e_ice=[1.0,1.0,1.0,1.0,1.0,1.0]

    #table 1
    #19v 19h 22v 22h 37v 37h
    c0=[240.58E+0, 240.58E+0, 242.04E+0, 242.04E+0, 239.55E+0, 239.55E+0]
    c1=[305.96E-2, 305.96E-2, 297.16E-2, 297.16E-2, 248.15E-2, 248.15E-2]
    c2=[-764.41E-4, -764.41E-4, -769.38E-4, -769.38E-4, -438.59E-4, -438.59E-4]
    c3=[885.95E-6, 885.95E-6, 931.80E-6, 931.80E-6, 278.71E-6, 278.71E-6]
    c4=[-40.80E-7, -40.80E-7, -44.85E-7, -44.85E-7, -3.23E-7, -3.23E-7]
    c5=[0.60E+0, 0.60E+0, 0.20E+0, 0.20E+0, 0.60E+0, 0.60E+0]
    c6=[-0.16E+0, -0.16E+0, -0.15E+0, -0.15E+0, -0.57E+0, -0.57E+0]
    c7=[-2.13E-2, -2.13E-2, -7.51E-2, -7.51E-2, -2.61E-2, -2.61E-2]
    a0=[11.80E+0, 11.80E+0, 13.01E+0, 13.01E+0, 28.10E+0, 28.10E+0]
    av1=[2.23E-3, 2.23E-3, 6.16E-3, 6.16E-3, 1.85E-3, 1.85E-3]
    av2=[0.00E-5, 0.00E-5, 0.67E-5, 0.67E-5, 0.17E-5, 0.17E-5]

    #table 2
    #19v 19h 22v 22h 37v 37h
    epsilon0=[162.53E+0, 83.88E+0, 166.99E+0, 86.98E+0, 186.31E+0, 101.42E+0]
    epsilon1=[-25.70E-2, -52.22E-2, -34.08E-2, -59.52E-2, -56.37E-2, -85.88E-2]
    epsilon2=[17.29E-3, 18.76E-3, 17.35E-3, 19.38E-3, 14.81E-3, 20.76E-3]
    epsilon3=[-11.77E-5, -9.25E-5, -10.36E-5, -8.99E-5, -2.96E-5, -7.07E-5]
    epsilon4=[21.62E-1, -14.72E-1, 21.64E-1, -15.15E-1, 21.23E-1, -17.01E-1]
    epsilon5=[0.70E-2, 0.21E-2, 0.75E-2, 0.30E-2, 1.17E-2, 0.55E-2]
    epsilon6=[0.45E-1, -0.16E-1, 0.45E-1, -0.16E-1, 0.41E-1, -0.19E-1]
    epsilon7=[0.14E-4, -1.10E-4, 0.02E-4, -1.17E-4, -0.71E-4, -1.27E-4]
    betha=[-0.81E-4, 0.81E-4, -0.87E-4, 0.87E-4, -1.19E-4, 1.05E-4]
    mu=[0.41E-5, -0.13E-5, 0.54E-5, -0.16E-5, 1.25E-5, -0.29E-5]
    m1=[0.46E-3, 3.01E-3, 0.34E-3, 3.20E-3, -0.09E-3, 3.91E-3]
    m2=[3.78E-3, 7.50E-3, 3.48E-3, 7.39E-3, 2.38E-3, 7.00E-3]
    Tboffsets=[0.78, 2.1, 0.78, 0.0, -1.68, 0.13]
    #19v 19h 22v 22h 37v 37h
    Xi=[0.688, 0.688, 0.739, 0.739, 1.0, 1.0]

    #if stetement in case of no arguments
    """--------------------------
    e_ice=np.array([0.95,0.90,0.95,0.90,0.93,0.88])
    #test data
    V=5.0
    W=10.0
    L=1.0
    Ts=275.0
    Ti=260.0
    c_ice=0.5
    -----------------------------"""
    theta=53.1
    
    if c_ice < 0.0: c_ice=0.0
    if c_ice > 1.0: c_ice=1.0
    
    #order of coefficients and channels
    channels=['19V', '19H', '22V', '22H', '37V', '37H']
    #the isotropic Tb radiative transfer equation eq 15
    #Tb=TBU+tau*(emissivity*Ts+(1-emissivity)*(omega*TBD+tau*TBC))
    #eq 18a
    if V <= 48.0: Tv=273.16+0.8337*V-3.029E-5*(V**3.33)
    #eq 18b
    else:  Tv=301.16
    #eq 22 + 23
    Tl=(Ts+273.0)/2.0
    Al37=0.208*(1-0.026*(Tl-283))*L
    Al=np.array([0.2858*Al37, 0.2858*Al37, 0.3751*Al37, 0.3751*Al37, Al37, Al37])
    #eq 25
    q=theta-51.0
    t=Ts-273.16
    #eq 26
    W1=7.0
    W2=12.0
    #cosmic background
    TBC=2.7
    #initialisation
    TD=np.array([1.0,1.0,1.0,1.0,1.0,1.0])
    TU=np.array([1.0,1.0,1.0,1.0,1.0,1.0])
    A0=np.array([1.0,1.0,1.0,1.0,1.0,1.0])
    A0=np.array([1.0,1.0,1.0,1.0,1.0,1.0])
    Av=np.array([1.0,1.0,1.0,1.0,1.0,1.0])
    tau=np.array([1.0,1.0,1.0,1.0,1.0,1.0])
    TBU=np.array([1.0,1.0,1.0,1.0,1.0,1.0])
    TBD=np.array([1.0,1.0,1.0,1.0,1.0,1.0])
    E0=np.array([1.0,1.0,1.0,1.0,1.0,1.0])
    Ew=np.array([1.0,1.0,1.0,1.0,1.0,1.0])
    M1=np.array([1.0,1.0,1.0,1.0,1.0,1.0])
    M2=np.array([1.0,1.0,1.0,1.0,1.0,1.0])
    sigma=np.array([1.0,1.0,1.0,1.0,1.0,1.0])
    omega=np.array([1.0,1.0,1.0,1.0,1.0,1.0])
    emissivity=np.array([1.0,1.0,1.0,1.0,1.0,1.0])
    Tb=np.array([1.0,1.0,1.0,1.0,1.0,1.0])
    #loop over the channels
    for i in range(0,6):
        #print "now processing: "+channels[i]
        #eq17a
        TD[i]=c0[i]+c1[i]*V+c2[i]*(V**2)+c3[i]*(V**3)+c4[i]*(V**4)+c5[i]*(Ts-Tv)
        #eq17b
        TU[i]=TD[i]+c6[i]+c7[i]*V
        #eq 20
        A0[i]=(a0[i]/TD[i])**1.4
        #eq 21
        Av[i]=av1[i]*V+av2[i]*(V**2)
            
        #eq 19 transmitance through the atmosphere
        tau[i]=math.exp((-1.0/math.cos(math.radians(theta)))*(A0[i]+Av[i]+Al[i]))
        #eq 16a up-welling Tb
        TBU[i]=TU[i]*(1.0-tau[i])
        #eq 16b down-welling Tb
        TBD[i]=TD[i]*(1.0-tau[i])

        #eq 25 reference emissivity
        E0[i]=(epsilon0[i]+epsilon1[i]*t+epsilon2[i]*(t**2)+epsilon3[i]*(t**3)+epsilon4[i]*q+epsilon5[i]*t*q+epsilon6[i]*(q**2)+epsilon7[i]*(t**2)*q)/Ts
        #eq 27
        M1[i]=m1[i]+betha[i]*(theta-53.0)+mu[i]*(Ts-288.0)
        M2[i]=m2[i]+betha[i]*(theta-53.0)+mu[i]*(Ts-288.0)
        #eq 26 wind induced emissivity
        if W<=W1:   Ew[i]=M1[i]*W
        elif W1<W<W2: Ew[i]=M1[i]*W+0.5*(M2[i]-M1[i])*((W-W1)**2)/(W2-W1)
        elif W>=W2:   Ew[i]=M2[i]*W-0.5*(M2[i]-M1[i])*(W2+W1)
        else:   Ew[i]=M2[i]*W-0.5*(M2[i]-M1[i])*(W2+W1)
        #eq 29 sea surface slope variance: sigma**2
        sigma[i]=math.sqrt(5.22E-3*Xi[i]*W)
        if sigma[i]**2 > 0.07: sigma[i]=math.sqrt(0.07)
        #eq 28 reflection reduction factor due to surface roughness
        if i%2 == 0: omega[i]=1.0+2.5*(sigma[i]**2-68.0*sigma[i]**6)*tau[i]**3
        if i%2 != 0: omega[i]=1.0+6.1*(sigma[i]**2-68.0*sigma[i]**6)*tau[i]**2
        #eq 24
        emissivity[i]=E0[i]+Ew[i]
        #eq 16
        #return tbu+tau*((1-c_ice)*e*t_s+c_ice*e_ice*t_ice
        #              +(1-c_ice)*(1-e)*(omega*tbd+tau*tbc)
        #              +c_ice*(1-e_ice)*(tbd+tau*tbc));
        Tb[i]=TBU[i]+tau[i]*((1.0-c_ice)*emissivity[i]*Ts+c_ice*e_ice[i]*Ti[i]+(1.0-c_ice)*(1.0-emissivity[i])*(omega[i]*TBD[i]+tau[i]*TBC)+c_ice*(1.0-e_ice[i])*(TBD[i]+tau[i]*TBC))
    return Tb

