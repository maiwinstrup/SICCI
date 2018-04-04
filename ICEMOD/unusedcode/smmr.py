#!/usr/bin/env python
# -*- coding: utf-8 -*-
#F. T. Wentz. A model function for ocean microwave brightness temperatures. Journal of Geophysical Research 88(C3), 1892-1908, 1983.
#Ts: SST [K]
#U:  windspeed [m/s]
#V:  columnar water vapour [g/cm2]
#L:  columnar liquid water [mg/cm2]
#Ta: surface air temperature [K]
#c:  ice concentration [1/100]
#Ei: ice emissivity at 18v 18h 36v 36h
import sys, os, string, math
import numpy as np

def smmr(Ts,U,V,L,Ta,Ti,c,Ei):
#def smmr():

    #Ts=275.0
    #U=0.0
    V=0.1*V
    L=0.1*L
    #Ta=275.0
    #c=0.0
    #Ei=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

    theta=50.0                  #Nimbus 7, SMMR +/- 0.5deg. [deg]
    thetar=50.0*math.pi/180.0   #in radians [rad]
    Uf=100.0*0.047*U            #covert U to Uf [cm/s]
    Tc=2.76                     #cosmic temperature [K].
    lamda=5.9                   #lapse rate [K/km]
    aTa=289.0                   #average global air temperature [K]
    T=Ts-273.16                 #SST in C [C]

    #sea ice emissivity
    #Ei=[0.95, 0.90, 0.93, 0.88]

    #Regression coefficients for specular brightness temperature, table 4 [K,K-1,K-2,K/deg]
    s=[[1.3759E2, 2.368E-1, 1.565E-2, -2.311E-4, 2.03],\
       [0.7107E2, 0.891E-1, 1.000E-2, -1.476E-4, -1.28],\
       [1.4452E2, -0.336E-1, 2.076E-2, -2.497E-4, 2.05],\
       [0.7559E2, -0.935E-1, 1.371E-2, -1.661E-4, -1.32],\
       [1.5750E2, -3.936E-1, 2.285E-2, -2.048E-4, 2.08],\
       [0.8444E2, -3.675E-1, 1.657E-2, -1.568E-4, -1.4],\
       [1.6252E2, -4.916E-1, 2.237E-2, -1.775E-4, 2.10],\
       [0.8802E2, -4.545E-1, 1.699E-2, -1.477E-4, 1.43],\
       [1.8493E2, -7.405E-1, 1.694E-2, -0.539E-4, 2.11],\
       [1.0524E2, -7.666E-1, 1.718E-2, -1.033E-4, -1.59]]

    #regression coefficients for wind-induced emissivity, table 5: m [s/cm]
    m=[[1.55E-4, 4.90E-4],\
       [4.58E-4, 6.02E-4],\
       [1.41E-4, 4.61E-4],\
       [5.16E-4, 7.09E-4],\
       [2.66E-4, 2.66E-4],\
       [7.05E-4, 7.05E-4],\
       [2.68E-4, 2.68E-4],\
       [7.60E-4, 7.60E-4],\
       [2.80E-4, 2.80E-4],\
      [10.51E-4,10.51E-4]]

    #regression coefficients for wind-induced emissivity, table 5: b [s/cm deg]
    b=[-0.94E-5, 0.88E-5, -1.34E-5, 1.39E-5, -1.68E-5, 1.63E-5, -1.79E-5, 1.82E-5, -2.54E-5, 2.24E-5] 

    #brightness temperature of diffusely scattered atmospheric radiation, table 1 [s/cm]
    w=[0.70E-3, 1.18E-3, 1.34E-3, 2.37E-3, 1.23E-3, 2.33E-3, 0.81E-3, 1.73E-3, 0.75E-3, 1.82E-3]

    #Atmospheric absorption coefficients, table3 [millinapers/(g/cm2)]
    ao=[8.29E-3, 8.29E-3, 8.59E-3, 8.59E-3, 9.72E-3, 9.72E-3, 10.78E-3, 10.78E-3, 29.04E-3, 29.04E-3]
    av=[1.05E-3, 1.05E-3, 2.47E-3,2.47E-3, 13.62E-3, 13.62E-3, 45.45E-3, 45.45E-3, 23.90E-3, 23.90E-3] 
    aL=[0.112, 0.112, 0.401, 0.401, 1.125, 1.125, 1.36, 1.36, 2.224, 2.224] 

    #effective height, table 3 [km]
    He=[7.4, 7.4, 6.0, 6.0, 4.4, 4.4, 4.5, 4.5, 4.5, 4.5]

    #Temperature sensitivities of absorption coefficients, table 2, [K-1]
    Qo=[-1.14E-2, -1.14E-2, -1.14E-2, -1.14E-2, -1.14E-2, -1.14E-2, -1.13E-2, -1.13E-2, -1.11E-2, -1.11E-2]
    Qv=[-0.65E-2, -0.65E-2, -0.61E-2, -0.61E-2, -0.36E-2, -0.36E-2, -0.06E-2, -0.06E-2, -0.65E-2, -0.65E-2]
    QL=[-2.85E-2, -2.85E-2, -2.82E-2, -2.82E-2, -2.73E-2, -2.73E-2, -2.68E-2, -2.68E-2, -2.33E-2, -2.33E-2]

    dEa=np.array([1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0])
    Es=np.array([1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0])
    dE=np.array([1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0])
    Ew=np.array([1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0])
    E=np.array([1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0])
    ETs=np.array([1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0])
    yo=np.array([1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0])
    yv=np.array([1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0])
    yL=np.array([1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0])
    tau=np.array([1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0])
    d=np.array([1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0])
    Tb=np.array([1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0])

    for i in range(0,10): #going through the channels 6.6v, 6.6h, 10v, 10h, 18v, 18h, 21v, 21h, 37v, 37h

        #the surface emissivity
        #eq. 40-42
        if Uf <= 65: dEa[i]=m[i][0]*Uf #Uf<65cm/s
        elif 65 < Uf <= 75: dEa[i]=m[i][0]*Uf+0.05*(m[i][1]-m[i][0])*(Uf-65)**2 #65 < Uf < 75cm/s
        else: dEa[i]=m[i][1]*Uf-70*(m[i][1]-m[i][0])
        
        Es[i]=(s[i][0]+s[i][1]*T+s[i][2]*T**2+s[i][3]*T**3+s[i][4]*(theta-49))/Ts #eq.38
        dE[i]=dEa[i]+b[i]*Uf*(theta-49) #eq.39
        Ew[i]=Es[i]+dE[i] #eq.35
        Tw=Ts
        #Ti=0.4*Ta+0.6*273.15

        E[i]=(1-c)*Ew[i]+c*Ei[i]
        ETs[i]=(1-c)*Ew[i]*Tw+c*Ei[i]*Ti[i]

        #the atmosphere
        #air temperature coefficient, eq.30
        yo[i]=1+Qo[i]*(Ta-aTa) 
        yv[i]=1+Qv[i]*(Ta-aTa)
        yL[i]=1+QL[i]*(Ta-aTa)

        #atmospheric transmittance, tau, eq.28
        tau[i]=math.exp((-1.0/math.cos(thetar))*(ao[i]*yo[i]+av[i]*yv[i]*V+aL[i]*yL[i]*L))

        #effective emission depth, d, eq.29
        d[i]=He[i]*((tau[i]-1-tau[i]*math.log(tau[i]))/(math.log(tau[i])-tau[i]*math.log(tau[i]))) 

        #brightness temperature, Tb, eq.27
        #Tb[i] = tau[i]*(Ew[i]*Ts + (1.0+w[i]*Uf)*(1.0-Ew[i])*((1.0-tau[i])*(Ta-lamda*d[i])+tau[i]*Tc)) + (1.0-tau[i])*(Ta-lamda*(He[i]-d[i]))
        Tb[i] = (c)*(tau[i]*(Ei[i]*Ti[i] + (1.0-Ei[i])*((1.0-tau[i])*(Ta-lamda*d[i])+tau[i]*Tc)) + (1.0-tau[i])*(Ta-lamda*(He[i]-d[i]))) + \
             (1.0-c)*(tau[i]*(Ew[i]*Ts    + (1.0-Ew[i])*((1.0-tau[i])*(Ta-lamda*d[i])+tau[i]*Tc)) + (1.0-tau[i])*(Ta-lamda*(He[i]-d[i])))
        #Tb[i]=(tau[i])*(ETs[i]+(1-c)*(1+w[i]*Uf)*(1-E[i])*((1-tau[i])*(Ta-lamda*d[i]))+c*(1-E[i])*((1-tau[i])*(Ta-lamda*d[i])+tau[i]*Tc))+(1-tau[i])*(Ta-lamda*(He[i]-d[i]))
    return Tb
