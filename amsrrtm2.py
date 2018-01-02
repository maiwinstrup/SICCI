#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""---------------------------------------------------------------------------
AMSR Ocean Algorithm; Frank J. Wentz, Thomas Meissner; Remote
Sensing Systems; Version 2; November 2; 2000.
Tb=f(V,W,L,Ts,Ti,c_ice)
V: columnar water vapor [mm]
W: windspeed over water [m/s]
L: columnar cloud liquid water [mm]
Ts: sea surface temperature [K]
Ti: ice surface temperature [K]
c_ice: ice concentration [0-1]
-------------------------------------------------------------------------------"""
import sys, os, string, math, cmath
import numpy as np

def amsr(V,W,L,Ta,Ts,Ti_amsrv,Ti_amsrh,c_ice,e_icev,e_iceh):
    theta=55.0
    c_ice=c_ice-1.0
    #V=0.0
    #W=0.0
    #L=0.0
    #Ts=275.0
    #Ti=260.0
    #c_ice=0.0
    #e_icev=np.array([1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0])
    #e_iceh=np.array([1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0])

    frequencies=np.array([6.93, 10.65, 18.70, 23.80, 36.50, 50.30, 52.80, 89.00])

    b0=np.array([239.50E+0,  239.51E+0,  240.24E+0,  241.69E+0,  239.45E+0,  242.10E+0,  245.87E+0,  242.58E+0])
    b1=np.array([213.92E-2,  225.19E-2,  298.88E-2,  310.32E-2,  254.41E-2,  229.17E-2,  250.61E-2,  302.33E-2])
    b2=np.array([-460.60E-4, -446.86E-4, -725.93E-4, -814.29E-4, -512.84E-4, -508.05E-4, -627.89E-4, -749.76E-4])
    b3=np.array([457.11E-6,  391.82E-6,  814.50E-6,  998.93E-6,  452.02E-6,  536.90E-6,  759.62E-6,  880.66E-6])
    b4=np.array([-16.84E-7,  -12.20E-7,  -36.07E-7,  -48.37E-7,  -14.36E-7,  -22.07E-7,  -36.06E-7,  -40.88E-7])
    b5=np.array([0.50E+0,     0.54E+0,    0.61E+0,    0.20E+0,    0.58E+0,    0.52E+0,    0.53E+0,    0.62E+0])
    b6=np.array([-0.11E+0,   -0.12E+0,   -0.16E+0,   -0.20E+0,   -0.57E+0,   -4.59E+0,  -12.52E+0,   -0.57E+0])
    b7=np.array([-0.21E-2,   -0.34E-2,   -1.69E-2,   -5.21E-2,   -2.38E-2,   -8.78E-2,  -23.26E-2,   -8.07E-2])
    ao1=np.array([8.34E-3,    9.08E-3,   12.15E-3,   15.75E-3,   40.06E-3,  353.72E-3, 1131.76E-3,   53.35E-3])
    ao2=np.array([-0.48E-4,  -0.47E-4,   -0.61E-4,   -0.87E-4,   -2.00E-4,  -13.79E-4,   -2.26E-4,   -1.18E-4])
    av1=np.array([0.07E-3,    0.18E-3,    1.73E-3,    5.14E-3,    1.88E-3,    2.91E-3,    3.17E-3,    8.78E-3])
    av2=np.array([0.00E-5,    0.00E-5,   -0.05E-5,    0.19E-5,    0.09E-5,    0.24E-5,    0.27E-5,    0.80E-5])

    aL1=np.array([0.0078, 0.0183, 0.0556, 0.0891,  0.2027,  0.3682,  0.4021,  0.9693])
    aL2=np.array([0.0303, 0.0298, 0.0288, 0.0281,  0.0261,  0.0236,  0.0231,  0.0146])
    aL3=np.array([0.0007, 0.0027, 0.0113, 0.0188,  0.0425,  0.0731,  0.0786,  0.1506])
    aL4=np.array([0.0000, 0.0060, 0.0040, 0.0020, -0.0020, -0.0020, -0.0020, -0.0020])
    aL5=np.array([1.2216, 1.1795, 1.0636, 1.0220,  0.9546,  0.8983,  0.8943,  0.7961])

    r0v=np.array([-0.27E-3,  -0.32E-3,  -0.49E-3,  -0.63E-3,  -1.01E-3, -1.20E-3, -1.23E-03, -1.53E-3])
    r0h=np.array([0.54E-3,   0.72E-3,   1.13E-3,   1.39E-3,   1.91E-3,  1.97E-3,  1.97E-03,  2.02E-3])
    r1v=np.array([-0.21E-4,  -0.29E-4,  -0.53E-4,  -0.70E-4,  -1.05E-4, -1.12E-4, -1.13E-04, -1.16E-4])
    r1h=np.array([0.32E-4,   0.44E-4,   0.70E-4,   0.85E-4,   1.12E-4,  1.18E-4,  1.19E-04,  1.30E-4])
    r2v=np.array([-2.10E-5,  -2.10E-5,  -2.10E-5,  -2.10E-5,  -2.10E-5, -2.10E-5, -2.10E-05, -2.10E-5])
    r2h=np.array([-25.26E-6, -28.94E-6, -36.90E-6, -41.95E-6, -54.51E-6, -5.50E-5, -5.50E-5,  -5.50E-5])
    r3v=np.array([0.00E-6,   0.08E-6,   0.31E-6,   0.41E-6,   0.45E-6,  0.35E-6,  0.32E-06, -0.09E-6])
    r3h=np.array([0.00E-6,  -0.02E-6,  -0.12E-6,  -0.20E-6,  -0.36E-6, -0.43E-6, -0.44E-06, -0.46E-6])

    m1v=np.array([0.00020, 0.00020, 0.00140, 0.00178, 0.00257, 0.00260, 0.00260, 0.00260])
    m1h=np.array([0.00200, 0.00200, 0.00293, 0.00308, 0.00329, 0.00330, 0.00330, 0.00330])
    m2v=np.array([0.00690, 0.00690, 0.00736, 0.00730, 0.00701, 0.00700, 0.00700, 0.00700])
    m2h=np.array([0.00600, 0.00600, 0.00656, 0.00660, 0.00660, 0.00660, 0.00660, 0.00660])

    TD=np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
    TU=np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
    AO=np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
    AV=np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
    AL=np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
    tau=np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
    TBU=np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
    TBD=np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
    llambda=np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
    epsilon=np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0], dtype=np.complex)
    rho_H=np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0], dtype=np.complex)
    rho_V=np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0], dtype=np.complex)
    R_0H=np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
    R_0V=np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
    R_geoH=np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
    R_geoV=np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
    F_H=np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
    F_V=np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
    R_H=np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
    R_V=np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
    OmegaH=np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
    OmegaV=np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
    T_BOmegaH=np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
    T_BOmegaV=np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
    emissivityh=np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
    emissivityv=np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
    term=np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
    Delta_S2=np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
    Tv=np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
    Th=np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])

    T_C=2.7

    if c_ice < 0.0: c_ice=0.0
    if c_ice > 1.0: c_ice=1.0

    Ts_mix=c_ice*Ta+(1.0-c_ice)*Ts

    Tl=(Ts_mix+273.0)/2.0

    #eq 27
    Tv=273.16+0.8337*V-3.029E-5*(V**3.33)
    if V > 48:  Tv=301.16
    G = 1.05*(Ts_mix-Tv)*(1-((Ts_mix-Tv)**2)/1200.0)
    if math.fabs(Ts_mix-Tv) > 20: G = (Ts_mix-Tv)*14/math.fabs(Ts_mix-Tv)

    epsilon_R=4.44 # this value is from wentz and meisner, 2000, p. 28
    s=35.0
    ny=0.012 # Klein and Swift is using 0.02 which is giving a higher epsilon_R (4.9)
    light_speed=3.00E10
    free_space_permittivity=8.854E-12
    #eq 43
    epsilon_S = (87.90*math.exp(-0.004585*(Ts-273.15)))*(math.exp(-3.45E-3*s+4.69E-6*s**2+1.36E-5*s*(Ts-273.15)))
    #eq 44
    lambda_R = (3.30*math.exp(-0.0346*(Ts-273.15)+0.00017*(Ts-273.15)**2))-(6.54E-3*(1-3.06E-2*(Ts-273.15)+2.0E-4*(Ts-273.15)**2)*s)
    #eq 41
    C=0.5536*s
    #eq 42
    delta_t=25.0-(Ts-273.15)
    #eq 40
    qsi=2.03E-2+1.27E-4*delta_t+2.46E-6*delta_t**2-C*(3.34E-5-4.60E-7*delta_t+4.60E-8*delta_t**2)
    #eq 39
    sigma=3.39E9*(C**0.892)*math.exp(-delta_t*qsi)

    for i in range(0,8):
        #eq26
        TD[i]=b0[i]+b1[i]*V+b2[i]*V**2+b3[i]*V**3+b4[i]*V**4+b5[i]*G
        TU[i]=TD[i]+b6[i]+b7[i]*V
        #eq 28
        AO[i]=ao1[i]+ao2[i]*(TD[i]-270.0)
        #eq 29
        AV[i]=av1[i]*V+av2[i]*V**2
        #eq 33
        AL[i]=aL1[i]*(1.0-aL2[i]*(Tl-283.0))*L
        #eq 22
        tau[i]=math.exp((-1.0/math.cos(math.radians(theta)))*(AO[i]+AV[i]+AL[i])) 
        #eq 24
        TBU[i]=TU[i]*(1.0-tau[i])
        TBD[i]=TD[i]*(1.0-tau[i])

        llambda[i]=(light_speed/(frequencies[i]*1E9))

        #eq 35
        epsilon[i]=epsilon_R+((epsilon_S-epsilon_R)/(1.0+((cmath.sqrt(-1)*lambda_R)/llambda[i])**(1.0-ny)))-((2.0*cmath.sqrt(-1)*sigma*llambda[i])/light_speed)
        #eq.45
        rho_H[i]=(math.cos(math.radians(theta))-cmath.sqrt(epsilon[i]-math.sin(math.radians(theta))**2))/\
                 (math.cos(math.radians(theta))+cmath.sqrt(epsilon[i]-math.sin(math.radians(theta))**2))
        rho_V[i]=(epsilon[i]*math.cos(math.radians(theta))-cmath.sqrt(epsilon[i]-math.sin(math.radians(theta))**2))/\
                 (epsilon[i]*math.cos(math.radians(theta))+cmath.sqrt(epsilon[i]-math.sin(math.radians(theta))**2))
        #eq46
        R_0H[i]=np.absolute(rho_H[i])**2
        R_0V[i]=np.absolute(rho_V[i])**2+(4.887E-8-6.108E-8*(Ts-273.0)**3)

        #eq 57
        R_geoH[i]=R_0H[i]-(r0h[i]+r1h[i]*(theta-53.0)+r2h[i]*(Ts-288.0)+r3h[i]*(theta-53.0)*(Ts-288.0))*W
        R_geoV[i]=R_0V[i]-(r0v[i]+r1v[i]*(theta-53.0)+r2v[i]*(Ts-288.0)+r3v[i]*(theta-53.0)*(Ts-288.0))*W
        #eq.60
        W_1=7.0
        W_2=12.0
        if W<W_1: F_H[i]=m1h[i]*W
        elif W_1<W<W_2: F_H[i]=m1h[i]*W+0.5*(m2h[i]-m1h[i])*((W-W_1)**2)/(W_2-W_1)
        else: F_H[i]=m2h[i]*W-0.5*(m2h[i]-m1h[i])*(W_2+W_1)
        W_1=3.0
        W_2=12.0
        if W<W_1: F_V[i]=m1v[i]*W
        elif W_1<W<W_2: F_V[i]=m1v[i]*W+0.5*(m2v[i]-m1v[i])*((W-W_1)**2)/(W_2-W_1)
        else: F_V[i]=m2v[i]*W-0.5*(m2v[i]-m1v[i])*(W_2+W_1)

        R_H[i]=(1-F_H[i])*R_geoH[i]
        R_V[i]=(1-F_V[i])*R_geoV[i]

        emissivityh[i]=1-R_H[i]
        emissivityv[i]=1-R_V[i]

        if i >= 4: Delta_S2[i]=5.22E-3*W
        else: Delta_S2[i]=5.22E-3*(1-0.00748*(37.0-frequencies[i])**1.3)*W
        if Delta_S2[i]>0.069: Delta_S2[i]=0.069
        #eq.62
        term[i]=Delta_S2[i]-70.0*Delta_S2[i]**3
        OmegaH[i]=(6.2-0.001*(37.0-frequencies[i])**2)*term[i]*tau[i]**2.0
        OmegaV[i]=(2.5+0.018*(37.0-frequencies[i]))*term[i]*tau[i]**3.4
        #eq.61
        T_BOmegaH[i]=((1+OmegaH[i])*(1-tau[i])*(TD[i]-T_C)+T_C)*R_H[i] 
        T_BOmegaV[i]=((1+OmegaV[i])*(1-tau[i])*(TD[i]-T_C)+T_C)*R_V[i]
        
        Th6=TBU[0]+tau[0]*((1.0 - c_ice)*emissivityh[0]*Ts+c_ice*e_iceh[0]*Ti_amsrh[0]+(1.0 - c_ice)*\
            (1.0 - emissivityh[0])*(T_BOmegaH[0]+tau[0]*T_C)+c_ice*(1.0 - e_iceh[0])*(TBD[0]+tau[0]*T_C))

        Tv6=TBU[0]+tau[0]*((1.0 - c_ice)*emissivityv[0]*Ts+c_ice*e_icev[0]*Ti_amsrv[0]+(1.0 - c_ice)*\
            (1.0 - emissivityv[0])*(T_BOmegaV[0]+tau[0]*T_C)+c_ice*(1.0 - e_icev[0])*(TBD[0]+tau[0]*T_C))
        
        Th10=TBU[1]+tau[1]*((1.0 - c_ice)*emissivityh[1]*Ts+c_ice*e_iceh[1]*Ti_amsrh[1]+(1.0 - c_ice)*\
             (1.0 - emissivityh[1])*(T_BOmegaH[1]+tau[1]*T_C)+c_ice*(1.0 - e_iceh[1])*(TBD[1]+tau[1]*T_C))

        Tv10=TBU[1]+tau[1]*((1.0 - c_ice)*emissivityv[1]*Ts+c_ice*e_icev[1]*Ti_amsrv[1]+(1.0 - c_ice)*\
             (1.0 - emissivityv[1])*(T_BOmegaV[1]+tau[1]*T_C)+c_ice*(1.0 - e_icev[1])*(TBD[1]+tau[1]*T_C))
        
        Th18=TBU[2]+tau[2]*((1.0 - c_ice)*emissivityh[2]*Ts+c_ice*e_iceh[2]*Ti_amsrh[2]+(1.0 - c_ice)*(1.0 - emissivityh[2])*\
             (T_BOmegaH[2]+tau[2]*T_C)+c_ice*(1.0 - e_iceh[2])*(TBD[2]+tau[2]*T_C))

        Tv18=TBU[2]+tau[2]*((1.0 - c_ice)*emissivityv[2]*Ts+c_ice*e_icev[2]*Ti_amsrv[2]+(1.0 - c_ice)*(1.0 - emissivityv[2])*\
             (T_BOmegaV[2]+tau[2]*T_C)+c_ice*(1.0 - e_icev[2])*(TBD[2]+tau[2]*T_C))
        
        Th23=TBU[3]+tau[3]*((1.0 - c_ice)*emissivityh[3]*Ts+c_ice*e_iceh[3]*Ti_amsrh[3]+(1.0 - c_ice)*(1.0 - emissivityh[3])*\
             (T_BOmegaH[3]+tau[3]*T_C)+c_ice*(1.0 - e_iceh[3])*(TBD[3]+tau[3]*T_C))

        Tv23=TBU[3]+tau[3]*((1.0 - c_ice)*emissivityv[3]*Ts+c_ice*e_icev[3]*Ti_amsrv[3]+(1.0 - c_ice)*(1.0 - emissivityv[3])*\
             (T_BOmegaV[3]+tau[3]*T_C)+c_ice*(1.0 - e_icev[3])*(TBD[3]+tau[3]*T_C))
        
        Th36=TBU[4]+tau[4]*((1.0 - c_ice)*emissivityh[4]*Ts+c_ice*e_iceh[4]*Ti_amsrh[4]+(1.0 - c_ice)*(1.0 - emissivityh[4])*\
             (T_BOmegaH[4]+tau[4]*T_C)+c_ice*(1.0 - e_iceh[4])*(TBD[4]+tau[4]*T_C))

        Tv36=TBU[4]+tau[4]*((1.0 - c_ice)*emissivityv[4]*Ts+c_ice*e_icev[4]*Ti_amsrv[4]+(1.0 - c_ice)*(1.0 - emissivityv[4])*\
             (T_BOmegaV[4]+tau[4]*T_C)+c_ice*(1.0 - e_icev[4])*(TBD[4]+tau[4]*T_C))
        
        Th50=TBU[5]+tau[5]*((1.0 - c_ice)*emissivityh[5]*Ts+c_ice*e_iceh[5]*Ti_amsrh[5]+(1.0 - c_ice)*(1.0 - emissivityh[5])*\
             (T_BOmegaH[5]+tau[5]*T_C)+c_ice*(1.0 - e_iceh[5])*(TBD[5]+tau[5]*T_C))

        Tv50=TBU[5]+tau[5]*((1.0 - c_ice)*emissivityv[5]*Ts+c_ice*e_icev[5]*Ti_amsrv[5]+(1.0 - c_ice)*(1.0 - emissivityv[5])*\
             (T_BOmegaV[5]+tau[5]*T_C)+c_ice*(1.0 - e_icev[5])*(TBD[5]+tau[5]*T_C))
        
        Th52=TBU[6]+tau[6]*((1.0 - c_ice)*emissivityh[6]*Ts+c_ice*e_iceh[6]*Ti_amsrh[6]+(1.0 - c_ice)*(1.0 - emissivityh[6])*\
             (T_BOmegaH[6]+tau[6]*T_C)+c_ice*(1.0 - e_iceh[6])*(TBD[6]+tau[6]*T_C))

        Tv52=TBU[6]+tau[6]*((1.0 - c_ice)*emissivityv[6]*Ts+c_ice*e_icev[6]*Ti_amsrv[6]+(1.0 - c_ice)*(1.0 - emissivityv[6])*\
             (T_BOmegaV[6]+tau[6]*T_C)+c_ice*(1.0 - e_icev[6])*(TBD[6]+tau[6]*T_C))
        
        Th89=TBU[7]+tau[7]*((1.0 - c_ice)*emissivityh[7]*Ts+c_ice*e_iceh[7]*Ti_amsrh[7]+(1.0 - c_ice)*(1.0 - emissivityh[7])*\
             (T_BOmegaH[7]+tau[7]*T_C)+c_ice*(1.0 - e_iceh[7])*(TBD[7]+tau[7]*T_C))

        Tv89=TBU[7]+tau[7]*((1.0 - c_ice)*emissivityv[7]*Ts+c_ice*e_icev[7]*Ti_amsrv[7]+(1.0 - c_ice)*(1.0 - emissivityv[7])*\
             (T_BOmegaV[7]+tau[7]*T_C)+c_ice*(1.0 - e_icev[7])*(TBD[7]+tau[7]*T_C))
        
        Tv=[Tv6,Tv10,Tv18,Tv23,Tv36,Tv50,Tv52,Tv89]
        Th=[Th6,Th10,Th18,Th23,Th36,Th50,Th52,Th89]
        Tb=np.array([Tv6,Th6,Tv10,Th10,Tv18,Th18,Tv23,Th23,Tv36,Th36,Tv89,Th89])
    #return Tv, Th
    return Tb
