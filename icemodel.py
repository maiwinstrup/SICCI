#!/usr/bin/python
"""
The sea ice model? Contains various functions for the sea ice model. 
Ts_amrs(<Tb>): Effective temperature (?) - vertical only?

Usage:
 >>>   
"""

import numpy as np

def liseequation(Tb): 
    # Lise's equations for snowdepth
    sd = 1.7701 + 0.017462*Tb[0] - 0.02801*Tb[4] + 0.0040926*Tb[8]
    #sd = max(0.0001,sd) # Never negative (and never 0)
    
    return sd

def icethickness(Tb):
    # Rasmus' model for ice thickness
    Tb6v = np.float(Tb[0]) #6v 
    Tb6h = np.float(Tb[1]) #6h
    Tb10v = np.float(Tb[2]) #10v
    Tb19v = np.float(Tb[4]) #19v
    Tb19h = np.float(Tb[5]) #19h
    Tb37v = np.float(Tb[8]) #37v

    # Spectral gradient at various frequencies:
    GR0610v=(Tb10v-Tb6v)/(Tb10v+Tb6v)
    GR0618v=(Tb19v-Tb6v)/(Tb19v+Tb6v)
    GR0618h=(Tb19h-Tb6h)/(Tb19h+Tb6h)
    GR0636v=(Tb37v-Tb6v)/(Tb37v+Tb6v)

    # Rasmus equation for ice thickness
    ITsim = 4.81691061767 - 47.3872663622*GR0610v - 52.8250933515*GR0618v + \
        49.9619055273*GR0618h + 33.4020392591*GR0636v + 0.142250481385*Tb6h - \
        0.115797387369*Tb19h - 0.0385743856297*Tb37v
    ITsim=max(0.3,ITsim)
    ITsim=min(3.5,ITsim)
    return ITsim

def Ts_amsr(Tb):
    # Calculating an array of effective temperatures at (approx) the following frequencies:
    #  [Teff6v,Teff10v,Teff19v,Teff23v,Teff37v,Teff50v,Teff52, Teff89v]
    # We assume that effective temperatures are identical for vertical and horizontal polarization.
    
    # Lise's equation for snow depth:
    SDsim = liseequation(Tb)
   
    # Snow-ice interface temperature, Tsi:
    # Can be estimated from both the 6 and 10Ghz channels. Using the average. 
    Tb6v = np.float(Tb[0]) #6v
    Tb10v = np.float(Tb[2]) #10v
    Tsi6 = 1.144*Tb6v - 0.815/SDsim - 27.08
    Tsi10 = 1.101*Tb10v - 0.999/SDsim - 14.28
    Tsi = 0.5*(Tsi6 + Tsi10)
#   Tsi=max(200.0,Tsi)
#   Tsi=min(273.15,Tsi)

    # Effective temperature:
    Teff6v  = 0.888*Tsi + 30.245
    Teff10v = 0.901*Tsi + 26.569
    Teff19v = 0.920*Tsi + 21.536
    Teff23v = 0.932*Tsi + 18.417
    Teff37v = 0.960*Tsi + 10.902
    Teff50v = 0.989*Tsi + 2.959
    Teff89v = 1.0604*Tsi - 16.384
    Teff=np.array([Teff6v,Teff10v,Teff19v,Teff23v,Teff37v,Teff50v,Teff50v,Teff89v])
    return Teff

def ev_ice(Tb,IST,SD):
    # Calculating vertical emissivity at approximately the following frequencies: 
    #   [e6v,e10v,e18v,e23v,e37v,e50v,e52v,e89v]
    # Some are replaced with neighboring frequencies, i.e. using instead:    
    #   [e6v,e10v,e18v,e18v,e37v,e37v,e37v,e89v]
    
    # Rasmus equation for ice thickness
    ITsim = icethickness(Tb)

    # Rasmus forward model for simulated brightness temperatures:
    Tbsim=simulateTb(IST,SD,ITsim)
    
    # Calculating the emissivity at various wavelengths (=Tb/Teff)
    Teff = Ts_amsr(Tb) # Includes channel 50GHz
        
    e6v=Tbsim[0]/Teff[0]
    e10v=Tbsim[2]/Teff[1]
    e18v=Tbsim[4]/Teff[2]
    e37v=Tbsim[6]/Teff[4]
    e89v=Tbsim[8]/Teff[7]
    ev=np.array([e6v,e10v,e18v,e18v,e37v,e37v,e37v,e89v]) 
    
    # Check limits:
    for i in range(len(ev)):
        ev[i]=max(0.0, ev[i])
        ev[i]=min(1.0, ev[i])
    return ev

def eh_ice(Tb,IST,SD):
    # Calculating horizontal emissivity at following frequencies: 
    #   [e6h,e10h,e18h,e23h,e37h,e50h,e52h,e89h]
    # Some are replaced with neighboring frequencies, i.e. using instead:    
    #   [e6h,e10h,e18h,e18h,e37h,e37h,e37h,e89h]
   
     # Rasmus equation for ice thickness
    ITsim = icethickness(Tb)

    # Rasmus forward model for simulated brightness temperatures:
    Tbsim=simulateTb(IST,SD,ITsim) # At following wavelengths: T6v,T6h,T10v,T10h,T18v,T18h,T36v,T36h,T89v,T89h

    # Calculating the emissivity at various wavelengths (=Tb/Teff)
    Teff = Ts_amsr(Tb) 
    e6h=Tbsim[1]/Teff[0]
    e10h=Tbsim[3]/Teff[1]
    e18h=Tbsim[5]/Teff[2]
    e37h=Tbsim[7]/Teff[4]
    e89h=Tbsim[9]/Teff[7]
    eh=np.array([e6h,e10h,e18h,e18h,e37h,e37h,e37h,e89h]) 
   
    # Check for limits:
    for i in range(len(eh)):
        eh[i]=max(0.0, eh[i])
        eh[i]=min(1.0, eh[i])
    return eh

def simulateTb(IST,SD,IT):
    # Rasmus forward model for simulated brightness temperatures (March+April 2013):
    # A regression model for sea ice surface emissions. 
    # Calculated from ice surface temperatures (IST), snow depth (SD), and ice thickness (IT).
    T6vsim = 151.981535394 + 0.39827296166*IST + 23.3600203008*SD - 3.03183834111*IT
    T6hsim = 55.2623240539 + 0.687577210357*IST + 12.9621301692*SD -1.66486943272*IT
    T10vsim = 145.878105173 + 0.435432823207*IST + 0.743658800361*SD - 4.20200228328*IT
    T10hsim = 45.1075848344 + 0.753868056619*IST - 18.7322202698*SD - 3.49097027752*IT
    T18vsim = 138.073033941 + 0.479944899424*IST - 71.814765826*SD - 5.57027693189*IT
    T18hsim = 78.4247537087 + 0.641671517779*IST - 85.1843265565*SD - 5.34143981869*IT
    T36vsim = 123.102065469 + 0.526904188123*IST - 216.727298154*SD - 4.03697464853*IT
    T36hsim = 131.862853412 + 0.429489868157*IST - 214.352191714*SD - 3.03524993537*IT
    T89vsim = 2.52567864217 + 0.902202528995*IST - 180.427137566*SD + 1.90480465092*IT
    T89hsim = 31.1206976877 + 0.743826485118*IST - 184.806381816*SD + 3.19723383624*IT
    Tb=np.array([T6vsim,T6hsim,T10vsim,T10hsim,T18vsim,T18hsim,T36vsim,T36hsim,T89vsim,T89hsim])
    return Tb

def icemod(Tb,IST,SD):
    # Rasmus equation for ice thickness
    ITsim = icethickness(Tb)    
 
    # Rasmus forward model for simulated brightness temperatures:
    Tbsim=simulateTb(IST,SD,ITsim)
    
    # Effective temperatures at the "correct" wavelengths:
    Teff = Ts_amsr(Tb) # Includes channel 23GHz and 50GHz (but 23 GHz exist in data!)
    e6v=Tbsim[0]/Teff[0]
    e6h=Tbsim[1]/Teff[0]
    e10v=Tbsim[2]/Teff[1]
    e10h=Tbsim[3]/Teff[1]
    e18v=Tbsim[4]/Teff[2]
    e18h=Tbsim[5]/Teff[2]
    e37v=Tbsim[6]/Teff[4]
    e37h=Tbsim[7]/Teff[4]
    e89v=Tbsim[8]/Teff[6]
    e89h=Tbsim[9]/Teff[6]
    e=np.array([e6v, e6h, e10v, e10h, e18v, e18h, e37v, e37h, e89v, e89h])

    # Check for limits:
    for i in range(len(e)):
        e[i]=max(0.0, e[i])
        e[i]=min(1.0, e[i])
        
    return IST, SD, Tbsim, e, Teff