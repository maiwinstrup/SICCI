#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Collection of subroutines for calculation of brightness temperatures, effective 
temperatures, and emissivities of the snow and ice pack. 

Based off the original MEMLS code (1997) by the Institute of Applied Physics, 
University of Bern, Switzerland.

@author: Mai Winstrup, Danish Meteorological Institute, Copenhagen, Denmark, 2018
"""
import numpy as np

def brightnesstemp(rih,riv,tih,tiv,sih,siv,Ti,Tground,Tsky):
    """
    Surface brightness temperatures, effective temperatures, and emissivities 
    of the snow and ice pack. 
        
    [Tbh,Tbv,eh,ev,Teffh,Teffv] = brightnesstemp(rih,riv,tih,tiv,sih,siv,T,Tgnd,Tsky)
        rih, riv:  internal reflectivity of layer j, h/v-polarization
        tih, tiv:  internal transmissivity of layer j, h/v-polarization
        sih, siv:  layer interface reflectivity on top of layer j, h/v-polarization
        Ti:    temperature of layer i [K]
        Tground:  "Ground" temperature (i.e. water temperature) [K] 
        Tsky:  Brightness temperature of downwelling sky radiation [K]
        
    Uses: upwellingTb
    """
    
    #%% Brightness temperatures:
    # Upwelling brightness temperatures at every layer in the stack, 
    # horisontal and vertical polarization:
    Dh = upwellingTb(rih,tih,sih,Ti,Tground,Tsky) 
    Dv = upwellingTb(riv,tiv,siv,Ti,Tground,Tsky)
    
    # Brightness temperatures at surface:
    Tbh = (1-sih[-1])*Dh[-1] + sih[-1]*Tsky
    Tbv = (1-siv[-1])*Dv[-1] + siv[-1]*Tsky
        
    #%% Total emissitivy of snow/ice and background:
    # Equations from e.g. Wiesmann and Mätzler, 1999.
    # Upwelling brightness temperatures when Tsky = 0:
    Dh0 = upwellingTb(rih,tih,sih,Ti,Tground,0)
    Dv0 = upwellingTb(riv,tiv,siv,Ti,Tground,0)
    Tbh0 = (1-sih[-1])*Dh0[-1]
    Tbv0 = (1-siv[-1])*Dv0[-1]
    
    # Upwelling brightness temperatures when Tsky = 100:
    Dh100 = upwellingTb(rih,tih,sih,Ti,Tground,100) 
    Dv100 = upwellingTb(riv,tiv,siv,Ti,Tground,100)
    Tbh100 = (1-sih[-1])*Dh100[-1] + sih[-1]*100
    Tbv100 = (1-siv[-1])*Dv100[-1] + siv[-1]*100
    
    # Total emissivities:
    ev = 1-(Tbv100 - Tbv0)/100.
    eh = 1-(Tbh100 - Tbh0)/100.
    
    #%% Effective temperatures:
    Teffv = Tbv0/ev
    Teffh = Tbh0/eh
    return Tbh, Tbv, eh, ev, Teffh, Teffv
        
def upwellingTb(ri,ti,si,Ti,Tground,Tsky):
    """
    Upwelling brightness temperatures for each layer in the snow and ice stack. 
    Calculated according to the equations in Matzler (1996), Note 6. 

    [D] = upwellingTb(ri,ti,si,Ti,Tground,Tsky)
        D:       Upwelling brightness temperatures for each layer in the stack
        ri:      Internal reflectivity of layer i
        ti:      Internal transmissivity of layer i
        si:      Layer interface reflectivity on top of layer i (incl. ground)
        Ti:      Temperature of layer i [K]
        Tground: Brightness temperature of the "ground" (aka. ice-water interface below the sea ice) 
        Tsky:    Downwelling brightness temperature of the sky
    """

    N = len(ri) # Number of ice/snow layers
    ei = 1 - ri - ti # Emissivity of layer j
    
    if N < 1: 
        print 'ERROR: No scattering layer'
        return

    elif N == 1: 
        k1 = (1-ri[0]*si[0])*(1-ri[0]*si[1]) - ti[0]*si[0]*ti[0]*si[1]
        D = ti[0]*si[0] * ((1-si[0])*ri[0]*Tground + (1-si[1])*Tsky*ti[0] + ei[0]*Ti[0]) / (k1) + \
            (1-ri[0]*si[0]) * ((1-si[0])*Tground*ti[0] + (1-si[1])*Tsky*ri[0] + ei[0]*Ti[0]) / (k1)
      
    else:
        # Forming matrices M1, M2, M3, M4, E, F, according to Matzler (1996),
        # Note 6s, p.30-31:
        M1 = np.diag(ri*si[0:N],k=0) + np.diag(ti[0:N-1]*(1-si[1:N]),k=1)
        M2 = np.diag(ti*si[1:N+1],k=0) + np.diag(ri[1:N]*(1-si[1:N]),k=-1)
        M3 = np.diag(ti*si[0:N],k=0) + np.diag(ri[0:N-1]*(1-si[1:N]),k=1)
        M4 = np.diag(ri*si[1:N+1],k=0) + np.diag(ti[1:N]*(1-si[1:N]),k=-1)
        
        E = ei*Ti
        E[0] = E[0] + ri[0]*(1-si[0])*Tground # First element
        E[-1] = E[-1] + ti[-1]*(1-si[-1])*Tsky # Last element
        
        F = ei*Ti
        F[0] = F[0] + ti[0]*(1-si[0])*Tground # First element
        F[-1] = F[-1] + ri[-1]*(1-si[-1])*Tsky # Last element
        
        I = np.eye(N)

        # Compute vector with upwelling brightness temperatures, D:
        try: 
            M5 = M3.dot(np.linalg.inv(I-M1).dot(M2)) + M4
            D = np.linalg.inv(I-M5).dot(M3.dot(np.linalg.inv(I-M1).dot(E))+F)
        except np.linalg.LinAlgError as err:
            if 'Singular matrix' in err.message:
                D = np.empty(Ti.shape)
                D[:] = np.nan
            else:
                raise
    return D