#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
MEMLS
    Physical model for the radiation transmission through a layered snow/ice 
    pack at microwave frequencies. 
    The original MEMLS model was produced by Institute of Applied Physics, 
    University of Bern, Switzerland (Copyright 1997), and extended to apply to 
    snow-covered sea ice by rtt@dmi.dk and mwi@dmi.dk (2018).
   
@author: mwi
"""
import numpy as np
import setupmodel
import permittivity
import scattering
import reflectivity
import brightnesstemp

def memls(P,freqs,theta_in,polarization,makeplots='No'): 
    """
    Physical forward model for brightness temperatures for microwave radiation 
    with incoming angle (theta_in) and polarization at various frequencies 
    (freqs) for a layered snow/ice-pack described by parameters P.
    """
    
    #%% Set up depth profiles of temperature, salinity, etc:
    modelgrid = setupmodel.setupmodel(P["hsnow"],P["hice"],P["nsnow"],P["nice"],P["Ta"],
                                      P["sal_snowice_surface"],P["sal_ice_bottom"],
                                      P["rho_snow"],P["pcc_snow_top"],
                                      P["pcc_snow_bottom"],P["pcc_ice"],makeplots)[0]
  
    #%% Parameters used in the model: 
    Tgnd=271.35  # "Ground" (here: water) temperature 
    Tsky=0 # Downwelling sky radiation, in terms of brightness temperature
    
    # Ground-snow interface reflectivity, s0: 
    s0h=0.75
    s0v=0.25
   
    #%% Initialize arrays:
    eh = np.zeros(len(freqs))
    ev = np.zeros(len(freqs))
    Tbh = np.zeros(len(freqs))
    Tbv = np.zeros(len(freqs))
    Teffv = np.zeros(len(freqs))
    Teffh = np.zeros(len(freqs))
    dpeneh = np.zeros(len(freqs))
    dpenev = np.zeros(len(freqs))
    
    #%% Calculate surface brightness temperatures for each frequency band
    for ifreq in range(len(freqs)):
        # Dielectric permittivity profile:
        epsi, epsii = permittivity.permittivity(modelgrid['rho'],modelgrid['Temperature'],modelgrid['Wetness'],modelgrid['Salinity'],modelgrid['SnowIceType'],freqs[ifreq])
        
        # Clean out coherent layers:
        #        Ti = T
        #        [rnum,rroi,repsi,repsii,rtei,rsih,rsiv,rdi,rdei,rTi,rpci,rWi,rgai,rsitype,rsal] = reducemodel.slred(num,rho,epsi,epsii,tei,sih,siv,di,dei,Ti,pci,freqs[ifreq],Wi,gai,sitype,sal)
       
        # Depth profile of scatter and absorption coefficients:
        gs6, gbih, gbiv, ga2i = scattering.scattercoeff(modelgrid,epsi,epsii,freqs[ifreq])
        # OBS: epsi, epsii bliver delvist udregnet igen i scatter coefficienterne...

        # Layer/interface reflectivities and layer transmissivity: 
        rih,riv,tih,tiv,sih,siv,theta_eff,dpeneh[ifreq],dpenev[ifreq] = reflectivity.reflectivity( \
            np.array(modelgrid['Depth']),np.array(modelgrid['LayerThickness']),
            theta_in,epsi,epsii,gs6,gbih,gbiv,ga2i,s0h,s0v)
        # OBS: penetration depth is not correct!
        
        # Surface brightness temperatures, and emissivity:
        Tbh[ifreq], Tbv[ifreq], eh[ifreq], ev[ifreq], Teffh[ifreq], Teffv[ifreq] = \
            brightnesstemp.brightnesstemp(rih,riv,tih,tiv,sih,siv,np.array(modelgrid.Temperature), Tgnd, Tsky)
       
    return Tbh, Tbv, eh, ev, Teffh, Teffv, dpeneh, dpenev