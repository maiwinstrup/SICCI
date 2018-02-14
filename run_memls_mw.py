#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 09:39:23 2018

@author: mwi
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import memls_functions as memls # den oprindelige, men med freqs som input
import setupmodel

#%% Set up model:
# Model grid: 
hsnow = 0.3 # snow thickness [m]
hice = 3.0 # ice thickness [m]
nsnow = 10 # number of layers in snow
nice = 20 # number of layers in ice

#modelgrid = setupmodel.definegrid(hsnow,hice,nsnow,nice)

# Temperature profile:
Ta = 273.0-10 # Surface temperature, K
#Tsi,T = setupmodel.tempdistr(Ta,hsnow,hice,modelgrid)
# Assuming thermal equilibrium, giving rise to piecewise linear temperature profile
# OBS: The thermal conductivity is VERY important for the retrieval

# Salinity distribution:
sal_snowice_surface = 0.0 # Salinity at the surface of the ice 
sal_ice_bottom = 3.0; # Salinity of the bottom ice [ppt]

# Estimate using equation from Cox&Weeks, 1974:
# As function of ice thicknesss, for cold sea ice sampled during growth season:
if hice <= 0.4: sal_mean = 14.24-19.39*hice
else: sal_mean = 7.88-1.59*hice
# Salinity at bottom will be twice the mean salinity (linear function)
sal_ice_bottom = sal_mean*2  # permil

#salinity = setupmodel.saldistr(sal_snowice_surface,sal_ice_bottom,hsnow,hice,modelgrid)
# Assuming linear increase in salinity with depth in the ice; salinity of snow 
# is set equal to zero.

# Density profile:
rho_snow = 360. # Estimate from paper
#rho = setupmodel.densdistr(rho_snow,T,salinity,modelgrid)

# Scattering profile in snow and ice:
pcc_snow_top=0.07; # Correlation length of new snow [mm]
pcc_snow_bottom=0.3; # Correlation length of old snow [mm] - not used???! 
pcc_ice=1.5; # Correlation length of multiyear ice [mm]
#pcc = setupmodel.pccdistr(pcc_snow_top,pcc_snow_bottom,pcc_ice,hsnow,modelgrid)
#Chech units: roi = roi/1000? di = di/100?

modelgrid, T, salinity, rho, pcc, sitype, si, Wi, num, di = setupmodel.setupmodel(hsnow,hice,nsnow,nice,Ta,sal_snowice_surface,sal_ice_bottom,rho_snow,pcc_snow_top,pcc_snow_bottom,pcc_ice)

#%%
plt.figure()
plt.plot(modelgrid['Depth'],rho)
ax = plt.gca()
ax.set_ylim([910, 930])
plt.title('rho')

plt.figure()
plt.plot(modelgrid['Depth'],salinity)
plt.title('salinity')

plt.figure()
plt.plot(modelgrid['Depth'],T)
plt.title('T')

plt.figure()
plt.plot(modelgrid['Depth'],pcc)
plt.title('pcc')

#%%
freqs= [1.4, 6.9, 10.7, 18.7, 23.8, 36.5, 50.0, 89.0]
frequency, eh, ev, TBH, TBV, TeffH, TeffV, pend_h, pend_v = memls.memls_mod(freqs,np.array(num),di,T,Wi,rho,pcc,salinity,sitype,si)

#%%
plt.figure()
plt.plot(frequency,TBH,'.-b')

plt.figure()
plt.plot(frequency,TBV,'.-r')
      
plt.figure()
plt.plot(eh)
plt.plot(ev)

plt.figure()
plt.plot(TeffH)
plt.plot(TeffV)

plt.figure()
plt.plot(pend_h)
plt.plot(pend_v)          