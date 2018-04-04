#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 09:39:23 2018

@author: mwi
"""

import numpy as np
import matplotlib.pyplot as plt

import MEMLS.setupmodel as setupmodel
import MEMLS.memls as memls

import MEMLS_Clara.memls_functions_clara as memls_functions_clara
import MEMLS_Clara.memls_functions_clara0 as memls_functions_clara0

#%% Model parameters:
# Model grid: 
hsnow = 0.3 # snow thickness [m]
hice = 3.0 # ice thickness [m]
nsnow = 10 # number of layers in snow (as slabs)
nice = 20 # number of layers in ice

# Temperature profile:
Ta = 273.0-30 # Surface temperature, K
#Tsi,T = setupmodel.tempdistr(Ta,hsnow,hice,modelgrid)
# Assuming thermal equilibrium, giving rise to piecewise linear temperature profile
# OBS: The thermal conductivity is VERY important for the retrieval

# Salinity distribution:
sal_snowice_surface = 0.0 # Salinity at the surface of the ice 
sal_ice_bottom = 3.0; # Salinity of the bottom ice [ppt]
# Assuming linear increase in salinity with depth in the ice; salinity of snow 
# is set equal to zero.

# Estimate using equation from Cox&Weeks, 1974:
# As function of ice thicknesss, for cold sea ice sampled during growth season:
#if hice <= 0.4: sal_mean = 14.24-19.39*hice
#else: sal_mean = 7.88-1.59*hice
# Salinity at bottom will be twice the mean salinity (linear function)
#sal_ice_bottom = sal_mean*2  # permil, [ppt/psu]

# Density profile:
rho_snow = 360./1000 # Estimate from paper [g/cm3]
#rho = setupmodel.densdistr(rho_snow,T,salinity,modelgrid)

# Scattering profile in snow and ice:
pcc_snow_top=0.07; # Correlation length of new snow [mm]
pcc_snow_bottom=0.3; # Correlation length of old snow [mm] - not used???! 
pcc_ice=1.5; # Correlation length of multiyear ice [mm]

# Model parameters:
P = dict([("hsnow",hsnow),("hice",hice),("nsnow",nsnow),("nice",nice),("Ta",Ta),
          ("sal_snowice_surface",sal_snowice_surface),("sal_ice_bottom",sal_ice_bottom),
          ("rho_snow",rho_snow),("pcc_snow_top",pcc_snow_top),
          ("pcc_snow_bottom",pcc_snow_bottom),("pcc_ice",pcc_ice)])

#%%
freqs= [1.4, 6.9, 10.7, 18.7, 23.8, 36.5, 50.0, 89.0]
theta= (55*np.pi)/180 # incidence angle in radians
polarization = 0  # p = polarization of incident beam

Tbh, Tbv, eh, ev, Teffh, Teffv, dpeneh, dpenev = memls.memls(P,freqs,theta,polarization,makeplots='yes')

# Test med flere lag:
#P1 = P
#P1["nsnow"]=nsnow*10
#P1["nice"]=nice*10
#Tbh0, Tbv0, eh0, ev0, Teffh0, Teffv0, dpeneh0, dpenev0 = memls.memls(P1,freqs,theta,polarization)

#%% Plot resultater:
plt.figure()
plt.plot(freqs,Tbh,'.-b',label='TbH')
plt.plot(freqs,Tbv,'.-r',label='TbV')
#plt.plot(freqs,Tbh0,'--g',label='TbH1')
#plt.plot(freqs,Tbv0,'--g',label='TbV1')
#plt.xlabel('Frequency')
#plt.title('Brightness temperatures')
#plt.legend()

#%% Results using Claras model:
modelgrid,T, salinity, rho, pcc, sitype, si, Wi, num, di = \
    setupmodel.setupmodel(P["hsnow"],P["hice"],P["nsnow"],P["nice"],P["Ta"],
                                      P["sal_snowice_surface"],P["sal_ice_bottom"],
                                      P["rho_snow"],P["pcc_snow_top"],
                                      P["pcc_snow_bottom"],P["pcc_ice"])
roi = np.array(modelgrid["rho"]*1000)
pci = np.array(modelgrid["pcc"])
sal = np.array(modelgrid["Salinity"])
Ti = np.array(modelgrid["Temperature"])

x, eh, ev, yh, yv, Teffh, Teffv, pend_h, pend_v = memls_functions_clara.memls_mod(freqs,num,di,Ti,Wi,roi,pci,sal,sitype,si)[0:9]
x, eh, ev, yh0, yv0, Teffh, Teffv, pend_h, pend_v = memls_functions_clara0.memls_mod(freqs,num,di,Ti,Wi,roi,pci,sal,sitype,si)[0:9]

#plt.figure()
plt.plot(freqs,yh,'.-c',label='TbH-clara')
plt.plot(freqs,yv,'.-m',label='TbV-clara')
plt.plot(freqs,yh0,'--k',label='TbH-clara-original')
plt.plot(freqs,yv0,'--k',label='TbV-clara-original')
#plt.xlabel('Frequency')
#plt.title('Claras brightness temperatures')
plt.legend()
plt.show()

##%% Import real data:
#year = 2015
#month = 3 # (3,
#day = 24 #19
#
#filename = './RRDP_v2.0/NERSC_OIB/QSCAT-vs-SMAP-vs-SMOS-vs-ASCAT-vs-AMSR2-vs-ERA-vs-NERSC-OIB-'+ str(year)+ str(month).zfill(2) + str(day).zfill(2) +'.text'
#
#df=pd.read_csv(filename,sep=',',skiprows=[0], na_values=['noval',' noval','  noval','   noval','    noval','     noval','      noval'])
#    # List all variables
#    #df.keys()
#    #list the variables where to drop nans
#nndf=df[[u'#latitude',u'longitude',u'6.9GHzH', u'6.9GHzV', u'7.3GHzH', u'7.3GHzV', u'10.7GHzH' , u'10.7GHzV',\
#        u'18.7GHzH', u'18.7GHzV', u'23.8GHzH', u'23.8GHzV',u'36.5GHzH', u'36.5GHzV', u'89.0GHzH', u'89.0GHzV',\
#        u'sd_mean',u'sd_std',u't2m']].dropna()
#    
#Tb7h0 = nndf['6.9GHzH'].as_matrix()
#Tb7v0 = nndf['6.9GHzV'].as_matrix()
#Tb11h0 = nndf['10.7GHzH'].as_matrix()    
#Tb11v0 = nndf['10.7GHzV'].as_matrix()
#Tb19h0 = nndf['18.7GHzH'].as_matrix()    
#Tb19v0 = nndf['18.7GHzV'].as_matrix()
#Tb24h0 = nndf['23.8GHzH'].as_matrix()    
#Tb24v0 = nndf['23.8GHzV'].as_matrix()
#Tb37h0 = nndf['36.5GHzH'].as_matrix()    
#Tb37v0 = nndf['36.5GHzV'].as_matrix()
#Tb89h0 = nndf['89.0GHzH'].as_matrix()    
#Tb89v0 = nndf['89.0GHzV'].as_matrix()          
#
#plt.figure()
#for x in xrange(0,len(Tb7h0)):  
#    freqs = [6.9, 10.7, 18.7, 23.8, 36.5, 89]
#    Tbh = np.array([Tb7h0[x],Tb11h0[x],Tb19h0[x],Tb24h0[x],Tb37h0[x],Tb89h0[x]])
#    Tbv = np.array([Tb7v0[x],Tb11v0[x],Tb19v0[x],Tb24v0[x],Tb37v0[x],Tb89v0[x]])
#      
#    plt.plot(freqs,Tbh,'.-b')
#    plt.plot(freqs,Tbv,'.-r')