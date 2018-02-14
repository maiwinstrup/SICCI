#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 11:04:14 2018

@author: mwi
"""
import pandas as pd
import numpy as np

def setupmodel(hsnow,hice,nsnow,nice,Ta,sal_snowice_surface,sal_ice_bottom,rho_snow,pcc_snow_top,pcc_snow_bottom,pcc_ice):
    # Set it all up:
    modelgrid = definegrid(hsnow,hice,nsnow,nice)
    Tsi,T = tempdistr(Ta,hsnow,hice,modelgrid)
    salinity = saldistr(sal_snowice_surface,sal_ice_bottom,hsnow,hice,modelgrid)
    rho = densdistr(rho_snow,T,salinity,modelgrid)
    pcc = pccdistr(pcc_snow_top,pcc_snow_bottom,pcc_ice,hsnow,modelgrid)   
    
    # Various parameters based on the above: (to be removed eventually)
    mask = modelgrid["Type"]=='snow'
    sitype = np.zeros(len(modelgrid))
    sitype[mask] = 2 # snow
    sitype[~mask] = 4  # ice
    si = np.zeros(len(modelgrid))
    si[mask] = 0 # snow
    si[~mask] = 1 # ice
    
    Wi = np.zeros(len(modelgrid)) # No wetness, water content

    num = range(len(modelgrid))
    di = np.array(modelgrid['LayerThickness'])

    return modelgrid, T, salinity, rho, pcc, sitype, si, Wi, num, di

def definegrid(hs,hi,nsnow,nice):
    # Define model grid:
    # Using a fixed number of layers within the snow and ice fraction
    
    # Gridsize within the snow and ice fraction, respectively: 
    ds=hs/nsnow # cell thickness in snow
    di=hi/nice # cell thicness in ice
    
    modelgrid = pd.DataFrame()
    # Set material type (snow/ice):
    modelgrid['Type'] = ['snow']*nsnow+['ice']*nice
    # Layer thickness:
    layerthickness = np.array([ds]*nsnow+[di]*nice)
    modelgrid['LayerThickness'] = layerthickness
    # Layer depth (centered):
    modelgrid['Depth'] = np.cumsum(layerthickness)-layerthickness/2
    return modelgrid

def tempdistr(Ta,hsnow,hice,modelgrid):
    # Compute temperature profile (T) and snow/ice interface temperature  (Tsi)
    # Assuming thermal equilibrium and no change in material with depth, which 
    # gives rise to a step-wise linear temperature profile within snow pack and 
    # ice. 
    
    kice=2.0 # Thermal conductivity of ice [W/mK]
    ksnow=0.31 # Thermal conductivity of snow [W/mK]
    Twater=271.35 # Temperature of water under ice

    # Temperature at snow/ice interface:
    f=(ksnow*hice)/(kice*hsnow)
    Tsi=(Twater+f*Ta)/(f+1.0)

    # Temperature gradients:
    dTsnow=(Tsi-Ta)/hsnow
    dTice=(Twater-Tsi)/hice

    # Temperature at gridpoints:
    depth = np.array(modelgrid['Depth'])
    mask = depth<hsnow
    T = np.zeros(len(depth))
    T[mask] = Ta+dTsnow*depth[mask]
    T[~mask] = Tsi+dTice*(depth[~mask]-hsnow)
    return Tsi, T

def saldistr(sal_ice_surface,sal_ice_bottom,hsnow,hice,modelgrid):
    # Compute the salinity distribution down the snow/ice pack
    # Salinity of snow is set to zero, and salinity is linearly increasing 
    # through the ice. 
    
    # Salinity gradient in ice:
    dsal=(sal_ice_bottom-sal_ice_surface)/hice;

    # Linear salinity profile:
    depth = modelgrid['Depth']
    salinity = np.zeros(len(depth))
    
    mask = modelgrid['Type']=='snow'
    salinity[mask] = 0
    salinity[~mask]=dsal*(depth[~mask]-hsnow)
    return salinity

def densdistr(rho_snow,T,salinity,modelgrid):
    # Compute the snow/ice density profile.
    
    # Snow density is assumed constant value (rho_snow)
    # Ice density is calculated based on T, salinity (and porosity). 
    # Porosity is here set to 0:
    porosity = 0.0
    
    # Snow density:
    depth = modelgrid['Depth']
    rho = np.zeros(len(depth))
    mask = modelgrid['Type']=='snow'
    rho[mask] = rho_snow
    
    # Sea ice density profile:
    # Using formulation from Eiken (chap 2), and Cox & Weeks (1983).
    Tice = T[~mask]
    salice = salinity[~mask]
    rho_ice = np.zeros(len(Tice))
    for i in xrange(0,len(Tice)):
        rho_ice[i] = icedensity(Tice[i],salice[i],porosity)
    rho[~mask] = rho_ice
    
    return rho
    
def icedensity(TK,salinity,porosity):
    # Sea ice density formulation from Eiken (chap 2), and Cox & Weeks (1983).
        
    # Temperature in Celcius:
    T = TK-273.15
    
    rhoi=0.917-1.403e-4*T
    Sb=1.725-18.756*T-0.3964*T**2
    rhob=1.0+0.0008*Sb
    
    if (T<0.0 and T>-2.0):
        ba=np.array([-4.1221e-2,-18.407,5.8402e-1,2.1454e-1])
        bb=np.array([9.0312e-2,-1.6111e-2,1.2291e-4,1.3603e-4])
        F1 = ba[0] + ba[1]*T + ba[2]*T**2 + ba[3]*T**3
        F2 = bb[0] + bb[1]*T + bb[2]*T**2 + bb[3]*T**3
        Vb = rhoi*salinity/(F1-rhoi*salinity*F2)
        
    else:
        if (T>-22.9 and T<-2): a = np.array([-4.732,-22.45,-6.397e-1,-1.074e-2])
        if T<=-22.9: a = np.array([9.899e3,1.309e3,5.527e1,7.160e-1])
        F1 = a[0] + a[1]*T + a[2]*T**2 + a[3]*T**3
        Vb = rhoi*salinity/(F1+salinity*rhoi-salinity*rhob)
        
    rhop=(rhoi*(1.0-Vb)+rhob*Vb)
    Va=(1.0-(rhop/rhoi))+(rhop*salinity/Sb)*(1.0/rhoi - 1.0/rhob)

    # Corresponding density:
    rho=1000.0*(rhoi*(1.0-Vb-Va-porosity)+rhob*Vb)
    return rho

def pccdistr(pcc_snow_top,pcc_snow_bottom,icepcc,hsnow,modelgrid):
    # Scattering profile:
    #snow_top=0.07; #the correlation length of new snow [mm]
    #snow_bottom=0.3; #the correlation length of old snow [mm] not used??
    #icepcc=1.5; #the correlation length of multiyear ice [mm]

    depth = modelgrid['Depth'] # Mid-layer depths
    mask = modelgrid['Type']=='snow'
    # Using some weird function in snowpack...
    pcc = np.zeros(len(depth))
    pcc[mask] = pcc_snow_top+3.0*hsnow/0.5*depth[mask]*(np.sqrt(depth[mask])-depth[mask])
    # And constant in the ice:
    pcc[~mask] = icepcc
    return pcc