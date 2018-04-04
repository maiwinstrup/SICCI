#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 11:04:14 2018

@author: mwi
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def setupmodel(hsnow,hice,nsnow,nice,Ta,sal_snowice_surface,sal_ice_bottom,rho_snow,pcc_snow_top,pcc_snow_bottom,pcc_ice,makeplots='no'):
    # Set it all up:
    modelgrid = definegrid(hsnow,hice,nsnow,nice)
    Tsi,T = tempdistr(Ta,hsnow,hice,modelgrid)
    #salinity = saldistr(sal_snowice_surface,sal_ice_bottom,hsnow,hice,modelgrid)
    salinity = saldistr_Cshape(sal_snowice_surface,sal_ice_bottom,hsnow,hice,modelgrid)
    rho = densdistr(rho_snow,T,salinity,modelgrid)
    pcc = pccdistr(pcc_snow_top,pcc_snow_bottom,pcc_ice,hsnow,modelgrid)   
    
    modelgrid['Temperature'] = T
    modelgrid['Salinity'] = salinity
    modelgrid['rho'] = rho
    modelgrid['pcc'] = pcc # should be in mm
    
    Wi = np.zeros(len(modelgrid)) # No wetness, water content
    modelgrid['Wetness'] = Wi
    
    # Various parameters based on the above: (to be removed eventually)
    # Sitype
    # 3= first year ice
    # 4: multi-year ice
    mask = modelgrid["Type"]=='snow'
    sitype = np.zeros(len(modelgrid))
    sitype[mask] = 2 # snow
    sitype[~mask] = 4  # ice
    modelgrid['SnowIceType'] = sitype

    si = np.zeros(len(modelgrid))
    si[mask] = 0 # snow
    si[~mask] = 1 # ice
       
    num = range(len(modelgrid))
    num = np.array(num)
    
    di = np.array(modelgrid['LayerThickness'])
    
    #%% Check direction: layer 1 is lowest: set index so increasing towards the surface
    if modelgrid['Depth'][0]<modelgrid['Depth'][1]:
        modelgrid = modelgrid.reindex(index=modelgrid.index[::-1])
        T = T[::-1]
        salinity = salinity[::-1]
        rho = rho[::-1]
        pcc=pcc[::-1]
        sitype = sitype[::-1]
        si = si[::-1]
        Wi = Wi[::-1]
        di=di[::-1]

    #%% Plot profiles:
    if makeplots == 'yes':
        plt.figure(figsize=(10,10))
        plt.subplot(2,2,1)
        plt.plot(modelgrid['Depth'],rho)
        ax = plt.gca()
        ax.set_ylim([0.910, 0.930])
        plt.title('rho')
        plt.xlabel('Depth')
        
        plt.subplot(2,2,2)
        plt.plot(modelgrid['Depth'],salinity)
        plt.title('salinity')
        plt.xlabel('Depth')
        
        plt.subplot(2,2,3)
        plt.plot(modelgrid['Depth'],T)
        plt.title('T')
        plt.xlabel('Depth')
    
        plt.subplot(2,2,4)
        plt.plot(modelgrid['Depth'],pcc)
        plt.title('pcc')
        plt.xlabel('Depth')
        plt.show()

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
    # Salinity of snow surface is set to zero, and salinity is linearly 
    # increasing through the snow and ice, but with different gradients. 
    
    # Salinity gradient in ice:
    dsal_ice =(sal_ice_bottom-sal_ice_surface)/hice
    # Salinity gradient in snow (salinity at snow surface is 0):
    dsal_snow=sal_ice_surface/hsnow
    
    # Linear salinity profile:
    depth = modelgrid['Depth']
    salinity = np.zeros(len(depth))
    
    mask = modelgrid['Type']=='snow'
    salinity[mask] = dsal_snow*(depth[mask])
    salinity[~mask]= sal_ice_surface + dsal_ice*(depth[~mask]-hsnow)
    return salinity

def saldistr_Cshape(sal_ice_surface,sal_ice_bottom,hsnow,hice,modelgrid):
    # Compute the salinity distribution down the snow/ice pack
    # Salinity of snow surface is set to zero, and salinity is linearly 
    # increasing through the snow. In the ice, the salinity displays a 
    # C-shaped profile. 
    
    # Salinity gradient in ice:
    #dsal_ice =(sal_ice_bottom-sal_ice_surface)/hice
    # Salinity gradient in snow (salinity at snow surface is 0):
    dsal_snow=sal_ice_surface/hsnow
    
    depth = modelgrid['Depth']
    salinity = np.zeros(len(depth))

    # Linear salinity profile in snow:
    mask = modelgrid['Type']=='snow'
    salinity[mask] = dsal_snow*(depth[mask])
    #salinity[~mask]= sal_ice_surface + dsal_ice*(depth[~mask]-hsnow)
    
    # C-shaped in ice:
    dtop = 0.2 #05 # Thickness of top high-salinity layer
    # Her bruges så sal_ice_bottom som bulk salinity of ice. Top layer is 
    # twice as high as middle layer.
    Smid = (hice*sal_ice_bottom)/(2*2*dtop+(hice-dtop))
    mask_zone1 = (depth>hsnow) & (depth<hsnow+dtop)
    salinity[mask_zone1] = Smid*2
    mask_zone2 = (depth>=hsnow+dtop) & (depth<=hsnow+hice-dtop)
    salinity[mask_zone2] = Smid
    mask_zone3 = (depth>hsnow+hice-dtop)
    salinity[mask_zone3] = Smid*2
        
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
    
    if T>-2.0: #(T<0.0 and T>-2.0):
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
    rho=rhoi*(1.0-Vb-Va-porosity)+rhob*Vb #[g/cm3]
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
    #pcc[mask] = pcc_snow_top+3.0*hsnow/0.5*depth[mask]*(np.sqrt(depth[mask])-depth[mask])
    
    # Linear function in snow:
    # Gradient:
    dpcc_snow =(pcc_snow_bottom-pcc_snow_top)/hsnow
    pcc[mask] = dpcc_snow*(depth[mask])
    
    # And constant in the ice:
    pcc[~mask] = icepcc
    return pcc