#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Subroutines for computing depth profiles of scatter and absorption coefficients 
in the snow and ice pack. 

@author: mwi
"""
import numpy as np
import permittivity
import brine
import matplotlib.pyplot as plt

def scattercoeff(modelgrid,epsi,epsii,freq):
    """
    Depth profiles of scatter and absorption coefficients in the snow and ice 
    pack.
    
    [sccoeff_h, sccoeff_v, abscoeff] = scattercoeff(modelgrid,freq)
        sccoeff_h: Profile of scatter coefficients, horizontal polarization
        sccoeff_v: Profile of scatter coefficients, vertical polarization
        abscoeff:  Profile of absorption coefficients
        modelgrid: Depth profiles of input and calculated depth profiles
        freq:      Frequency [GHz]
        
        pci should be given in mm, is used as an effective size parameter for inclusions
    Uses: ...
    
    gb is scattering coefficient, ga2i absorption coefficient (?)
    #ga2i is two-flux absorption coefficient
    """
    
    rho = modelgrid['rho'] # Density [g/cm³]
    pci = modelgrid['pcc'] # 
    T = modelgrid['Temperature'] # Temperature [K]
    W = modelgrid['Wetness'] # Wetness []
    sal = modelgrid['Salinity']
    
    # Complex refraction index, n, from permittivity profile:
    #n, gamma = refractionconstants(epsi,epsii,freq) # not used
    
    # Absorption coefficient, based on the permittivity profile:
    # 6-flux absorption coefficients:
    gai = absorption(epsi,epsii,freq)  # this is the POWER absorption coefrfiicent, kappa_alpha
    # is gamma_a in note 2, Matzler.     
    
    #  Compute the 6-flux scatter coefficients:
    # New snow: 
    # Using scattering coefficient corresponding to partly recrystallized snow, 
    # a linear combination of the improved? born approximation and xx??
    sitype = modelgrid['SnowIceType']
    mask = (sitype == 1)
    scattermodel = 'linearmix' 
    gs6 = np.zeros(len(T)) 
    gs6[mask] = sixfluxsnow(rho[mask],T[mask],W[mask],pci[mask],epsi[mask],epsii[mask],freq,scattermodel)
    
    # Recrystallized snow:
    mask = (sitype == 2)
    scattermodel = 'shells' # Scattering coefficient of fully recrystalized snow (iborn)
    gs6[mask] = sixfluxsnow(rho[mask],T[mask],W[mask],pci[mask],epsi[mask],epsii[mask],freq,scattermodel)

    # First-year ice: 
    lambda0=2.99793/(freq*10) # Free space wavelength [m]
    k = 2*np.pi/lambda0 # Free space wave number
    
    mask = (sitype == 3)
    gs6[mask] = firstyearice(T[mask],sal[mask],pci[mask],k,freq)
    
    # Multi-year ice:  
    mask = (sitype == 4)
    gs6[mask] = multiyearice(rho[mask],pci[mask],epsi[mask],epsii[mask],k,freq)
    
    # Compute the 2-flux absorption and scattering coefficients:
    gbih,gbiv,ga2i = twofluxcomponents(gs6,epsi,gai)
    return gs6, gbih, gbiv, ga2i

def sixfluxsnow(rho,T,W,pci,epsi,epsii,freq,scattermodel):
    """
    Scattering and absorption coefficients for dry snow computed based on 
    structural parameters. Several options available for the scatter model.  
    Partly recalculating epsi in this function :-/
        
    [gs6] = snow(rho,T,W,pci,epsi,freq,scattermodel)
        gs6:    6-flux scattering coefficient
        rho:    density [g/cm^3]
        T:      temperature [K]
        W:      wetness
        pci:    correlation length [mm]
        epsi:   relative permittivity (real part)
        freq:   frequency [GHz]
        scattermodel: scattering coefficient algorithm
    
    Version history:
    1.0b    wi 15.7.95
    1.0     wi 23.9.97 bug fixed
    1.1     wi 26.9.97 latest fit on experimental data was added (option 7)
    1.2     wi 13.10.97 option 8 added, adapted scattering of a shell/sphere to note 9/ver2 
    1.3     wi  4.11.97 option 9, 10 and 11 added 
    1.4     wi 27.05.98 born approximation added (borna.m)
    12.03.07 bug in iborn shperes fixed (rtt)  
    """ 
    
    epsi_ice,epsii_ice = permittivity.epureice(T,freq)
    
    rho_ice = 0.9167
    vi = rho/rho_ice
 
    lambda0=2.99793/(freq*10) # Free space wavelength [m]
    k = 2*np.pi/lambda0 # Free space wave number
   
    #%% Estimate the 6-flux scattering coefficients:
    if scattermodel  == '1':
        gs6 = ((130 * (freq/50)**2.7) * pci**3) / (rho**1.3 + 0.001)
  
    elif scattermodel == '2': # fit from 26/8-97 on v-polarized data > 11 GHz (?)
        gs6 = 0.0704*(pci**2.32)*(freq**1.68)*rho**(-0.63)

    elif scattermodel == 'spheres': # No. "4"
        # From Mätzler, J. Appl. Phys. 83(11) 6111-6117 eqs 27+32 iborn
        # Host material is assumed to be pure ice, with air bubbles as inclusions. 
        # This is for dry snow, no impurities. 
        epsi_air = 1
        epsii_air = 0
        #        [epsi_eff,epsii_eff] = permittivity.sphericalinclusions(epsi_air,epsii_air,epsi_ice,epsii_ice,vi,'ImprovedBorn')
        # This is similar to already calculated epsi, epsii - IF calculated in this way)
        # Using precalculated values:
        epsi_eff = epsi
        epsii_eff = epsii
        gs6 = improvedbornscattering(epsi_air,epsii_air,epsi_ice,epsii_ice,epsi_eff,epsii_eff,0.001*pci,vi,k,freq,'spheres') #gs6 = (3./32) * (0.001*pci)**3 * k**4 * vi*(1-vi) * abs((2*epsi_eff+epsi_air)*(epsi_ice-epsi_air)/(2*epsi_eff+epsi_ice))**2  # eq. 33, Matzler (1998)
        
    elif scattermodel == 'shells': # No. "5" - bliver brugt
        # For shells (new and recrystalized snow): 
        # Mätzler, J. Appl. Phys. 83(11) 6111-6117 eq 39 iborn
        epsi_air = 1
        epsii_air = 0
        [epsi_eff,epsii_eff] = permittivity.sphericalinclusions(epsi_air,epsii_air,epsi_ice,epsii_ice,vi,'ImprovedBorn_shells')# Matzler 1998, eq. 39
        gs6 = improvedbornscattering(epsi_air,epsii_air,epsi_ice,epsii_ice,epsi_eff,epsii_eff,0.001*pci,vi,k,freq,'shells')
        
        #gs6 = abs(2./3 + 1./(3.*epsi_ice**2))*(0.001*pci)*k**2*vi*(1-vi)*(epsi_ice-1)**2/(16*epsi_eff) # eq. 42
  
    elif scattermodel == 'linearmix': # 6 - bliver brugt
        # Scattering coefficient of only partly recrystalized snow, linear combination of iborn shells and spheres:
        epsi_air = 1
        epsii_air = 0
        [epsi_spheres,epsii_spheres] = permittivity.sphericalinclusions(epsi_air,epsii_air,epsi_ice,epsii_ice,vi,'ImprovedBorn')
        # gs6_spheres = (3./32) * (0.001*pci)**3 * k**4 * vi*(1-vi) * abs((2*epsi_eff+epsi_air)*(epsi_ice-epsi_air)/(2*epsi_eff+epsi_ice))**2  # eq. 33, Matzler (1998)
        gs6_spheres = improvedbornscattering(epsi_air,epsii_air,epsi_ice,epsii_ice,epsi_spheres,epsii_spheres,0.001*pci,vi,k,freq,'spheres')
        
        [epsi_shells,epsii_shells] = permittivity.sphericalinclusions(epsi_air,epsii_air,epsi_ice,epsii_ice,vi,'ImprovedBorn_shells')# Matzler 1998, eq. 39
        #gs6_shells = abs(2./3 + 1./(3.*epsi_ice**2))*(0.001*pci)*k**2*vi*(1-vi)*(epsi_ice-1)**2/(16*epsi_eff)
        gs6_shells = improvedbornscattering(epsi_air,epsii_air,epsi_ice,epsii_ice,epsi_shells,epsii_shells,0.001*pci,vi,k,freq,'shells')

        a = 0.1664
        b = 0.2545
        gs6 = a*gs6_spheres+b*gs6_shells
        
    elif scattermodel == 'fit7': # fit from 26/9-97
        gs6 = 73.21 * (pci**3)*((freq/50)**2.68)*rho**(-1)
        
    elif scattermodel == 'fit8': # fit from 13/10-97
        gs6 = 136 * (pci**2.85) * ((freq/50)**2.5) / (rho + 0.001)
  
    elif scattermodel == 'fit9': # fit from 4/11-97 (without density)
        gs6 = 564 * (pci**3.0)* ((freq/50)**2.5)
  
    elif scattermodel == 'fit10': # fit from 4/11-97 (without density, uses corr. length from exp. fit) 
        gs6 = (3.16*pci + 295*(pci**2.5)) * (freq/50)**2.5  # muligvis deq. 60 i MEMLS paper af Wiesmann, 1998
        
    elif scattermodel == 'fit11': # fit from 4/11-97 (with density, uses corr. length from exp. fit!)
        gs6 = (9.20*pci - 1.23*rho + 0.54)**2.5 * (freq/50)**2.5  # do, eq. 59
        # valid for densities up to 0.4 g/cm3, pc from 0.05-0.3mm, freq from 5-100 GHz. 
        # for pcc < 0.05 mm, scattermodel 10 should be used. 
    return gs6

def improvedbornscattering(epsi_host,epsii_host,epsi_inc,epsii_inc,epsi_eff,epsii_eff,pci,vi,k,freq,method):
    # pci should here be used in m for the calculations    
    if method=='spheres':
        gs6 = (3./32) * pci**3 * k**4 * vi*(1-vi) * abs((2*epsi_eff+epsi_host)*(epsi_inc-epsi_host)/(2*epsi_eff+epsi_inc))**2  # eq. 33, Matzler (1998)
    elif method == 'shells':
        gs6 = abs(2./3 + 1./(3.*epsi_inc**2))*pci*k**2*vi*(1-vi)*(epsi_inc-1)**2/(16*epsi_eff)
    return gs6
    

def firstyearice(T,sal,pcc,k,freq):
    # Scattering coefficient of ice calculated from structural parameters
    
    # OBS: RECALCULATING epsi, epsii....
    [epsi_ice, epsii_ice] = permittivity.epureice(T,freq) 
    [epsi_brine, epsii_brine] = permittivity.ebrine(T,freq)
    
    # Improved born approximation by C. Mätzler (1998). J. Appl. Phys. 83(11),6111-7
    # Polder/VanSanten mixing formulae for spheical inclusions
    # effective dielectric constant of medium consisting of e1 and e2
    v = brine.volume(T,sal)
    # Assuming spherical brine inclusions:
    epsi_eff, epsii_eff = permittivity.sphericalinclusions(epsi_ice,epsii_ice,epsi_brine,epsii_brine,v,'ImprovedBorn')
        
    # Improved born approximation by C. M�tzler (1998). J. Appl. Phys. 83(11),6111-7
    # Scattering coefficient of a collection of spherical particles with 
    # correlation length pcc using improved born approximation:
    epsi_host = epsi_ice
    epsii_host = epsii_ice
    epsi_inc = epsi_brine
    epsii_inc = epsii_brine
    gs6 = improvedbornscattering(epsi_host,epsii_host,epsi_inc,epsii_inc,epsi_eff,epsii_eff,pcc*0.001,v,k,freq,'spheres')
    return gs6

def multiyearice(rho,pci,epsi,epsii,k,freq):
    epsi_ice = epsi # permittivity of the saline ice - and with air bubbles!!!
    epsii_ice = epsii 
    
    # Air volume:    
    volair=(0.926-rho)/0.926
    volair[volair<0] = 0.
    
    # Born approximation: mixing in air bubbles with the saline ice: already done??
    epsi_air = 1.0
    epsii_air = 0.0
    epsi_eff, epsii_eff = permittivity.sphericalinclusions(epsi_ice,epsii_ice,epsi_air,epsii_air,volair) 
    # Permittivity - er det ikke allerede udregnet??? Has already been done once?
    
    # scattering coefficient of a collection of spherical particles with 
    # correlation length pcc using improved born approximation
    epsi_host = epsi_ice
    epsii_host = epsii_ice
    epsi_inc = epsi_air
    epsii_inc = epsii_air
    gs6=improvedbornscattering(epsi_host,epsii_host,epsi_inc,epsii_inc,epsi_eff,epsii_eff,pci*0.001,volair,k,freq,'spheres') 
    return gs6

def refractionconstants(epsi,epsii,freq): # is not actually used
    # The complex index of refraction, n:
    ni = (np.sqrt(epsi-1j*epsii)).real
    nii = np.abs((np.sqrt(epsi-1j*epsii)).imag)
    n = ni-1j*nii

    # Absorption coefficient:    
    c = 2.99793*10**8 # Speed of light in vacuum [m/s]
    lambda0=c/(freq*10**9) # Free space wavelength [m]
    alpha = 2*np.pi/lambda0*nii
    # Phase constant:
    beta = 2*np.pi/lambda0*ni
    # Propagation constant:
    gamma = alpha + 1j*beta
    
    return n, gamma

def absorption(epsi,epsii,freq):
    """
    Computes the absorption coefficient for snow from the dielectric properties.
    [gai] = absorption(epsi,epsii,freq)
        gai:   absorption coefficient [m^-1]
        eps:   real part of dielectric permittivity
        epsi:  imaginary part of dielectric permittivity
        freq:  frequency [GHz]
        
        This is the 6-flux absorption coefficients, gamma_a
    """ 
    
    # Absorption coefficient suitable for snow.
    # From Ulaby et al. 1981 vol. 1 E4.140b
    
    lambda0=2.99793/(freq*10) # Free space wavelength [m]
    #gai=(4.0*np.pi/lambda0)*(np.sqrt(epsi+epsii*1j)).imag
    k=2.0*np.pi/lambda0
    gai=k*np.abs((np.sqrt(epsi-epsii*1j)).imag)
   
#    ni = (np.sqrt(epsi-1j*epsii)).real
#    nii = np.abs((np.sqrt(epsi-1j*epsii)).imag)
#    alpha = 2*np.pi/lambda0*nii
    gai = 2*gai #alpha
    
    # Not suitable for saline ice > about 10psu
    #gai = ((2*pi*10*freq).*epsii)./(c.*sqrt(epsi - (epsii.^2./4.*epsi)))
    return gai

def twofluxcomponents(gs6,epsi,ga6):
    """
    The 2-flux scattering and absorption coefficients calculated from 6-flux 
    components. Assuming no specular scattering.
    Equations from Matzler note 2. 
    
    [gb2h,gb2v,ga2] = twofluxcomponents(gs6,epsi,ga6)
        gb2h:  2-flux scattering coefficient at h polarization
        gb2v:  2-flux scattering coefficient at v polarization
        ga2:   2-flux absorption coefficient
        gs6:   6-flux scattering coefficient
        ga6:   6-flux absorption coefficient
        epsi:  relative permittivity of material
    """
    
    # Specular component of scattering coefficient.
    # Usually 0, but can be important in new snow. 
    dgb0h = 0
    dgb0v = 0
    
    # Convert to 2-flux absorption and scattering coefficients:
    # This is done according to the equations in Matzler note 2, eq. 12+13
    x = np.sqrt((epsi - 1)/epsi)
    gb6 = 0.5 * gs6 * (1-x)
    gc6 = 0.25 * gs6 * x # Coefficient for coupling between horizontal and vertical fluxes

    # The two-flux absorption coefficient: 
    gtr = (4*gc6) / (ga6 + 2*gc6) 
    ga2 = ga6 * (1 + gtr) # Matzler, note 2, eq. 12
    
    # The two-flux scattering coefficients:
    gb2h = (gb6 + dgb0h) + gtr*gc6 # Matzler, note 2, eq. 13, and including specular scattering
    gb2v = (gb6 + dgb0v) + gtr*gc6
    return gb2h, gb2v, ga2