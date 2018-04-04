#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Collection of subroutines for calculation of reflection and transmissivity 
coefficients for the snow and ice stack, while accounting for diffuse 
scattering of beam. Penetration depth is also calculated. 

Based off the original MEMLS code (1997) by the Institute of Applied Physics, 
University of Bern, Switzerland.

@author: Mai Winstrup, Danish Meteorological Institute, Copenhagen, Denmark, 2018
"""

import numpy as np
import matplotlib.pyplot as plt

def reflectivity(depth,layerthickness,theta_in,epsi,epsii,gs6,gbih,gbiv,gai,s0h,s0v):
    """
    [rih, riv, tih, tiv, sih, siv, theta_eff] = reflectivity(layerthickness,theta,epsi,epsii,gs6,gbih,gbiv,gai,s0h,s0v)
        rih, riv:       layer reflectivity, h/v polarization
        tih, tiv:       layer transmissivity, h/v polarization
        sih, siv:       layer interface reflectivity, h/v polarization
        theta_eff:      effective propagation angle [rad]
        layerthickness: layer thickness [m]
        theta_in:       beam incidence angle [rad]
        epsi, epsii:    permittivity, real and imaginary parts
        gs6:            6-flux scattering coefficients
        gbih, gbiv:     2-flux scattering coefficients
        gai:            2-flux absorption coefficients
        s0h, s0v:       ground-snow interface reflectivity (with ground = water)
        
    Uses: fresnelreflectivity, effpropagationangle, polarizationmixing, 
          reflectivitytransmissivity, penetrationdepth
    """
    
    # Reflectivities from fresnel coefficients and assuming Snell's 
    # approximation for propagation angle (=theta_snell):
    [sih,siv,theta_snell] = fresnelreflectivity(theta_in,epsi,s0h,s0v)
    
    # Effective propagation angle, path length, and one-way transmittivity, 
    # when accounting for diffuse scattering:
    theta_eff,pathlength_eff,tscat,pathlength = effpropagationangle(theta_snell,epsi,gs6,layerthickness)
    
    # Polarization mixing of interface reflectivities (calculated from 
    # first-order reflectivities):
    sih,siv = polarizationmixing(sih,siv,tscat)
    
    # Improved layer reflectivity and transmissivity:
#    [rih,tih] = reflectivitytransmissivity(gai,gbih,pathlength_eff)
#    [riv,tiv] = reflectivitytransmissivity(gai,gbiv,pathlength_eff)
    [rih,tih] = reflectivitytransmissivity(gai,gbih,pathlength)
    [riv,tiv] = reflectivitytransmissivity(gai,gbiv,pathlength)
    # We here use the updated pathlength, whereas the original MEMLS files 
    # use the Snell-derived pathlength for this computation. 
    
    # Penetration depth:
    dpenh = penetrationdepth(gs6,gai,tih,pathlength,layerthickness)
    dpenv = penetrationdepth(gs6,gai,tiv,pathlength,layerthickness)
    
    return rih, riv, tih, tiv, sih, siv, theta_eff, dpenh, dpenv

def fresnelreflectivity(theta,epsi,s0h,s0v):
    """ 
    Fresnel reflection coefficients, when assuming a loss-less medium (i.e. 
    eps'' = 0). Layer n+1 is the air above the snowpack. 

    [sih,siv,theta_snell] = fresnelreflectivity(theta,epsi,s0h,s0v)
        sih, siv:    interface reflectivity, h/v polarization (incl. ground layer)
        theta_snell: Snell propagation angle through snowpack (incl. air layer)
        theta:       local incidence angle
        epsi:        real part of dielectric permittivity
        s0h, s0v:    ground reflectivity, h/v-polarization
    """
    
    # Propagation angle with respect to the surface normal, when assuming the 
    # propagation through the snowpack is governed by Snells law, both for 
    # vertically and horizontally polarized radiation.
    ns = np.sqrt(epsi) # Refraction index
    theta_snell = np.append(np.arcsin(np.sin(theta)/ns), theta)
    # Appended with ray incidence angle for layer n+1 (the air above the snowpack). 

    # This is Snell's approximation for propagation angle in layer j.
    # Deviations in propagation angle occur due to volume scattering, as 
    # accounted for later. 
    
    # Fresnel coefficients:
    N = len(epsi) # Number of snow/ice layers
    epsi = np.append(epsi,1) # Using air as layer n+1
    depsn = epsi[0:N]/epsi[1:N+1]
    tein = theta_snell[1:N+1]
    part1 = np.sqrt(depsn-(np.sin(tein))**2)
    frh = (np.cos(tein)-part1)/(np.cos(tein)+part1)
    frv = (part1-depsn*np.cos(tein))/(part1+depsn*np.cos(tein))
    
    # Interface reflectivity between two regular layers is the square of the 
    # fresnel reflection coefficients:
    sih = frh**2
    siv = frv**2
    
    # Append with "ground" reflectivities:
    # Ground is here the ice-water interface.
    sih = np.append(s0h,sih) 
    siv = np.append(s0v,siv)
    
    return sih, siv, theta_snell

def effpropagationangle(theta_snell,epsi,gs6,layerthickness):
    """
    Effective propagation angle and effective path length after correction for 
    diffuse scattering. Also the total one-way transmittivity from snow surface 
    to bottom of layer (tscat) is calculated. 
    Equations from Wiessmann & Matzler (1999), p.9-10, E38-41. 
    
    theta_eff, pathlength_eff, tscat = effpropagationangle(theta_snell,epsi,gs6,layerthickness)        
        theta_eff:      Effective propagation angle in each layer after correction [rad]
        pathlength_eff: Effective pathlength within each layer after correction [m]
        tscat:          Total one-way transmittivity from snow surface to bottom of layer
        theta_snell:    Snells propagation angle in layer [rad]
        epsi:           Relative permittivity
        gs6:            6-flux scattering coefficient
        layerthickness: Layer thickness [m]
    """

    # Maximum propagation angle:
    ns = np.sqrt(epsi) # Refraction index, incl air layer
    costheta_max = 0.5*(1+np.sqrt(1-(1/ns)**2))
    # Actual value is between max and Snell values for propagation angle.
    # Weighting according to tscat; the transmissivity for scattering radiation 
    # from snow surface to layer surface. 
    
    # Transmittivity from surface and down the snow/ice stack for 
    # non-scattering radiation:
    # Pathlength through layer, using snell approximation for angle:
    pathlength = layerthickness/np.cos(theta_snell[:-1]) # Not including air
    dscat = 0.5*gs6*pathlength # Value for each layer
    # Summing from surface and down the snow/ice stack:
    tscat = np.exp(-np.cumsum(dscat[::-1])) 
    # And reverse to usual direction:
    tscat = tscat[::-1] # One way transmittivity, to the bottom of each layer

    # Effective propagation angle:
    costheta_eff = tscat*np.cos(theta_snell[:-1])+(1-tscat)*costheta_max
    theta_eff = np.arccos(costheta_eff)

    # Updated value for pathlength:
    pathlength_eff = layerthickness/np.cos(theta_eff)
    
    return theta_eff, pathlength_eff, tscat, pathlength

def polarizationmixing(sih,siv,tscat):
    """
    Polarization mixing of the interface reflectivities of each layer due to 
    diffuse scattering of the photons (calculated from the first-order 
    reflectivities). Wiesmann & Matzler (1999), E42-44.

    [sih_adj,siv_adj] = polarizationmixing(sih,siv,tscat)
        sih_eff: effective interface reflectivity at h-polarization after adjustments for mixing
        siv_eff: effective interface reflectivity at v-polarization after adjustments for mixing
        sih:     interface reflectivity at h-polarization
        siv:     interface reflectivity at v-polarization
        tscat:   total transmissivity from snow surface to bottom of layer
    """
    # Assuming negligible scattering in the air, which corresponds to tscat = 1
    # (tscat is calculated for lower interface of layer).
    tscat = np.append(tscat,1) 
    
    smean = 0.5*(sih+siv)
    deltas = 0.5*tscat*(sih-siv)
    sih_eff = smean + deltas
    siv_eff = smean - deltas
    
    return sih_eff,siv_eff

def reflectivitytransmissivity(gai,gbi,pathlength):
    """
    Layer reflectivity and transmissivity. 
    From Wiesmann & Matzler (1999), E22-26.
    
    [ri,ti] = reflectivitytransmissivity(gai,gbi,pathlength)
        ri:   layer reflectivity
        ti:   layer transmissivity
        gai:  absorption coefficient
        gbi:  scattering coefficient 
        pathlength:  effective path length    
    """ 
   
    gamma = np.sqrt(gai * (gai + 2*gbi))
    t0 = np.exp(-gamma*pathlength)
    
    r0 = np.zeros(len(t0))
    mask = gbi>0.00001 
    r0[mask] = gbi[mask] / (gai[mask]+gbi[mask]+gamma[mask])
    
    t02 = t0**2
    r02 = r0**2
    ri = r0*(1-t02)/(1-t02*r02)
    ti = t0*(1-r02)/(1-t02*r02)
    
    return ri, ti

def penetrationdepth(gs6,gai,ti,pathlength,layerthickness):
    """
    Beam penetration depth, i.e. the depth where beam is attenuated to 1/e. 
    
    [dpene] = penetrationdepth(layerthickness,ke,theta,trans)
    
    dpene:          penetration depth [m]
    layerthickness: layer thicknesses [m]
    theta:          propagation angle in snow pack [rad]
    
    
    #L=exp(ke*d*sec(theta))
    #sec = 1/cos
    # num : layers
    # di : layer thickness
    # ke : gs6+ga2i (?) - extinction coefficient
    # theta : rtei (?) - propagation angle
    # trans : ti (?)
    """
    
    # Extinction coefficient for each layer (with contribution from absorption 
    # and scattering):
    ge = gs6+gai
    # Fraction of incoming radiation lost during traverse through each layer:
    lossfactor=np.exp(-ge*pathlength) 
    
   # dpene_ice = 1/ge[0]
   # dpene_snow = 1/ge[-1]
   # print 'ice est:', dpene_ice 
   # print 'snow est:', dpene_snow
    
    # Attenuation of beam at bottom of each layer:
    # Summing contributions from top and down the snow/ice stack:
    tifromtop = np.append(1,ti[::-1])
    totalloss = np.cumsum(tifromtop[0:-1]*lossfactor[::-1])  # må være korrekt
    #totalloss = np.cumsum(ti[::-1]*lossfactor[::-1]) # clara
    
    # Linear interpolation between layers to obtain the depth where beam is 
    # attenuated to 1/e (=0.37) of incoming radiation:
    depth = np.cumprod(layerthickness[::-1]) # Depth at bottom of layer
    pdepth = np.interp(1.0/np.e,np.append(0,totalloss),np.append(0,depth))
    
    return pdepth
    
    # Total loss:
#    totalloss = np.cumsum(loss)
#    np.prod(1.0/loss[idx:nn])*np.prod(trans[idx:nn])
    
#    num = range(0,len(di))
    
#    nn = len(num)
 #   inten=np.ones(len(num))
  #  pene=np.zeros(len(num)) 
    
   # einv=0.367879 #1/e
    #x=0.0
    
    #loss=np.exp(ke*di*(1/np.cos(theta)))
    
   # trans=np.array(trans) #probably ti
#    #from Rasmus
#    for idx in range(nn):
#      inten[idx]=np.prod(1.0/loss[0:idx])*np.prod(trans[0:idx])
#      pene[idx]=np.sum(di[0:idx])
    #modified by Clara
    #for idx in range(len(num)-1,-1,-1):
     # #print idx,'/',nn
      #inten[idx]=np.prod(1.0/loss[idx:nn])*np.prod(trans[idx:nn]) #looks at the radiation
      #pene[idx]=sum(di[idx:nn]) #depth of the interface where 1/e is reached
      #print inten
      
      # looks more precisely into the layer, where this is (by assuming a linear profile)
      #if inten[idx] <= einv:
          #print inten[idx]
          #print idx,'/',nn
          #print pene[idx],'/',sum(di)
       #   if idx==nn-1:
#              #print 'loop 1'
#              slope=(1.0-inten[idx])/di[idx]
#              x=(1.0-einv)/slope
#              #print x
#              break
#          else:
#              print 'loop 2'
#              slope=(inten[idx]-inten[idx+1])/di[idx]
#              x=pene[idx+1]+(inten[idx+1]-einv)/slope
#              #print x
#              break
#    return x
#
#### by Clara, using old deffinition (without scattering loss)
#def pendep(num,di,abscoeff):
#    # num : layers
#    # di : layer thickness
#    # abscoeff : ga2i
#    nn = len(num)
#    dep = np.zeros(len(num))
#    pene=np.zeros(len(num))
#    einv=0.367879 #1/e
#    x=0.0
#    for idx in range(len(num)-1,-1,-1):
#        print idx
#        print inten
#        dep[idx] = 1./abscoeff[idx]
#        pene[idx]=sum(di[idx:nn]) #depth of the interface where 1/e is reached
#        # looks more precisely into the layer, where this is (by assuming a linear profile)
##        if inten[idx] <= einv: 
##            if idx==nn:
##                slope=(1.0-inten[idx])/di[idx]
##                x=(1.0-einv)/slope
##                break
##            else:
##                slope=(inten[idx-1]-inten[idx])/di[idx]
##                x=pene[idx-1]+(inten[idx-1]-einv)/slope
##                break
#    return pene  