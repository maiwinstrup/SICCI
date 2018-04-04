#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 17:44:54 2017

Functions needed for MEMLS itself

@author: m300411
"""

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

def epsice(Ti,freq):
#%   calculates the dielectric loss factor of pure ice 
#%   After Hufford, Mitzima and Matzler
#%
#%   eice = epsice(Ti,freq)
#%      eice:  dielectric permittivity of ice (imaginary part)
#%      Ti:    temperature in K 
#%      freq:  frequency in GHz
#%
#%   Version history:
#%      1.0    wi 15.7.95
#%   
#%   Uses:
#%
#%
#%
#%   Copyright (c) 1997 by the Institute of Applied Physics,
#%   University of Bern, Switzerland
  pp = (300.0/Ti)-1;
  B1 = 0.0207
  b = 335.25
  B2 = 1.16e-11
  db = np.exp(-10.02+0.0364*(Ti-273.15)) #replaced 273 by 273.15 by Clara
  beta = ((B1*np.exp(b/Ti))/(Ti*(np.exp(b/Ti)-1)**2))+B2*freq**2+db
  alpha = ((0.00504 + 0.0062*pp)* np.exp(-22.1*pp))
  eice = (alpha/freq) + (beta*freq)
  return eice
  
def epsr(roi):
#%   calculates the dielectric permittivity for dry snow from 
#%   density .
#%
#%   epsi = epsr(roi)
#%       epsi:  real part of dielectric permittivity
#%       roi:   density g/cm^3
#%
#%   Version history:
#%      1.0    wi 15.7.95
#%      1.1    wi 23.9.97 added Looyenga for snow denser than 0.4 g/cm^3
#%
#%   Uses:
#%       epsice
#%
#%
#%
#%   Copyright (c) 1997 by the Institute of Applied Physics,
#%   University of Bern, Switzerland

  epsi = np.zeros(len(roi))
  ilist = np.where(roi<=0.4)
  jlist = np.where(roi>0.4)
  vfi = roi/0.917
  ehb = 0.99913
  esb = 1.4759
  
  for i in  ilist:
    epsi[i] = 1 + 1.5995 * roi[i] + 1.861 * roi[i]**3
  for j in jlist:
    epsi[j] = ((1 - vfi[j]) * ehb + vfi[j] * esb)**3
  return epsi
  
  
def ro2epsd(roi,Ti,freq):
#%   calculates the dielectric permittivity from 
#%   density for dry snow.
#%
#%   [epsi,epsii] = ro2epsd(roi,Ti,freq)
#%       epsi:  real part of dielectric permittivity
#%       epsii: imaginary part of dielectric permittivity
#%       roi:   density
#%       Ti:    snow temperature in Kelvin
#%       freq:  frequency
#%
#%   Version history:
#%      1.0    wi 15.7.95
#%      2.0    wi 12.11.97  enhanced with Polder and van Santen Equations (see Polder.m)
#%
#%   Uses:
#%       epsice, epsr, polder
#%
#%
#%
#%   Copyright (c) 1997 by the Institute of Applied Physics, 
#%   University of Bern, Switzerland

  eice = epsice(Ti,freq)
  epsi = epsr(roi)
  
  #imaginary part after Tiuri 84
  #epsii = eice.*(0.52.*roi + 0.62.*(roi.^2));
  
  #imaginary part after Polder and van Santen 1946 (Effective-Medium Approx)
  f = roi/0.917;
  ei = 3.185;
  N = len(roi)
  A = np.zeros((N))
  epsp = np.zeros((N))
  A = A + 0.3;
  for i in range(N):
     if f[i] < 0.55:
        A[i] = 0.476 - 0.64 * f[i]
     if f[i] <= 0.333:
        A[i] = 0.1 + 0.5 * f[i]

     #epsp(i) = polder(A(i),ei,f(i));

  epsp = epsi
  A3 = 1-2*A
  ea = (epsp*(1-A))+A;
  ea3 = epsp*(1-A3)+A3;
  K1 = (ea/(ea+A*(ei-1)))**2
  K3 = (ea3/(ea3+A3*(ei-1)))**2;
  Ksq = (2*K1+K3)/3
  epsii = np.sqrt(epsi)*eice*Ksq*f 
  return epsi, epsii

def mixmod(f,Ti,Wi,epsi,epsii):  
#%   calculates the permittivity for Wetness > 0
#%      Physical Mixing Model Weise 97 after Matzler 1987 (corrected)
#%      water temperature is assumed constant at 273.15 K
#%
#%   [epsi,epsii] = mixmod(f,Ti,Wi,epsi,epsii)
#%       epsi:  real part of the permittivity
#%       epsii: imaginary part of the permittivity
#%       f:     frequency [GHz]
#%       Ti:    physical snow temperature
#%       Wi:    wetness [%], no! in vol. frac. 0-1. according to mail 06/01/05 Mätzler -rtt
#%       epsi:  real part of dry snow perm.
#%       epsii: imaginary part of dry snow perm.
#%
#%   Version history:
#%      1.0    wi 15.7.95
#%   
#%   Uses: -
#%
#%
#%   Copyright (c) 1997 by the Institute of Applied Physics,
#%   University of Bern, Switzerland

  Aa = 0.005
  Ab = 0.4975
  Ac = 0.4975
  euw = 4.9
  esw = 88.045 
  frw = 0.11109 # inverse relaxation frequency of water
  
  esa = (esw - epsi)/(3*(1+Aa*(esw/epsi-1)))
  esb = (esw - epsi)/(3*(1+Ab*(esw/epsi-1)))
  esc = (esw - epsi)/(3*(1+Ac*(esw/epsi-1)))
  eua = (euw - epsi)/(3*(1+Aa*(euw/epsi-1)))
  eub = (euw - epsi)/(3*(1+Ab*(euw/epsi-1)))
  euc = (euw - epsi)/(3*(1+Ac*(euw/epsi-1)))
  
  fa = 1 + Aa * (esw-euw)/(epsi+Aa*(euw-epsi))
  fb = 1 + Ab * (esw-euw)/(epsi+Ab*(euw-epsi))
  fc = 1 + Ac * (esw-euw)/(epsi+Ac*(euw-epsi))
  
  eea = esa - eua
  eeb = esb - eub
  eec = esc - euc
  
  fwa = frw/fa
  fwb = frw/fb
  fwc = frw/fc
  
  depsia = eua + eea / (1+(fwa*f)**2)
  depsib = eub + eeb / (1+(fwb*f)**2)
  depsic = euc + eec / (1+(fwc*f)**2)
  depsi = Wi * (depsia+depsib+depsic)
  
  depsiia = fwa*f*eea / (1+(fwa*f)**2);
  depsiib = fwb*f*eeb / (1+(fwb*f)**2);
  depsiic = fwc*f*eec / (1+(fwc*f)**2);
  depsii = Wi * (depsiia + depsiib + depsiic);
  
  epsi = epsi + depsi
  epsii = epsii + depsii
  return epsi, epsii

def epice(T,freq):
#% Dielectric constant of pure ice: Mätzler, 1998, Microwave properties of ice and snow, in:
#% B. Smitht et al. (eds.) solar system ices, 241-257, Kluwer.
#
#% T: thermometric temperature in K
#% freq: Frequency in GHz
  if np.any(T < 100):
      T = T+273.15
  epi=3.1884+9.1e-4*(T-273.15) #modified from 273 to 273.15 
  
  #The Hufford model for the imagenary part.
  
#  theta=300/T
#  alpha=(0.00504+0.0062*theta)*np.exp(-22.1*theta)
#  beta=(0.502-0.131*theta/(1+theta))*1e-4 +(0.542e-6*((1+theta)/(theta+0.0073))**2)
#  epii=(alpha/freq) + (beta*freq)
  Ti = T
  pp = (300.0/Ti)-1;
  B1 = 0.0207
  b = 335.25
  B2 = 1.16e-11
  db = np.exp(-10.02+0.0364*(Ti-273.15)) #replaced 273 by 273.15 by Clara
  beta = ((B1*np.exp(b/Ti))/(Ti*(np.exp(b/Ti)-1)**2))+B2*freq**2+db
  alpha = ((0.00504 + 0.0062*pp)* np.exp(-22.1*pp))
  epii = (alpha/freq) + (beta*freq)
    
  epui=epi+epii*1j
  return epui

def Sb(T):
  if np.any(T>100.):
    T = T-273.15
  Sb = np.ones(len(T))*np.nan  
  for n,Ti in enumerate(T):
      if Ti >= -43.2 and Ti <= -36.8:
          Sb[n]=508.18 + 14.535*Ti + 0.2018*Ti**2
      elif Ti>=-36.8 and Ti<=-22.9:
          Sb[n]= 242.94 + 1.5299*Ti + 0.0429*Ti**2
      elif Ti>-22.9 and Ti<-8:
          Sb[n] = -1.20 - 21.8*Ti - 0.919*Ti**2 - 0.0178*Ti**3
      elif Ti>=-8 and Ti<0.:
          Sb[n] = 1./(0.001-(0.05411/Ti)) 
      elif Ti>=0:
          Sb[n] = 0.
  return Sb
  #return T

def Vb(T,S):
    #liquid water fraction or volume brine fraction
    #if everything is nan
    if np.any(T[~np.isnan(T)])==False:
      lwf = T
    else:
      if np.any(T>100.):
        T = T-273.15
      Sbr = Sb(T)
      lwf=np.ones(len(S))*np.nan
      for n in range(len(S)):
          if Sbr[n]>0.:
              lwf[n] = S[n]/Sbr[n]
          elif Sbr[n]<=0:
              lwf[n]=1.
      lwf[lwf>1.] = 1
    return lwf


def Nsw(Ssw):
#% normality of sea water or brine Ulaby et al. 1986, E20
#% Ssw: salinity of brine or sea water
  N = 0.9141*Ssw*(1.707e-2 +1.205e-5*Ssw+4.058e-9*(Ssw**2))
  return N

def condbrine(T):
  #conductivity of brine  
  #from Stogryn and Desargant, 1985, Eq.7
  if any(T>100):
      T=T-273.15
  condbrine=np.zeros(len(T))    
  for i,Ti in enumerate(T):
      if Ti >= -22.9:
          condbrine[i] = -Ti*np.exp(0.5193 + 0.8755*0.1*Ti)
      elif Ti < -22.9:
          condbrine[i] = -Ti*np.exp(1.0334 + 0.1100*Ti)
  return condbrine        
          
def relaxt(T):
  #relaxation time  
  #from Stogryn and Desargant, 1985, Eq.12
  #fit up to -25°C
  #in nanoseconds  
  # not really close to the MEMLS version
  if any(T>100):
      T=T-273.15
  relax = (1./(2*np.pi))*(0.10990 + 0.13603*0.01*T + 0.20894*0.001*T**2 + 0.28167*0.00001*T**3)
  relax[T==0] = (1./(2*np.pi))*0.1121
  relax = relax/1e9
  return relax

def epsib0(T):
  #static dielectric const. of brin
  #from Stogryn and Desargant, 1985, Eq.10
  #fit up to -25°C
  if any(T>100):
      T=T-273.15    
  epsib0 = (939.66 - 19.068*T)/(10.737 - T)
  epsib0[T==0] = 87.92
  return epsib0

def ebrine(T,freq):
    # brine permittivity
    #from Stogryn and Desargant, 1985, Eq.1
    #fit up to -25°C
    if any(T>100):
      T=T-273.15   
      
    f = freq*1e9
    e0 = 8.85419*1e-12
    
    epsiwoo = (82.79 + 8.19*T**2)/(15.68 + T**2)
    epsiwoo[T==0] = 5.28
    
    epsis = epsib0(T)
    tau = relaxt(T)
    sig = condbrine(T)
    
    ebr = epsiwoo + (epsis - epsiwoo)/(1.-2*np.pi*f*tau*1j) + (1j*sig)/(2*np.pi*e0*f)
    
    return ebr.real, ebr.imag


def eice_s2p(e1,e2,v):
#% improved born approximation by C. Mätzler (1998). J. Appl. Phys. 83(11),6111-7
#% Polder/VanSanten mixing formulae for spheical inclusions
#%effective dielectric constant of medium consisting of e1 and e2
#%e1: dielectric constant of background
#%e2: dielectric constant of sherical inclusion
#%v: fraction of inclusions

  eeff=0.25*(2.*e1-e2+3.*v*(e2-e1)+np.sqrt((2.*e1-e2+3.*v*(e2-e1))**2 +8.*e1*e2))
  return eeff
  
  
def sie(si,sal,Ti,freq,epsi,epsii):  
  #% computes the dielectric constant of ice if it is an ice layer
  #% Background information: Ulaby et al. 1986 vol. III.
  #% i.e. if si == 1
  #% si: sea ice/snow [1/0]
  #% sal: salinity [ppt/psu]
  #% Ti: Thermometric temperature of layer [K]
  #% freq: Frequency in GHz
  #% epsi: initial permittivity (of snow)
  #% epsii: initial loss (of snow)
  
  T=Ti-273.15                          #get therm. temp. in C
  eice=epice(Ti,freq)                  #fresh ice dielectric const
  volb=Vb(T,sal)                       #volume of brine                        #salinity of brine                     
  [eb,ebi]=ebrine(T,freq)             #dielectric constant of brine
  #emis=eice_rn2p(eice,eb+ebi*i,volb)  #dielectric constant of sea ice (random needles)
  emis=eice_s2p(eice,eb+ebi*1j,volb)    #dielectric constant of sea ice (spherical inclusions)
  
  aepsi=emis.real
  aepsii=emis.imag
  ###added by Clara
  #aepsii[aepsii>1.] = 1.
  ###
  
  epsi=epsi-epsi*si+aepsi*si
  epsii=epsii-epsii*si+aepsii*si
  return epsi, epsii

def mysie(si,rho,Ti,sal,freq,epsi,epsii):
    #% computes the dielectric constant of ice if it is an ice layer
    #% Background information: Ulaby et al. 1986 vol. III.
    #% i.e. if si == 1
    #% si: sea ice/snow [1/0]
    #% rho: density of icelayer [kg/m3]
    #% Ti: Thermometric temperature of layer [K]
    #% sal: salinity of ice [psu]
    #% freq: Frequency in GHz
    #% epsi: initial permittivity (of snow)
    #% epsii: initial loss (of snow)

    #permittivity of saline ice
    [sepsi, sepsii] = sie(si,sal,Ti,freq,epsi,epsii)

    eice=epice(Ti,freq)                            #fresh ice dielectric const
    vola=(0.926 - rho)/0.926                      #volume of air
    ###modified by Clara
    vola[vola<0]=0.   
    ###
    emis=eice_s2p(sepsi + 1j*sepsii,1.0 + 0.0j,vola)   #dielectric constant of sea ice (spherical inclusions)

    aepsi=emis.real
    aepsii=emis.imag

    epsi=epsi-epsi*si+aepsi*si
    epsii=epsii-epsii*si+aepsii*si
    
    ##added by Clara
    #epsii[epsii>1] = 1.
    ###
    
    return epsi, epsii
  
def abscoeff(epsi,epsii,Ti,freq,Wi):
#  %   computes the absorption coefficient from the dielectric properties
#  %
#  %
#  %   [gai] = abscoeff(epsi,epsii,Ti,freq,Wi)
#  %       gai:   absorption coefficient [m^-1]
#  %       epsi:  real part diel
#  %       epsii: imaginary part diel
#  %       Ti:    physical temperature
#  %       freq:  frequency [GHz]
#  %       Wi:    volumetric liquid water content
#  %
#  %   Version history:
#  %      1.0    wi 15.7.95
#  %      1.1    wi 12.11.97 more precise formula for gai used
#  %   
#  %   Uses:
#  %
#  %
#  %
#  %   Copyright (c) 1997 by the Institute of Applied Physics, 
#  %   University of Bern, Switzerland

  # constants
  c = 2.99793

  lamd=c/(10*freq)
  gai=(4.0*np.pi/lamd)*(np.sqrt(epsi+epsii*1j)).imag
  # Absorption coefficient, suitable for snow but not saline ice > about 10psu
  #gai = ((2*pi*10*freq).*epsii)./(c.*sqrt(epsi - (epsii.^2./4.*epsi)));
  return gai
  
def pfadi(tei,di):
#%   calculates the effective path length in a layer
#%
#%   dei = pfadi(tei,di)
#%       dei:  effective path length [m]
#%       tei:  local incidence angle
#%       di:   thickness [m]
#%
#%   Version history:
#%      1.0    wi 15.7.95
#%   
#%   Uses:
#%
#%
#%
#%   Copyright (c) 1997 by the Institute of Applied Physics,
#%   University of Bern, Switzerland
  N = len(di)
  dei = di/np.cos(tei[0:N])
  return dei
  
def fresnelc0(tei,epsi):
#%   fresnel reflection coefficients (assuming eps'' = 0)
#%     (layer n+1 is the air above the snowpack)
#%
#%   [sih,siv] = fresnel(tei,roi)
#%       sih:  interface reflectivity at h pol
#%       siv:  interface reflectivity at v pol
#%       tei:  local incidence angle
#%       epsi: real part of dielectric permittivity
#%
#%   Version history:
#%      1.0    wi 15.7.97
#%   
#%   Uses:
#%       epsr
#%
#%
#%   Copyright (c) 1997 by the Institute of Applied Physics,
#%   University of Bern, Switzerland
  N = len(epsi)-1;
  siv = np.zeros(len(epsi[0:N]))
  sih = np.zeros(len(epsi[0:N]))
  
  for n in range(N):
    epso = epsi[n+1]
    epsu = epsi[n]
    #tein = tei[n] #replaced n+1 by n by Clara # commented on 22.01.2018 because results are too different to previous ones
    tein = tei[n+1]
    # these coefficients are the square of the Fresnel coefficients
    sih[n] = ((np.sqrt(epso)*np.cos(tein) - np.sqrt(epsu - epso * np.sin(tein)**2))/(np.sqrt(epso)*np.cos(tein) + np.sqrt(epsu - epso * np.sin(tein)**2)))**2
    siv[n] = ((epsu*np.cos(tein) - np.sqrt(epso)*np.sqrt(epsu - epso*np.sin(tein)**2))/(epsu*np.cos(tein) + np.sqrt(epso)*np.sqrt(epsu - epso*np.sin(tein)**2)))**2
  
  return sih, siv

def fresnelrc(tei,epsr):
#%   fresnel reflection coefficients (assuming eps'' = 0)
#%     (layer n+1 is the air above the snowpack)
#%
#%   [FH,FV] = fresnelrc(tei,epsr)
#%       FH:   Fresnel reflection coefficient at h pol
#%       FV:   Fresnel reflection coefficient at v pol
#%       tei:  local incidence angle
#%       epsr: (real part) dielectric permittivity
#%
#%   Version history:
#%      1.0    wi 15.7.95
#%   
#%   Uses:
#%
#%
#%
#%   Copyright (c) 1997 by the Institute of Applied Physics, 
#%   University of Bern, Switzerland
  N = len(epsr)
  
  FH = np.zeros(len(epsr[0:N-1]))
  FV = np.zeros(len(epsr[0:N-1]))
  
  #% epsi = [epsi;0];
  
  for n in range(0,N-1):
    epsn = epsr[n]/epsr[n+1]
    #tein = tei[n] #replaced n+1 by n by Clara  # commented on 22.01.2018 because results are too different to previous ones
    tein = tei[n+1]
    sinq = np.sin(tein)**2
    qeps = sinq/epsn
    wurz = np.sqrt(1-qeps)
    wsub = epsn-sinq
    nd = np.sqrt(epsn)
    
    FH[n] = ((nd*wurz-np.cos(tein))/(nd*wurz+np.cos(tein)))
    FV[n] = ((wurz-nd*np.cos(tein))/(wurz+nd*np.cos(tein)))
    
  return FH,FV  
  
def slred(num,roi,epsi,epsii,tei,sih,siv,di,dei,Ti,pci,freq,Wi,gai,si,sal):
#%   locates and treats coherent layers in a snowpack
#%     see Technote 11
#%     with repo=1 the snow layer table is printed after substeps
#%
#%   [num,roi,epsi,epsii,tei,sih,siv,di,dei,Ti,pci,Wi,gai,si,sal] = slred(num,roi,epsi,epsii,tei,sih,siv,di,dei,Ti,pci,freq,Wi,gai,si,sal)
#%       num:  index of the layer in the original snowpack
#%       roi:  density [g/cm^3]
#%       epsi:
#%       epsii:
#%       tei:  local incidence angle
#%       sih:  layer reflectivity at h pol
#%       siv:  layer reflectivity at v pol
#%       di:   layer thickness
#%       dei:  local path length [m]
#%       Ti:   physical snow temperature [K]
#%       pci:  correlation length [mm]
#%       Wi:   wetness  
#%       gai:  absorption coefficient
#%       freq: frequency
#%       si: sea ice layer (0/1)?
#%       sal: salinity [ppt]
#%
#%      1.0    wi 21.8.95
#%      2.0    wi 13.8.98 completely rewritten
#%   
#%   Uses:
#%       fresnelrc
#%
#%
#%
#%   Copyright (c) 1997 by the Institute of Applied Physics, 
#%   University of Bern, Switzerland
  #%  constants
  cc = 0.299793
  FIC = 4 * np.pi * freq / cc
  fc = 4.712
  repo = 0
  #   is there any layer -> checked in main
  N = len(roi)
  theta = tei[N]
  ns = np.sqrt(epsi)
  fi = FIC*di*ns*np.cos(tei[0:N])
    
  if repo == 1:
     print num,di,roi,Ti,pci,ns,gai,si,sal,epsi,epsii,fi,tei[0:N]*180/np.pi

  # %   find thin layers 
  #ilist = where(fi < fc)
  i = np.where(fi < fc)
  A = 0.* roi
  #for i in ilist:
  A[i] = A[i] + 1;
  #%   bottom layer is assumed noncoherent
  A[0] = 0
  
  if len(i[0]) > 0:
        #%  disp ([num2str(max(size(i))),' coherent layers of ',num2str(N),' detected: ',num2str(freq),'GHz']);
        # %   identify succeeding coherent layers and mark the packs from 2 ... scmax
        pl = 0
        sc = 1
        scmax = 0
        ml = 0
        mlo = 0
        for m in range(1,N):
          if A[m] == 1 and pl == 1:
              if ml == 0:
                sc = sc + 1
                ml = 1
                A[m-1] = sc
                # A[m] = 1, pl = 1
                A[m] = sc
                scmax = sc
          else: 
            if pl == 1:
              #% A(m) = 0, pl = 1 -> non coherent
              pl = 0
            else: 
              if A[m] == 1:
                  # A(m) = 1, pl = 0 -> first coherent
                  pl = 1
                  ml = 0

    
        #%  combine succeeding coherent layers by weighting with the phase
        if scmax > 0:
          for m in range(1,scmax):
            B = np.where(A == m+1)
            fitot = sum(fi[B])
            fitv = fi[B] / fitot
            tal = max(B[0])
            di[tal] = sum(di[B])
            dei[tal] = sum(dei[B])
            roi[tal] = sum(roi[B]* fitv)
            Ti[tal] = sum(Ti[B]* fitv)
            Wi[tal] = sum(Wi[B]* fitv)
            pci[tal] = sum(pci[B]* fitv)
            ns[tal] = sum(ns[B]* fitv)
            gai[tal] = sum(gai[B]* fitv)
            si[tal] = sum(si[B]* fitv)
            sal[tal] = sum(sal[B]* fitv)
            epsi[tal] = sum(epsi[B]* fitv)
            epsii[tal] = sum(epsii[B]* fitv)
            tei[tal] = sum(tei[B]* fitv)
            fi[tal] = fitot
            A[tal] = 1
        i = np.where(A < 2)
        if len(i[0]) > 0:
          num=num[i]
          roi=roi[i]
          tei=np.append(tei[i],tei[N])
          di=di[i]
          Ti=Ti[i]
          pci=pci[i]
          Wi=Wi[i]
          gai=gai[i]
          si=si[i]
          sal=sal[i]
          ns=ns[i]
          epsi=epsi[i]
          epsii=epsii[i]
          fi=fi[i]
          dei=dei[i]
          N = len(roi)
        if (repo == 1):
          print [num,di,roi,Ti,pci,ns,gai,si,sal,epsi,epsii,fi,tei[0:N]*180/np.pi]

        
        i = np.where(fi < fc)
        A = roi * 0
        A[i] = A[i] + 1
        A[0] = 0
        sih = np.zeros(len(roi))
        siv = np.zeros(len(roi))
        X   = np.zeros(len(roi))
        
        #%   calculate interface reflection coefficients
        [FH,FV] = fresnelrc(tei,np.append(epsi,1))
    
        #%   reduction on layers of type 0 (coherent layer effects are 
        #%     taken into account in the layer reflectivities)
        #%     for layers of type 0 shi = FH^2
        i = np.where(A == 0)
        #s are square of Fresnel coefficients
        sih[i] = FH[i]**2
        siv[i] = FV[i]**2
    
    
        #%     for layers of type 1 shi-1 = ...
        ilist = np.where(A == 1)
        for i in ilist:
          #Eq. 35 of the paper
          X[i] =  2*FH[i]*FH[i-1]*np.cos(fi[i])
          sih[i-1] = (FH[i]**2+FH[i-1]**2+X[i])/(1+FH[i]**2*FH[i-1]**2+X[i])
          X[i] =  2*FV[i]*FV[i-1]*np.cos(fi[i])
          siv[i-1] = (FV[i]**2+FV[i-1]**2+X[i])/(1+FV[i]**2*FV[i-1]**2+X[i])
          N = len(di)
        
        if (repo == 1):
          print num,di,roi,Ti,pci,ns,gai,si,sal,epsi,epsii,fi,tei[0:N]*180/np.pi,FH,FV,sih,siv
        
        #%     remove layers of type 1
        i = np.where(A == 0)
        ns=ns[i]
        num=num[i]
        roi=roi[i]
        tei=tei[i]
        di=di[i]
        Ti=Ti[i]
        pci=pci[i]
        Wi=Wi[i]
        gai=gai[i]
        si=si[i]
        sal=sal[i]
        FH=FH[i]
        FV=FV[i]
        sih=sih[i]
        siv=siv[i]
        fi=fi[i]
        epsi=epsi[i]
        epsii=epsii[i]
        dei=dei[i]
        N = len(di)
        tei = np.append(tei,theta)
        if (repo == 1):
          print num,di,roi,Ti,pci,ns,gai,si,sal,epsi[0:N],epsii,fi,tei[0:N]*180/np.pi,FH,FV,sih,siv
     
  else:
        #%disp(['I am inside slred']) 
        #%disp (['no coherent layers detected: ',num2str(freq),'GHz']);
        if (repo == 1):
          print [num,di,roi,Ti,pci,ns,gai,si,sal,epsi,epsii,fi,tei[0:N]*180/np.pi,sih,siv]

  return  num,roi,epsi,epsii,tei,sih,siv,di,dei,Ti,pci,Wi,gai,si,sal

def sccoeff(roi,Ti,pci,freq,Wi,gai,sccho):
#%   calculates the scattering coefficient from structural parameters
#%     different algorithms can be chosen, by changing "sccho"
#%
#%   [gbih,gbiv,gs6,ga2i] = sccoeff(roi,Ti,pci,freq,Wi,gai,sccho)
#%       gbih:  2-flux scattering coefficient at h pol
#%       gbiv:  2-flux scattering coefficient at v pol
#%       gs6:   6-flux scattering coefficient
#%       ga2i:  2-flux absorption coefficient
#%       roi:   density
#%       Ti:    physical temperature
#%       pci:   correlation length
#%       freq:  frequency
#%       Wi:    wetness
#%       gai:   absorption coefficient
#%       sccho: scattering coefficient algorithm chosen
#%
#%   Version history:
#%      1.0b    wi 15.7.95
#%      1.0     wi 23.9.97 bug fixed
#%      1.1     wi 26.9.97 latest fit on experimental data was added (option 7)
#%      1.2     wi 13.10.97 option 8 added, adapted scattering of a shell/sphere to note 9/ver2 
#%      1.3     wi  4.11.97 option 9, 10 and 11 added 
#%      1.4     wi 27.05.98 born approximation added (borna.m)
#%                 12.03.07 bug in iborn shperes fixed (rtt)
#%
#%   Uses:
#%       borna, ro2epsd, mixmod
#%
#%
#%
#%   Copyright (c) 1997 by the Institute of Applied Physics,
#%   University of Bern, Switzerland

  #% constants
  c = 2.99
  roair = 0.001293
  roice = 0.917
  #% specular component of scattering coefficient
  #% usually 0 can be important in new snow!
  dgb0h = 0
  dgb0v = 0
  #% aus der Theorie scattering coefficient
  k = freq*(2*np.pi/0.299793)
  eice = 3.18
  vfi = roi/roice
  
  #% choose the scattering algorithm that should be used
  wahl = sccho
  
  
  [epsi,epsii] = ro2epsd(roi,Ti,freq)
  [epsi,epsii] = mixmod(freq,Ti,Wi,epsi,epsii)
  
  
  #% 6-flux scattering coefficient
  if wahl == 1:
     gs6 = ((130 * ((freq/50)**2.7)) * pci**3) / (roi**1.3 + 0.001)
  
  
  #%fit vom 26.8.97 auf alle Daten v-pol, > 11 GHz
  if wahl == 2:
     gs6 = 0.0704 * (pci**2.32)*(freq**1.68)*roi**(-0.63)

  
  #% for spheres: Mätzler, J. Appl. Phys. 83(11) 6111-6117 eqs 27+32 iborn
  epseff = (2.-eice+3.*vfi*(eice-1)+ np.sqrt((2.-eice+3.*vfi*(eice-1))**2+8.*eice))/4.
  sphe = (3./32)*(0.001*pci)**3*k**4*vfi*(1-vfi)*abs((2.*epseff+1)*(eice-1)/(2.*epseff+eice))**2
  if wahl == 4:
     gs6 = sphe
  
  #% for shells(new and recrystalized snow): 
  #% Mätzler, J. Appl. Phys. 83(11) 6111-6117 eq 39 iborn
  epseff = 1.+(vfi*(eice-1)*(2.+1/eice))/(3.-vfi*(1-1./eice))
  shel = abs(2./3 + 1./(3.*eice**2))*(0.001*pci)*k**2*vfi*(1-vfi)*(eice-1)**2./(16.*epseff)
  if wahl == 5:
     gs6 = shel
  
  #% as linearcombination
  if wahl == 6:
     a = 0.1664
     b = 0.2545
     gs6 = a*sphe+b*shel
  
  #%fit vom 26.9.97
  if wahl == 7:
     gs6 = 73.21 * (pci**3)*((freq/50)**2.68)*roi**(-1)

  #%fit vom 13.10.97
  if wahl == 8:
     gs6 = 136 * (pci**2.85) * ((freq/50)**2.5) / (roi + 0.001)
  
  #%fit vom 4.11.97 (without density)
  if wahl == 9:
     gs6 = 564 * (pci**3.0)* ((freq/50)**2.5)
  
  #%fit vom 4.11.97 (without density, uses corr. length from exp. fit!)
  if wahl == 10:
     gs6 = (3.16 * pci + 295 * (pci**2.5))* ((freq/50)**2.5)
  
  #%fit vom 4.11.97 (with density, uses corr. length from exp. fit!)
  if wahl == 11:
     gs6 = (9.20 * pci - 1.23 * roi + 0.54)**2.5 * ((freq/50)**2.5)
  
  omega = np.sqrt((epsi - 1.)/epsi)
  
  
  #%Born Approximation
  if wahl == 12:
    print 'Born approximation, missing the functions!'
#    kp = bornsnk(roi,0)
#    [gb6,gc6,gf6,gs6] = borna(k,vfi,pci,epsi,eice,epseff,kp)
  else:
    gb6 = 0.5 * gs6 * (1.-omega)
    gc6 = 0.25 * gs6 * omega
  
  gbiv = np.zeros(len(gs6))
  gbih = np.zeros(len(gs6))
  gtr = np.zeros(len(gs6))
  ga2 = np.zeros(len(gai))
  #%dgb2h = zeros(len(gs6))
  #%dgb2h = dgb2h + 0.25;
  #%dgb2v = zeros(size(gs6));
  #%dgb2v = dgb2v + 0.1;
  
  #% -> 2 Flux
  gtr = (4. * gc6) / (gai + 2. * gc6) #gc is coefficient for coupling between horiz and vert fluxes
  ga2i = gai * (1. + gtr) #two-flux absorption coefficient
  
  gbih = (gb6 + dgb0h) + gtr * gc6
  gbiv = (gb6 + dgb0v) + gtr * gc6
    
  return gbih,gbiv,gs6,ga2i

def iborn_s2p(e1,e2,eeff,v,k,pcc):
#  % improved born approximation by C. M�tzler (1998). J. Appl. Phys. 83(11),6111-7
#  %scattering coefficient of a collection of spherical
#  %particles with correlation length pcc using improved born approximation
  ss=(3. *pcc**3 *k**4 /32.) *v *(1-v) *abs(((e2-e1)*(2. *eeff+e1)) /(2. *eeff+e2))**2
  return ss
  
def scice(si,gbih,gbiv,gs6,ga2i,Ti,sal,freq,pci):
#%   calculates the scattering coefficient from structural parameters
  k=(2*3.14159) /(0.3 /freq)
  eice=3.15+0.002*1j
  T=Ti-273.15
  
  volb=Vb(T,sal)
  [eb,ebi]=ebrine(T,freq)
  ebri=eb-ebi*1j
  #%emis=eice_rn2p(eice,eb+ebi*i,volb);
  emis=eice_s2p(eice,eb+ebi*1j,volb)
  ags6=iborn_s2p(eice,ebri,emis,volb,k,pci*0.001)
  gs6=gs6-gs6*si+ags6*si
  return gbih,gbiv,gs6,ga2i

def scice_my(si,gbih,gbiv,gs6,ga2i,Ti,dens,freq,pci,sal):
    #%calculates the scattering coefficient of MY ice from structural parameters

    k=(2.*3.14159)/(0.3/freq)
    eice=3.15+0.002j
    #T=Ti-273.15;
    epsi=eice.real
    epsii=eice.imag

    #permittivity of saline ice
    [sepsi, sepsii] = sie(si,sal,Ti,freq,epsi,epsii)

    eice=sepsi+sepsii*1j

    vola=(0.926-dens)/0.926
    ### added by Clara , commented on 22.09., uncommented on 17.10.
    vola[vola<0] = 0.
    ###
    emis=eice_s2p(eice,1.0+0.0j,vola)
    ags6=iborn_s2p(eice,1.0+0.0j,emis,vola,k,pci*0.001)
    gs6=gs6-gs6*si+ags6*si
    return gbih,gbiv,gs6,ga2i

def pfadc(teta,di,epsi,gs6):
#%   calculates the effective path length in a layer
#%
#%   [dei,tei,tscat] = pfadc(teta,di,epsi,gs6)
#%       dei:  effective path length [m]
#%       tei:  local incidence angle
#%       tscat: scattering 
#%       teta: incidence angle at snow air interface
#%       di:   thickness [m]
#%       epsi: dielectric permittivity
#%       gs6:  6-flux scattering coefficient
#%
#%   Version history:
#%      1.0    wi 15.10.97
#%
#% 
#%   Uses: -
#%
#%
#%
#%   Copyright (c) 1997 by the Institute of Applied Physics,
#%   University of Bern, Switzerland
  N = len(di)
  ns = np.sqrt(epsi)
  costetasn = np.sqrt(1-(np.sin(teta)/ns)**2)
  cosc = np.sqrt(1-(1/ns)**2)
  costetasc = 0.5 * (1 + cosc)
  dei = di/costetasn
  
  #tauscat = zeros(len(epsi)+array([1,0]))
  tauscat = np.zeros(len(epsi)+1)
  tscat = np.zeros(len(epsi))
  costeta = np.zeros(len(epsi))
  
  for m in range(N-1,-1,-1):
     tauscat[m] = tauscat[m+1] + dei[m] * gs6[m]/2
     tscat[m] = np.exp(-1 * tauscat[m])
     costeta[m] = tscat[m] * costetasn[m] + (1-tscat[m]) * costetasc[m]
  
  tei = np.arccos(costeta)
  tei*180/np.pi
#  costeta;
#  tauscat;
#  tscat;  
  return dei,tei,tscat

def polmix(tscat,sih,siv):
#%   calculates the polarization mixing of the interface reflectivities
#%       of each layer (taking into account the first order scattering)
#%
#%   [sih,siv] = polmix(tscat,sih,siv)
#%       sih:   interface reflectivity at h-pol
#%       siv:   interface reflectivity at v-pol
#%       tscat: tau scat
#%
#%   Version history:
#%      1.0    wi 14.10.97
#%      1.1    wi  4.11.97  bug fix (layer numbering problem)
#%
#%   Uses: -
#%
#%
#%
#%   Copyright (c) 1997 by the Institute of Applied Physics, 
#%   University of Bern, Switzerland
  tscat = np.append(tscat,1)
  
  smean = 0.5 * (sih + siv)
  deltas = 0.5 * tscat * (sih - siv)
  sih = smean + deltas
  siv = smean - deltas
  return sih,siv

def rt(gai,gbi,dei):
    #%   calculates the layer reflectivity and transmissivity
    #%
    #%
    #%   [ri,ti] = rt(gai,gbi,dei)
    #%       ri:   layer reflectivity
    #%       ti:   layer transmissivity
    #%       gai:  absorption coefficient
    #%       gbi:  scattering coefficient
    #%       dei:  path length
    #%
    #%   Version history:
    #%      1.0    wi 15.7.95
    #%   
    #%   Uses: -
    #%
    #%
    #%
    #%   Copyright (c) 1997 by the Institute of Applied Physics,
    #%   University of Bern, Switzerland
    gamma = np.sqrt(gai * (gai + 2 * gbi))
    t0i = np.exp(gamma * dei * (-1))
    r0i = np.zeros(len(t0i))
    i = np.where(gbi > 0.00001)
    r0i[i] = gbi[i] / (gai[i] + gbi[i] + gamma[i])
    
    t02 = t0i**2
    r02 = r0i**2
    ri = r0i * (1 - t02)/ (1 - t02 * r02)
    ti = t0i * (1 - r02) / (1 - t02 * r02)
    return ri, ti

def layer(ri,s_i,ti,Ti,Tgnd,Tsky):
#%   calculates the upwelling brightness temperatures D (see Note 6)
#%
#%   D = layer(ri,s_i,ti,Ti,Tgnd,Tsky)
#%       D:    upwelling brightness temperature
#%       ri:   layer reflectivity
#%       s_i:   interface reflectivity
#%       ti:   layer transmissivity
#%       Ti:   physical temperature [K]
#%       Tgnd: brightness temperature of the soil below the snowpack
#%       Tsky: brightness temperature of the sky
#%
#%   Version history:
#%      1.0    wi 15.7.95
#%      1.1    wi 26.9.97  handles also the special case of a single layer now
#%      1.2    wi 02.03.99 fixed error in 1 layer handling 
#%  
#%   Uses: -
#%
#%
#%
#%   Copyright (c) 1997 by the Institute of Applied Physics, 
#%   University of Bern, Switzerland

    N = len(ri)
    ei = 1 - ri - ti
    if N < 1: 
      print 'ERROR: No scattering layer'
      return

    if N == 1: 
      #%   k1 = (ri(1)*(1-s_i(1))*Tgnd+ti(1)*(1-si(2))*Tsky+ei(1)*Ti(1))/(1-ri(1)*s_i(1));
      #%   k2 = ti(1)*s_i(1)/(1-ri(1)*s_i(1));
      #%   k3 = 1-ri(1)*s_i(2)-ti(1)*s_i(1)*k2;
      #%   D = (ti(1)*s_i(1)*k1+ti(1)*(1-si(1))*Tgnd+ri(1)*(1-s_i(2))*Tsky+ei(1)*Ti(1))/k3;
      k1 = (1-ri[0]*s_i[0])*(1-ri[0]*s_i[1]) - ti[0]*s_i[0]*ti[0]*s_i[1]
      D = ti[0]*s_i[0] * ((1-s_i[0])*ri[0]*Tgnd + (1-s_i[1])*Tsky*ti[0] + ei[0]*Ti[0]) / (k1) + (1-ri[0]*s_i[0]) * ((1-s_i[0])*Tgnd*ti[0] + (1-s_i[1])*Tsky*ri[0] + ei[0]*Ti[0]) / (k1);    
    else:
      M1 = np.diag(ri*s_i[0:N])
      H = np.diag(ti[0:N-1]*(1-s_i[1:N]))
      M1[0:N-1,1:N] = M1[0:N-1,1:N] + H
  
      I = np.eye(len(M1))
      
      M2 = np.diag(ti*s_i[1:N+1])
      H = np.diag(ri[1:N]*(1-s_i[1:N]))
      M2[1:N,0:N-1] = M2[1:N,0:N-1] + H
    
      M3 = np.diag(ti*s_i[0:N])
      H = np.diag(ri[0:N-1]*(1-s_i[1:N]))
      M3[0:N-1,1:N] = M3[0:N-1,1:N] + H
      
      M4 = np.diag(ri*s_i[1:N+1])
      H = np.diag(ti[1:N]*(1-s_i[1:N]))
      M4[1:N,0:N-1] = M4[1:N,0:N-1] + H
      
      E = ei * Ti
      E[0] = E[0] + ri[0]* (1 - s_i[0]) * Tgnd
      E[N-1] = E[N-1] + ti[N-1]* (1 - s_i[N])* Tsky
      
      F = ei * Ti
      F[0] = F[0]+ ti[0]* (1 - s_i[0]) * Tgnd
      F[N-1] = F[N-1] + ri[N-1]* (1 - s_i[N]) * Tsky
      
      try:
        M5 = M3.dot(np.linalg.inv(I - M1).dot(M2)) + M4      
        D = np.linalg.inv(I - M5).dot(M3.dot(np.linalg.inv(I - M1).dot(E)) + F)    
      except np.linalg.LinAlgError as err:
        if 'Singular matrix' in err.message:
          D = np.zeros(Ti.shape)
          D[:] = np.nan
        else:
          raise

    return D

def transang(ia,e1,e2):
    #conputes the transmission angle using Snellius
    ta=np.asin((np.sqrt(e1.real)/np.sqrt(e2.real))*np.sin(ia))
    return ta

def extloss(d,ang,ke):
    #computes the loss coefficient in a layer with extinction, ke, and thickness, d.
    loss=(ke*d/np.cos(ang))
    return abs(loss)

def absloss(d,ang,ka):
    #computes the absorption loss coefficient in a layer with the absorption coefficient ka, and thickness, d.
    absorp=(ka*d/np.cos(ang))
    return abs(absorp)

def eice_rn2p(eh,ei,vi):
    # Polder / Van Santen dielectric constant of sea ice with random brine needles
    # e.g. Shokr (1998) Arctic sea ice in the microwave c-band, IEEE TGRS 36(2), 463-78.
    measureA=1.2
    measureB=1.2
    # initial estimate: Vant et al. (1978), J. Appl. Phys. 49(3), 1264-80.
    #emi=complex((3.05+0.72*vi),(0.024+0.33*vi))
    
    #inserted by Clara
    emi = (3.05+0.72*vi) + 1j * (0.024+0.33*vi)
    ####
    
    m=0
    # iterate to solve for emi: the dielectric const. of the mixture
    for n in range(len(emi)):
        measureA=1.2 #added by Clara
        measureB=1.2 #added by Clara
        while measureA > 0.001 and measureB > 0.0001:
            f1=(ei[n]-eh[n])/(ei[n]+emi[n])
            f2=5.0*emi[n]+ei[n]
            est=eh[n]+((vi[n]/3.0)*f1*f2)
            measureA=abs(emi[n].real-est.real)
            measureB=abs(emi[n].imag-est.imag)
            emi[n]=est
            m=m+1
            if (m>30): break
    return emi

def epicesal(T,freq,sal):
    # Dielectric constant of pure ice: Matzler, 1998, Microwave properties of ice and snow, in:
    # B. Smitht et al. (eds.) solar system ices, 241-257, Kluwer.
    # T: thermometric temperature in K
    # freq: Frequency in GHz
    epi=3.1884+9.1E-4*(T-273.15)

    # The Hufford model for the imagenary part.
    if T == 0: T=273.15
    theta=300.0/T
    alpha=(0.00504+0.0062*theta)*np.exp(-22.1*theta)
    beta=(0.502-0.131*theta/(1.+theta))*1e-4+(0.542e-6*((1.+theta)/(theta+0.0073))**2)

    epii=(alpha/freq) + (beta*freq)

    epui=complex(epi,epii)

    eh=complex(1.0,0.0)
    v=Vb((T-273.15),sal)
    eb=ebrine2(freq,(T-273.15),sal)
    e=eice_rn2p(epui,eb,v)
    if sal > 0.0: e=epui
    return epui

def e_wet_snow(k,Dens,wc):
    # beregner dielektrisitetskonstanten af vad sne
    fws=k/20.94
    denws=Dens/1000.0
    wc=wc*100.0
    A1ws=1.0
    A2ws=1.0
    B1ws=0.0
    Aws=1.0+1.83*denws+0.02*A1ws*(wc**1.015)+B1ws
    Bws=0.073*A1ws
    Cws=0.073*A2ws
    f0=9.07 #relaxation frequency

    ews1=Aws+((Bws*wc)/(1.0+((fws/f0)**2)))
    ews2=(Cws*(fws/f0)*wc)/(1.0+((fws/f0)**2))
    e_ws=complex(ews1,ews2)
    return e_ws

def FT_(ia,ta,e1,e2):
    #computes the Transmission coefficient (HH), using angles and dielectrics
    nn1=1.0/np.sqrt(e1.real)
    nn2=1.0/np.sqrt(e2.real)
    T_=(2.0*nn2*np.cos(ia))/(nn2*np.cos(ta)+nn1*np.cos(ta))
    return abs(T_)

def FTll(ia,ta,e1,e2):
    #computes the transmission coefficient (VV) using the angle and dielectrics
    nn1=1.0/np.sqrt((e1.real))
    nn2=1.0/np.sqrt((e2.real))
    Tll=(2.0*nn2*np.cos(ia))/(nn2*np.cos(ta)+nn1*np.cos(ia))
    return abs(Tll)

def FRll(ia,ta,e1,e2):
    #computes the Fresnell reflection coefficient for vertical polarisation
    nn1=1.0/np.sqrt(e1.real)
    nn2=1.0/np.sqrt(e2.real)
    r=abs((nn1*np.cos(ia)-nn2*np.cos(ta))/(nn1*np.cos(ia)+nn2*np.cos(ta)))
    return r

def FR_(ia,ta,e1,e2):
    #computes the Fresnell reflection coefficient for horisontal polarisation
    nn1=1.0/np.sqrt(e1.real)
    nn2=1.0/np.sqrt(e2.real)
    r=abs((nn2*np.cos(ia)-nn1*np.cos(ta))/(nn2*np.cos(ia)+nn1*np.cos(ta)))
    return r

def HWabsorp(e,k):
    kk=np.sqrt(e)
    #print e, k*abs(cmath.sqrt(kk)), abs(kk.imag)
    #a=2.0*k*abs(kk)*abs(kk.imag)
    a=2.0*k*abs(kk.imag)
    return a

def iborn_s2p_sc(e1,e2,eeff,v,k,pcc):
    # improved born approximation by C. Matzler (1998). J. Appl. Phys. 83(11),6111-7
    #scattering coefficient of a collection of spherical
    #particles with correlation length pc using improved born approximation
    ss=(3.0*(pcc**3.0)*(k**4)/32.0)*v*(1.0-v)*abs(((e2-e1)*(2.0*eeff+e1))/(2.0*eeff+e2))**2
    return ss

def iborn_s2p_bc(e1,e2,eeff,v,k,pcc):
    # improved born approximation by C. Matzler (1998). J. Appl. Phys. 83(11),6111-7
    #backscattering coefficient of a collection of spherical
    #particles with correlation length pcc using improved born approximation
    R=3.0*pcc/4.0
    ss=((R**3.0)*(k**4.0/3.0))*v*(1.0-v)*abs(((e2-e1)*(2.0*eeff+e1))/(2.0*eeff+e2))**2.0
    return ss

def meteo_sc(rroi,rTi,rpci,freq,rWi,rgai,rmeteo,gbih,gbiv,gs6,ga2i):
    #the scattering coefficient of only partly recrystalized snow, linear combination of iborn, wahl==6
    [dumgbih,dumgbiv,ags6,dumga2i] = sccoeff(rroi,rTi,rpci,freq,rWi,rgai,6)
    gs6=gs6-gs6*rmeteo+ags6*rmeteo
    return gbih,gbiv,gs6,ga2i 

def recry_sc(rroi,rTi,rpci,freq,rWi,rgai,rrecry,gbih,gbiv,gs6,ga2i):
    #the scattering coefficient of fully recrystalized snow, iborn, wahl==5
    [dumgbih,dumgbiv,ags6,dumga2i] = sccoeff(rroi,rTi,rpci,freq,rWi,rgai,5)
    gs6=gs6-gs6*rrecry+ags6*rrecry
    return gbih,gbiv,gs6,ga2i

def absorp2f(gbih,gbiv,gs6,ga2i,epsi,epsii,roi,Ti,pci,freq,Wi,gai):
    #%   calculates the scattering coefficient from structural parameters
    #%
    #%   [gbih,gbiv,gs6,ga2i] = absorp2f(gbih,gbiv,gs6,ga2i,roi,Ti,pci,freq,Wi,gai)
    #%       gbih:  2-flux scattering coefficient at h pol
    #%       gbiv:  2-flux scattering coefficient at v pol
    #%       gs6:   6-flux scattering coefficient
    #%       ga2i:  2-flux absorption coefficient
    #%       roi:   density
    #%       Ti:    physical temperature
    #%       pci:   correlation length
    #%       freq:  frequency
    #%       Wi:    wetness
    #%       gai:   absorption coefficient
    #%       sccho: scattering coefficient algorithm chosen
    #%
    if roi[0]>10:
      roi=roi/1000.0
    #% constants
    c = 2.99
    roair = 0.001293
    roice = 0.926
    #% specular component of scattering coefficient
    #% usually 0 can be important in new snow!
    dgb0h = 0
    dgb0v = 0
    #% aus der Theorie scattering coefficient
    k = freq*(2.*np.pi/0.299793)
    eice = 3.18
    vfi = roi/roice

    omega = np.sqrt((epsi - 1)/epsi)

    gb6 = 0.5* gs6 * (1-omega)
    gc6 = 0.25* gs6 * omega

    gbiv = np.zeros(np.size(gs6))
    gbih = np.zeros(np.size(gs6))
    gtr = np.zeros(np.size(gs6))
    ga2 = np.zeros(np.size(gai))

    #% -> 2 Flux
    gtr = (4. * gc6) / (gai + 2.* gc6)
    ga2i = gai * (1 + gtr)

    gbih = (gb6 + dgb0h) + gtr * gc6
    gbiv = (gb6 + dgb0v) + gtr * gc6

    return gbih,gbiv,gs6,ga2i

def penetration_depth(num,di,ke,theta,trans):
    #L=exp(ke*d*sec(theta))
    #sec = 1/cos
    # num : layers
    # di : layer thickness
    # ke : gs6+ga2i (?) - extinction coefficient
    # theta : rtei (?) - propagation angle
    # trans : ti (?)
    nn = len(num)
    inten=np.ones(len(num))
    pene=np.zeros(len(num)) 
    einv=0.367879 #1/e
    x=0.0
    loss=np.exp(ke*di*(1/np.cos(theta)))
    trans=np.array(trans) #probably ti
#    #from Rasmus
#    for idx in range(nn):
#      inten[idx]=np.prod(1.0/loss[0:idx])*np.prod(trans[0:idx])
#      pene[idx]=np.sum(di[0:idx])
    #modified by Clara
    for idx in range(len(num)-1,-1,-1):
      #print idx,'/',nn
      inten[idx]=np.prod(1.0/loss[idx:nn])*np.prod(trans[idx:nn]) #looks at the radiation
      pene[idx]=sum(di[idx:nn]) #depth of the interface where 1/e is reached
      #print inten
      
      # looks more precisely into the layer, where this is (by assuming a linear profile)
      if inten[idx] <= einv:
          #print inten[idx]
          #print idx,'/',nn
          #print pene[idx],'/',sum(di)
          if idx==nn-1:
              #print 'loop 1'
              slope=(1.0-inten[idx])/di[idx]
              x=(1.0-einv)/slope
              #print x
              break
          else:
              #print 'loop 2'
              slope=(inten[idx]-inten[idx+1])/di[idx]
              x=pene[idx+1]+(inten[idx+1]-einv)/slope
              #print x
              break
    return x

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
    
def icemain(f1,f2,teta,s0h,s0v,ifile,Tsky,Tgnd,sccho):

    #%   MEMLS
    #%   main program, plots the brightness temperature
    #%   of a snowpack/and packice by frequency.
    #%   
    #%   mainly designed for 5 to 100 GHz
    #%
    #%   icemain(f1,f2,teta,s0h,s0v,ifile,Tsky,Tgnd,sccho)
    #%     f1:    start frequency [GHz]
    #%     f2:    stop frequency [GHz]
    #%     teta:  incidence angle [deg]
    #%     s0h:   snow-ground reflectivity 
    #%     s0v:   snow-ground reflectivity
    #%     ifile: Inputfile with layer- number, temp [K], density [kg/m3],
    #%                                  thickness [cm], cor. length [mm], exp. corr. length [mm] 
    #%     Tsky:  sky brightness temperature [K]
    #%     Tgnd:  ground temperature [K]
    #%     sccho: type of scattering coefficient (11 recommended)
    #%
    #%   Version history:
    #%      1.0    wi 19.03.97
    #%      2.0    wi 25.09.97 adapted to the new subroutines
    #%      2.1    wi 15.10.97 adapted to the new subroutines
    #%      2.2    wi 04.11.97 adapted to bug in polmix
    #%   Uses:
    #%       ro2epsd, mixmod, abscoeff, pfadi, pfadc, polmix, fresnelc, slred, sccoeff, rt, layer, (signatur)
    #%
    #%
    #%
    #%   Copyright (c) 1997 by the Institute of Applied Physics, 
    #%   University of Bern, Switzerland
    
    #%  aux. input parameters:

  teta = (teta*np.pi)/180.
  
  #ifile ='/work/mh0033/m300411/SatSim/MEMLS/INPUT/70N00W-p2/timestep1234_orig1lay.dat'
  #ifile = '/work/mh0033/m300411/SatSim/MEMLS/INPUT/syntes.dat'
  data = pd.read_table(ifile, delimiter=" ",names= ['num','Ti','Wi','roi','di','pci','epci','sal','si'])
  #data.columns = ['num','Ti','Wi','roi','di','pci','epci','sal','si']
  num = data['num'].values.copy()
  Ti = data['Ti'].values.copy()
  Wi = data['Wi'].values.copy()
  roi = data['roi'].values.copy()
  di = data['di'].values.copy()
  pci = data['pci'].values.copy()
  epci = data['epci'].values.copy()
  sal = data['sal'].values.copy()
  si = data['si'].values.copy()

  
  
  #check for di>0 - NEEDED?
  
  N = len(data['num'])
  
  #print 'N= ',N
  #print 'roi= ',roi
  #print 'pci= ',pci
  #print 'sal= ',sal
  
  #if no layers, finish
  if N == 0: 
      print 'leaving the function' 
      return
  
  
  # convert d in m and roi in kg/L (?) 
  roi = roi/1000.;
  di = di/100.;
  
  #must be added if layers are given by height
  #di = zz - [0,zz(1:N-1)]
  #eval (['zz = ',ifile(1:6),'(i,5);'])
  
  #create indices over frequency range
  x  = range(f1,f2+1)
  yh = range(f1,f2+1)
  yv = range(f1,f2+1)
  
  #loop over frequency
  for freq in range(f1,f2+1):
    #print freq
    
    num = data['num'].values.copy()
    Ti = data['Ti'].values.copy()
    Wi = data['Wi'].values.copy()
    roi = data['roi'].values.copy()
    di = data['di'].values.copy()
    pci = data['pci'].values.copy()
    epci = data['epci'].values.copy()
    sal = data['sal'].values.copy()
    si = data['si'].values.copy()
    roi = roi/1000.;
    di = di/100.;
 
    ### CHECK
    #print roi
    #print sal
    
    #epsi and epsii if everything was dry snow  
    [epsi,epsii] = ro2epsd(roi,Ti,freq)
    #epsi and epsii if wet
    [epsi,epsii] = mixmod(freq,Ti,Wi,epsi,epsii)
    #epsi and epsii for the ice layers
    [epsi,epsii] = sie(si,sal,Ti,freq,epsi,epsii)
    
    plt.figure()
    plt.plot(epsi)
    
    #absorption coefficient
    gai = abscoeff(epsi,epsii,Ti,freq,Wi)
    #local incidence angle
    ns = np.sqrt(epsi)
    tei = np.append(np.array([np.arcsin(np.sin(teta)/ns)]),teta)
    #effective path angle
    dei = pfadi(tei,di)
    #add a 1 to epsi to have same dimensions for tei and epsi, get fresnel
    #coefficients (modified to fresnelc0 because conflict with other functions)
    [sih,siv] = fresnelc0(tei,np.append(epsi,1))
    #clean out coherent layers (?)
    [rnum,rroi,repsi,repsii,rtei,rsih,rsiv,rdi,rdei,rTi,rpci,rWi,rgai,rsi,rsal] = slred(num,roi,epsi,epsii,tei,sih,siv,di,dei,Ti,pci,freq,Wi,gai,si,sal)
    
    #added/adapted in Version 2.1
    #scattering coefficients for snow
    [gbih,gbiv,gs6,ga2i] = sccoeff(rroi,rTi,rpci,freq,rWi,rgai,sccho)
    #scattering coefficients for ice
    [gbih,gbiv,gs6,ga2i] = scice(rsi,gbih,gbiv,gs6,ga2i,rTi,rsal,freq,rpci)
    #effective path length (again?) with tau scat and new local incidence
    #angle
    [rdei,rtei,tscat] = pfadc(teta,rdi,repsi,gs6)
    #combine interface and surface reflectivities
    rsih = np.append([s0h],[rsih])
    rsiv = np.append([s0v],[rsiv])
    
    #effect of polarization mixing
    [rsih,rsiv] = polmix(tscat,rsih,rsiv)
    # ????
    rtei*180./np.pi
    #layer reflectivity and transmissivity (horizontal)
    [ri,ti]  = rt(ga2i,gbih,rdei)
    
    #Brightness temperature at every layer (horizontal)
    Dh   = layer(ri,rsih,ti,rTi,Tgnd,Tsky)
    # amount of layers
    N = len(rroi)
    #TBH at top
    if N > 1:
        Tbh = (1-rsih[N])*Dh[N-1] + rsih[N]*Tsky
    else:
         Tbh = (1-rsih[N])*Dh + rsih[N]*Tsky
    #store it in an arry for the different frequencies
    yh[freq-f1] = round(Tbh,2)
        
    #layer reflectivity and transmissivity (vertical)
    [ri,ti]  = rt(ga2i,gbiv,rdei)
    #tscat
    #rsih
    #rsiv
    
    #Brightness temperature at every layer (vertical)
    Dv   = layer(ri,rsiv,ti,rTi,Tgnd,Tsky)
    #TBV at top
    if N > 1:
        Tbv = (1-rsiv[N])*Dv[N-1] + rsiv[N]*Tsky
    else:
        Tbv = (1-rsiv[N])*Dv + rsiv[N]*Tsky
    #store it in an array for the different frequencies
    yv[freq-f1] = round(Tbv,2)
    #[freq,N,Tbv]
    #[Dv;Tbv]

#  print x
#  print yv
#  print yh
  
  return  x,yv,yh

"""  
ifile = 'syntes_modif.dat'
workdata = pd.read_table(ifile, delimiter=" ",names= ['num','Ti','Wi','roi','di','pci','epci','sal','sitype','si'])
#data.columns = ['num','Ti','Wi','roi','di','pci','epci','sal','si']
num = workdata['num'].values.copy()
Ti = workdata['Ti'].values.copy()
Wi = workdata['Wi'].values.copy()
roi = workdata['roi'].values.copy()
di = workdata['di'].values.copy()
pci = workdata['pci'].values.copy()
epci = workdata['epci'].values.copy()
sal = workdata['sal'].values.copy()
sitype = workdata['sitype'].values.copy()
si = workdata['si'].values.copy()
"""

def memls_mod(freqs,num,di,Ti,Wi,roi,pci,sal,sitype,si):
    #### CHANGE DEPENDING ON WHAT INPUT WE WANT
    teta=55
    s0h=0.75
    s0v=0.25
    Tsky=0
    Tgnd=271.35
    sccho=2 #4: iborn for spheres, 5: iborn for shells and spheres, now a function of type.
    #################
    
    NUM = num.copy()
    DI = di.copy()
    TI = Ti.copy()
    WI = Wi.copy()
    ROI = roi.copy()
    PCI = pci.copy()
    SAL = sal.copy()
    SITYPE = sitype.copy()
    SI = si.copy()
    
    teta = (teta * np.pi) / 180
    
    N = len(num)
    if N == 0: 
      print 'leaving the function' 
      return
    
#    roi = roi/1000.
#    di = di/100. #??? needed ??? - I think yes! missing in the matlab script, maybe not actually, we have it already in cm in our files

    #loop over frequency
    #freqs= [1.4, 6.9, 10.7, 18.7, 23.8, 36.5, 50.0, 89.0]#, 157.0, 183.0]
    #freqs= [6.9, 10.7, 18.7, 23.8, 36.5, 50.0, 89.0]#, 157.0, 183.0]

    x  = np.array(freqs,dtype='float')
    yh = np.array(freqs,dtype='float')
    yv = np.array(freqs,dtype='float')
    yh0 = np.array(freqs,dtype='float')
    yv0 = np.array(freqs,dtype='float')
    yh100 = np.array(freqs,dtype='float')
    yv100 = np.array(freqs,dtype='float')
    tb = np.array([freqs,freqs],dtype='float').flatten() #To hold the tbs in runalg sequence
    emis = np.array([freqs,freqs],dtype='float').flatten() #To hold the tbs in runalg sequence
    pend_h = np.array(freqs,'float')
    pend_v = np.array(freqs,'float')

    for ifreq in range(len(freqs)):
          freq = freqs[ifreq]
          #print freq
          num = NUM.copy()
          di = DI.copy()#/100.
          Ti = TI.copy()
          Wi = WI.copy()
          roi = ROI.copy()/1000.
          pci = PCI.copy()
          sal = SAL.copy()
          sitype = SITYPE.copy()
          si = SI.copy()

          #roi = roi/1000.
          #di = di/100.

          ### CHECKING EVERYTHING IS FINE
#          print ifreq
#          print 'roi', roi
#          print 'di', di
          ###
          #epsi and epsii if everything was dry snow  
          [epsi,epsii] = ro2epsd(roi,Ti,freq)    
          #epsi and epsii if wet
          [epsi,epsii] = mixmod(freq,Ti,Wi,epsi,epsii)
          #select dielectric scheme for FY or MY ice
          fy = np.zeros(len(num))
          fy[sitype==3] = 1
          my = np.zeros(len(num))
          my[sitype==4] = 1
          #epsi and epsii for the ice layers
          [epsi,epsii] = sie(fy,sal,Ti,freq,epsi,epsii)
          [epsi,epsii] = mysie(my,roi,Ti,sal,freq,epsi,epsii)
          
#          plt.figure()
#          plt.plot(epsi)
#          plt.title('clara new')
#          plt.ylim(0, 0.4)
#          plt.grid('on')
          
          #absorption coefficient
          gai = abscoeff(epsi,epsii,Ti,freq,Wi) #ga is absorption coefficient
          
          #local incidence angle
          ns = np.sqrt(epsi) #real part of the refractive index of the slab
          tei = np.append(np.array([np.arcsin(np.sin(teta)/ns)]),teta) #np.arcsin(...) is the critical angle for total reflection
          #effective path angle
          dei = pfadi(tei,di)
          #add a 1 to epsi to have same dimensions for tei and epsi, get fresnel
          #coefficients (modified to fresnelc0 because conflict with other functions)
          [sih,siv] = fresnelc0(tei,np.append(epsi,1))
          #clean out coherent layers (?)
          [rnum,rroi,repsi,repsii,rtei,rsih,rsiv,rdi,rdei,rTi,rpci,rWi,rgai,rsitype,rsal] = slred(num,roi,epsi,epsii,tei,sih,siv,di,dei,Ti,pci,freq,Wi,gai,sitype,sal)
          
          #added/adapted in Version 2.1
          #scattering stuff: gb is scattering coefficient, ga2i absorption coefficient (?)
          #scattering coefficients for snow
          [gbih,gbiv,gs6,ga2i] = sccoeff(rroi,rTi,rpci,freq,rWi,rgai,sccho)
          rmeteo = np.zeros(len(rnum))
          rmeteo[rsitype==1] = 1
          rrecry = np.zeros(len(rnum))
          rrecry[rsitype==2] = 1
          #scattering coefficients for snow
          [gbih,gbiv,gs6,ga2i] = meteo_sc(rroi,rTi,rpci,freq,rWi,rgai,rmeteo,gbih,gbiv,gs6,ga2i)
          [gbih,gbiv,gs6,ga2i] = recry_sc(rroi,rTi,rpci,freq,rWi,rgai,rrecry,gbih,gbiv,gs6,ga2i)

         # plt.figure()
         # plt.plot(ga2i,'-k')
#         plt.plot(gbiv,'--b')
         # plt.title('claras version 1') - som om alt var sne
         # plt.grid('on')

          #select MY or FY ice scattering modules
          rfy = np.zeros(len(rnum))
          rfy[rsitype==3] = 1
          rmy = np.zeros(len(rnum))
          #rmy[rsitype==4] = 1
          rmy[rsitype>3] = 1 #weird bug, changed on 19.12.2017
          #scattering coefficients for ice
          [gbih,gbiv,gs6,ga2i] = scice(rfy,gbih,gbiv,gs6,ga2i,rTi,rsal,freq,rpci)
          [gbih,gbiv,gs6,ga2i] = scice_my(rmy,gbih,gbiv,gs6,ga2i,rTi,rroi,freq,rpci,rsal)
          #recompute the 2flux absorption coefficient 
          [gbih,gbiv,gs6,ga2i] = absorp2f(gbih,gbiv,gs6,ga2i,repsi,repsii,rroi,rTi,rpci,freq,rWi,rgai)
          
     
#             plt.figure()
#          plt.plot(rsih)
#          plt.title('claras version')
#          plt.grid('on')
             
          #effective path length (again?) with tau scat and new local incidence
          #angle
          [rdei,rtei,tscat] = pfadc(teta,rdi,repsi,gs6)
          ### added by Clara
          #rtei[np.isnan(rtei)==True] = 0.
          ###
          #combine interface and surface reflectivities
          rsih = np.append([s0h],[rsih])
          rsiv = np.append([s0v],[rsiv])
          
          #effect of polarization mixing
          [rsih,rsiv] = polmix(tscat,rsih,rsiv)
         
#          ### added by Clara
#          rsih[rsih>1.] = 1.
#          rsiv[rsiv>1.] = 1.
#          rsih[rsih<0] = 0.
#          rsiv[rsiv<0] = 0.
#          ###
          # ????
          rtei*180./np.pi
          
#          plt.figure()
#          plt.plot(rsih)
#          plt.title('claras version')
#          plt.grid('on')
          
          #layer reflectivity and transmissivity (horizontal)
          [ri,ti]  = rt(ga2i,gbih,rdei)
          ### added by Clara , commented on the 22.09., uncommented 17.10.
          if len(ri) > 1:
            for i in range(len(ri)-2,-1,-1):
              if np.isnan(ri[i])==True:
                ri[i] = ri[i+1]
              if np.isnan(ti[i])==True:
                ti[i] = ti[i+1]
          ###    
          rih = ri
          tih= ti
          #Brightness temperature at every layer (horizontal)
          Dh   = layer(ri,rsih,ti,rTi,Tgnd,Tsky)
         
          # amount of layers
          N = len(rroi)
          #TBH at top
          if N > 1:
              Tbh = (1-rsih[N])*Dh[N-1] + rsih[N]*Tsky
          else:
               Tbh = (1-rsih[N])*Dh + rsih[N]*Tsky    
          #store it in an arry for the different frequencies
          yh[ifreq] = Tbh #round(Tbh,2)
          print yh[ifreq]
          tb[2*ifreq+1] = Tbh #round(Tbh,2)
          
          #print rtei
          if (len(rnum) > 1) or (len(rnum)==1 and np.isnan(rdi[0])==False):
            pend_h[ifreq] = penetration_depth(rnum,rdi,ga2i+gs6,rtei,ti)
          else:
            pend_h[ifreq] = np.nan
          #print freq, pend[ifreq]
          
          #layer reflectivity and transmissivity (vertical)
          [ri,ti]  = rt(ga2i,gbiv,rdei)
          ### added by Clara , commented on the 22.09., uncommented 17.10.
          if len(ri) > 1 :
            for i in range(len(ri)-2,-1,-1):
              if np.isnan(ri[i])==True:
                ri[i] = ri[i+1]
              if np.isnan(ti[i])==True:
                ti[i] = ti[i+1]
          riv = ri
          tiv = ti
                
#          plt.figure()
#          plt.plot(ti)
#          plt.title('claras version')
#          plt.grid('on')
          
          ###
          #tscat
          #rsih
          #rsiv        
          #Brightness temperature at every layer (vertical)
          Dv   = layer(ri,rsiv,ti,rTi,Tgnd,Tsky)
          
          
          # amount of layers
          N = len(rroi)
          #TBV at top
          if N > 1:
              Tbv = (1-rsiv[N])*Dv[N-1] + rsiv[N]*Tsky
          else:
              Tbv = (1-rsiv[N])*Dv + rsiv[N]*Tsky
          #store it in an array for the different frequencies
          yv[ifreq] = Tbv #round(Tbv,2)
          tb[2*ifreq] = Tbv #round(Tbv,2)

          #print rtei
          if (len(rnum) > 1) or (len(rnum)==1 and np.isnan(rdi[0])==False):
            pend_v[ifreq] = penetration_depth(rnum,rdi,ga2i+gs6,rtei,ti)
          else:
            pend_v[ifreq] = np.nan
          #print freq, pend[ifreq]
              
          #Emissivities
          eDh   = layer(ri,rsih,ti,rTi,Tgnd,0)
          if N>1:
            eTbh = (1-rsih[N])*eDh[N-1]
          else:
            eTbh = (1-rsih[N])*eDh
          yh0[ifreq] = eTbh
          
          eDh = layer(ri,rsih,ti,rTi,Tgnd,100)
          if N>1:
            eTbh = (1-rsih[N])*eDh[N-1] + rsih[N]*100.
          else:
            eTbh = (1-rsih[N])*eDh + rsih[N]*100.
          yh100[ifreq] = eTbh
          
          eDv   = layer(ri,rsiv,ti,rTi,Tgnd,0)
          if N>1:
            eTbv = (1-rsiv[N])*eDv[N-1]
          else:
            eTbv = (1-rsiv[N])*eDv
          yv0[ifreq] = eTbv
          
          eDv   = layer(ri,rsiv,ti,rTi,Tgnd,100)
          if N>1:
            eTbv = (1-rsiv[N])*eDv[N-1] + rsiv[N]*100.
          else:
            eTbv = (1-rsiv[N])*eDv + rsiv[N]*100.
          yv100[ifreq] = eTbv
          

    
    #Calculate emissivities, formulas can be found in Wiesmann and Mätzler, 1999
    rv = (yv100 - yv0)/100.
    rh = (yh100 - yh0)/100.
    ev = 1 - rv
    eh = 1 - rh
    Teffv = yv0/ev
    Teffh = yh0/eh
    
    return x, eh, ev, yh, yv, Teffh, Teffv, pend_h, pend_v,epsi, epsii, gs6, gbih, gbiv, ga2i,rih,riv,tih,tiv,rsih,rsiv,rtei,rdei

def Vb_for_ndims(T,S):
    dims = T.shape
    T_flat = T.values.flatten()
    S_flat = S.values.flatten()
    lwf_flat = Vb(T_flat,S_flat)
    lwf = lwf_flat.reshape(dims)
    lwf2 = xr.DataArray(lwf,dims=T.dims)
    lwf2['time'] = T.time
    lwf2['depth'] = T.depth
    return lwf2

