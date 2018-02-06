#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 12:53:52 2018

@author: mwi
"""

import numpy as np
import pandas as pd
import optimalest 
import optimale_old
import matplotlib.pyplot as plt
from matplotlib import cm
import pdb
import numdifftools as nd
import sett_input

def validate_res(filename):
    #tic = time.clock()
    
    #%% Import Round Robin datafile
    df=pd.read_csv(filename,sep=',',skiprows=[0], na_values=['noval',' noval','  noval','   noval','    noval','     noval','      noval'])
    # List all variables
    #df.keys()
    #list the variables where to drop nans
    nndf=df[[u'6.9GHzH', u'6.9GHzV', u'7.3GHzH', u'7.3GHzV', u'10.7GHzH' , u'10.7GHzV',\
        u'18.7GHzH', u'18.7GHzV', u'23.8GHzH', u'23.8GHzV',u'36.5GHzH', u'36.5GHzV', u'89.0GHzH', u'89.0GHzV',u'sd_mean',u'sd_std',u't2m']].dropna()
    
    #%% Extract data:
    Tb7h0 = nndf['6.9GHzH'].as_matrix()
    Tb7v0 = nndf['6.9GHzV'].as_matrix()
    Tb11h0 = nndf['10.7GHzH'].as_matrix()    
    Tb11v0 = nndf['10.7GHzV'].as_matrix()
    Tb19h0 = nndf['18.7GHzH'].as_matrix()    
    Tb19v0 = nndf['18.7GHzV'].as_matrix()
    Tb24h0 = nndf['23.8GHzH'].as_matrix()    
    Tb24v0 = nndf['23.8GHzV'].as_matrix()
    Tb37h0 = nndf['36.5GHzH'].as_matrix()    
    Tb37v0 = nndf['36.5GHzV'].as_matrix()
    Tb89h0 = nndf['89.0GHzH'].as_matrix()    
    Tb89v0 = nndf['89.0GHzV'].as_matrix()
    
    sd0 = nndf['sd_mean'].as_matrix() # Snow depth
    sd0_std = nndf['sd_std'].as_matrix() # Snow depth std
    
    T2m0 = nndf['t2m'].as_matrix() # Temperate at 2m, from ECMWF

    #%% Initialize arrays
    P_array = np.ones((len(Tb7h0),6))
    S_array = np.ones((len(Tb7h0),6))
    cost_array = np.ones((len(Tb7h0),1))
    sd_array = np.ones((len(Tb7h0),1))
    sd_std_array = np.ones((len(Tb7h0),1))
    SDsim = np.ones((len(Tb7h0),1))
    T_array = np.ones((len(Tb7h0),1))

    #P_array_old = np.ones((len(Tb7h0),6))        
    
    #%% Run loops:
    for x in xrange(0,len(Tb7h0)):  
         Tb6h = np.float(Tb7h0[x]) #6h
         Tb6v = np.float(Tb7v0[x]) #6v
         Tb10v = np.float(Tb11v0[x]) #10v
         Tb10h = np.float(Tb11h0[x]) #10h
         Tb19v = np.float(Tb19v0[x]) #19v
         Tb19h = np.float(Tb19h0[x]) #19v
         Tb23v = np.float(Tb24v0[x]) #23v
         Tb23h = np.float(Tb24h0[x]) #23v
         Tb37v = np.float(Tb37v0[x]) #37v
         Tb37h = np.float(Tb37h0[x]) #37v
         Tb89v = np.float(Tb89v0[x]) #37v!
         Tb89h = np.float(Tb89h0[x]) #37n!
         
         Tb = np.array([Tb6v,Tb6h,Tb10v,Tb10h,Tb19v,Tb19h,Tb23v,Tb23h,Tb37v,Tb37h,Tb89v,Tb89h])
         
         #%%
         if np.isnan(Tb[0]):
            continue
            
         #%% Optimize:
         #pdb.set_trace()
         # Set input values for parameters:
         inputpar=sett_input.setinputpar(Tb)
         P_final, S_diag, cost = optimalest.optimize(Tb,inputpar)

         #P_final_mid, S_diag_mid, cost_mid = optimale.optimal(Tb)
         P_final_old, S_diag_old, cost_old = optimale_old.optimal(Tb)
         
         # Add to array:
         P_array[x,:]=P_final
         S_array[x,:]=S_diag
         cost_array[x,0]=cost
         
         sd_array[x,0] = np.float(sd0[x]) 
         sd_std_array[x,0] = np.float(sd0_std[x])
         
         # First guess        
         SDsim[x,0] = 1.7701 + 0.017462*Tb6v - 0.02801*Tb19v + 0.0040926*Tb37v

         # Surface temperature (almost)
         T_array[x,0] = np.float(T2m0[x])
    
    #%% Append to file
    #sd_array_tot.append(sd_array)
    #sdcalc_array_tot.append(P_array[:,5])

    #%% Time runtime
    #toc=time.clock()
    #print toc-tic
    
    #%% Plot differences: #
#    plt.figure()
#    plt.plot(P_array[:,5],sd_array[:,0],'.k')
#   plt.xlabel('Calculated snow depth')
#    plt.ylabel('Observed snow depth')
#    plt.show()
#    plt.plot([0.05, 0.35], [0.05, 0.35],'k-')
    
#    plt.figure()
#    plt.scatter(P_array[:,5],sd_array[:,0],c=sd_std_array.ravel(),cmap=cm.jet)
#    plt.scatter(SDsim,sd_array[:,0],c='k')
#    plt.xlabel('Calculated snow depth')
#    plt.ylabel('Observed snow depth')
#    plt.title('Colored according to std(sd_obs)')
#    plt.plot([0.05, 0.35], [0.05, 0.35],'k-')
#    
#    plt.figure()
#    plt.scatter(P_array[:,5],sd_array[:,0],c=S_array[:,5].ravel(),cmap=cm.jet)
#    plt.scatter(SDsim,sd_array[:,0],c='k')
#    plt.xlabel('Calculated snow depth')
#    plt.ylabel('Observed snow depth')
#    plt.title('Colored according to S(sd_calc)')
#    plt.plot([0.05, 0.35], [0.05, 0.35],'k-')       
#    
#    plt.figure()
#    plt.scatter(P_array[:,0],T_array)
#    plt.xlabel('Calculated Ts')
#    plt.ylabel('ECMWF Ts')
    
    return sd_array, sd_std_array, P_array[:,5], S_array[:,5], SDsim, T_array, P_array[:,0], P_array 

#%% Run it for multiple days:
#a= [12,13,14,15,17,19,21]

a= [12,13,14,15,17,18,19,21,24,25,26,28,31] # also 25
sd_meas_tot = []
sd_calc_tot = []
sd_std_meas_tot = [] 
sd_S_calc_tot = []
sd_init_tot = []
P1_tot = []

#%%
for day in a:
    #filename = './RRDP_v2.0/NERSC_OIB/QSCAT-vs-SMAP-vs-SMOS-vs-ASCAT-vs-AMSR2-vs-ERA-vs-NERSC-OIB-201403' + str(day) +'.text'
    filename = './RRDP_v2.0/NERSC_OIB/QSCAT-vs-SMAP-vs-SMOS-vs-ASCAT-vs-AMSR2-vs-ERA-vs-NERSC-OIB-201403' + str(day) +'.text'
    sd_meas, sd_std_array, sd_calc, sd_S_calc, SDsim, T_meas, T_calc, P_array = validate_res(filename)

    sd_meas = np.squeeze(np.asarray(sd_meas))
    sd_std_array = np.squeeze(np.asarray(sd_std_array))
    
    sd_meas_tot.extend(sd_meas)
    sd_std_meas_tot.extend(sd_std_array)
    sd_calc_tot.extend(sd_calc)
    sd_S_calc_tot.extend(sd_S_calc)
    sd_init_tot.extend(SDsim) 
    
    P1_tot.extend(P_array[:,2])
    
#%%
plt.figure()
plt.scatter(sd_calc_tot,sd_meas_tot,c=sd_std_meas_tot,cmap=cm.jet)
plt.scatter(sd_init_tot,sd_meas_tot,c='k')
plt.xlabel('Calculated snow depth')
plt.ylabel('Observed snow depth')
plt.title('Colored according to std(sd_obs)')
plt.plot([0.05, 0.35], [0.05, 0.35],'k-')
 
#%%   
plt.figure()
plt.scatter(sd_calc_tot,sd_meas_tot,c=sd_S_calc_tot,cmap=cm.jet)
#plt.scatter(sd_init_tot,sd_meas_tot,c='k')
plt.xlabel('Calculated snow depth')
plt.ylabel('Observed snow depth')
plt.title('Colored according to S(sd_calc)')
plt.plot([0.05, 0.35], [0.05, 0.35],'k-')

#%% Histogram:
plt.figure()
bins = np.linspace(0,0.4,40)
plt.hist(sd_calc_tot,bins,facecolor='green',alpha=0.5,label='calculated')
plt.xlabel('Snow depth')
plt.ylabel('Occurrence')
plt.grid(True)
plt.hist(sd_meas_tot,bins,facecolor='blue',alpha=0.5,label='measured')
plt.legend(loc='upper right')

# Or next to each other:
data = np.vstack(np.asmatrix([sd_calc_tot, sd_meas_tot]).T)
plt.figure()
plt.hist(data,bins,alpha=0.7)