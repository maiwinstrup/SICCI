#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 12:53:52 2018

@author: mwi
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
import numdifftools as nd
import multiprocessing as mp
import cartopy.crs as ccrs
import time
import optimale_old
import optimalest 
import sett_input
from icemodel import liseequation
from datetime import date

def validate_res(filename):
    
    #%% Import Round Robin datafile
    df=pd.read_csv(filename,sep=',',skiprows=[0], na_values=['noval',' noval','  noval','   noval','    noval','     noval','      noval'])
    # List all variables
    #df.keys()
    #list the variables where to drop nans
    nndf=df[[u'#latitude',u'longitude',u'6.9GHzH', u'6.9GHzV', u'7.3GHzH', u'7.3GHzV', u'10.7GHzH' , u'10.7GHzV',\
        u'18.7GHzH', u'18.7GHzV', u'23.8GHzH', u'23.8GHzV',u'36.5GHzH', u'36.5GHzV', u'89.0GHzH', u'89.0GHzV',\
        u'sd_mean',u'sd_std',u't2m']].dropna()
    
    #%% Extract data:
    lat = nndf['#latitude'].as_matrix()
    lon = nndf['longitude'].as_matrix()
    
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
    
    # Comparison parameters:
    sd_meas = nndf['sd_mean'].as_matrix() # Snow depth, measured
    sd_std_meas = nndf['sd_std'].as_matrix() # STD of snow depth, measured
    T2m = nndf['t2m'].as_matrix() # 2m air temperature, from ECMWF

    #%% Initialize arrays
    P_array = np.ones((len(Tb7h0),6))
    S_array = np.ones((len(Tb7h0),6))
    cost_array = np.ones((len(Tb7h0),1))
    SDsim = np.ones((len(Tb7h0),1))
    
    #%% Run loops:
    for x in xrange(0,len(Tb7h0)):  
        Tb7h = np.float(Tb7h0[x]) #6h
        Tb7v = np.float(Tb7v0[x]) #6v
        Tb11v = np.float(Tb11v0[x]) #10v
        Tb11h = np.float(Tb11h0[x]) #10h
        Tb19v = np.float(Tb19v0[x]) #19v
        Tb19h = np.float(Tb19h0[x]) #19v
        Tb24v = np.float(Tb24v0[x]) #23v
        Tb24h = np.float(Tb24h0[x]) #23v
        Tb37v = np.float(Tb37v0[x]) #37v
        Tb37h = np.float(Tb37h0[x]) #37v
        Tb89v = np.float(Tb89v0[x]) #37v!
        Tb89h = np.float(Tb89h0[x]) #37n!
         
        Tb = np.array([Tb7v,Tb7h,Tb11v,Tb11h,Tb19v,Tb19h,Tb24v,Tb24h,Tb37v,Tb37h,Tb89v,Tb89h])
        freq = np.array([7,11,19,24,37,89])
         
        if np.isnan(Tb[0]): continue
            
        #%% Optimize:
        # Set input values for parameters:
        inputpar=sett_input.setinputpar(Tb)
        P_final, S_diag, cost = optimalest.optimize(Tb,inputpar)
        
        # Original version of OE scheme:
        #P_final_old, S_diag_old, cost_old = optimale_old.optimal(Tb)
         
        # Add to array:
        P_array[x,:]=P_final
        S_array[x,:]=np.sqrt(S_diag) # Standard deviation (from variance)
        cost_array[x,0]=cost
        
        # First guess of snow depth:
        SDsim[x,0] = liseequation(Tb)
    
    # Save results in a dataframe for easy comparison:
    results = pd.DataFrame()
    results['lat'] = lat
    results['lon'] = lon
    results['sd_meas'] = sd_meas
    results['sd_std_meas'] =sd_std_meas
    results['sd'] = P_array[:,5]
    results['sd_std'] = S_array[:,5]
    results['sd_sim'] = SDsim
    results['T2m']=T2m
    results['Ts'] = P_array[:,0]
    results['Ts_std'] = S_array[:,0]
    results['V']=P_array[:,1]
    results['V_std']=S_array[:,1]
    results['W']=P_array[:,2]
    results['W_std']=S_array[:,2]
    results['L']=P_array[:,3]
    results['L_std']=S_array[:,3]
    results['c_ice']=P_array[:,4]
    results['c_ice_std']=S_array[:,4]
    return results
   
#%% Run for multiple days: 
if __name__ == "__main__":
    year = 2015
    month_day = (3,[19,24,25,26,27,29,30],4,[1,3])

    pool = mp.Pool(processes=4) # Allow 4 processes running at the same time 
    start_time = time.time()
    
    results_tot = {}
    i=0
    for month in month_day[::2]:
        i=i+1
        results = [pool.apply_async(
                validate_res,args=('./RRDP_v2.0/NERSC_OIB/QSCAT-vs-SMAP-vs-SMOS-vs-ASCAT-vs-AMSR2-vs-ERA-vs-NERSC-OIB-'+ 
                                   str(year)+ str(month).zfill(2) + str(day).zfill(2) +'.text',)) for day in month_day[2*i-1]]
        output = [p.get() for p in results] # Results are not sorted by date
        results_tot[i-1] = pd.concat(output)
    # Combine results in a dataframe:
    results_tot = pd.concat(results_tot)

    end_time = time.time()
    dt = end_time-start_time
    print 'dt', dt

    # Parallelization: try to avoid reading and writing in each - this slows it down