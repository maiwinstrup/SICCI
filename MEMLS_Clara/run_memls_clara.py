#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 15:19:41 2017

run MEMLS with the functions being in memls_functions (clalib)

@author: m300411
"""

import numpy as np
import pandas as pd
import glob
import xarray as xr
import clalib.memls_functions as memls


##########################################################################

exp_nb='018'

#went once through all experiments
#ee='75N00W-p2'
ee='NorthPole-p2'
ee2=ee.split("-")[0]


#on_mistral='yes'
on_mistral='no'

home_path='/home/mpim/m300411'
#home_path='/Users/claraburgard'

if on_mistral=='no':
    subfolder=glob.glob(home_path+'/mistral_home/SatSim/git_scripts/Documentation/Experiments/README_exp'+exp_nb+'*.txt')
    datet = subfolder[0].split("_")[5].split(".")[0]
    inputpath = home_path+'/mistral_work/SatSim/MEMLS_exp/exp'+exp_nb+'_'+ee2+'_'+datet+'/INPUT/netcdf_files/'
    outputpath = home_path+'/mistral_work/SatSim/MEMLS_exp/exp'+exp_nb+'_'+ee2+'_'+datet+'/OUTPUT/original_files/'
elif on_mistral=='yes':
    subfolder=glob.glob('/home/mpim/m300411/SatSim/git_scripts/Documentation/Experiments/README_exp'+exp_nb+'*.txt')
    datet = subfolder[0].split("_")[4].split(".")[0]
    inputpath = '/work/mh0033/m300411/SatSim/MEMLS_exp/exp'+exp_nb+'_'+ee2+'_'+datet+'/INPUT/netcdf_files/'
    outputpath = '/work/mh0033/m300411/SatSim//MEMLS_exp/exp'+exp_nb+'_'+ee2+'_'+datet+'/OUTPUT/original_files'


####################################################




sens_exp = ['complex','simpleallfunc','simpletemp','simplesalfunc','simpledens','simplesalconst','simpleallconst']

for sens in sens_exp:
  print sens+' started'
  netcdf_input = xr.open_dataset(inputpath+'inputMEMLS_'+ee2+'_'+sens_exp+'.nc')
  for tt,timet in enumerate(netcdf_input['time']):
    print str(tt)+' '+sens+' '+ee
    
    workdata=netcdf_input.sel(sens_exp=sens).isel(time=tt)
    
    num = range(1,len(workdata['temperature'].dropna('top_100_to_bottom_0').values)+1)

    if not num:
      frequency = [1.4, 6.9, 10.7, 18.7, 23.8, 36.5, 50.0, 89.0]
      eh = np.empty(len(frequency))
      eh[:] = np.nan
      ev = np.empty(len(frequency))
      ev[:] = np.nan
      TBH = np.empty(len(frequency))
      TBH[:] = np.nan
      TBV = np.empty(len(frequency))
      TBV[:] = np.nan
      TeffH = np.empty(len(frequency))
      TeffH[:] = np.nan
      TeffV = np.empty(len(frequency))
      TeffV[:] = np.nan
      pend_h = np.empty(len(frequency))
      pend_h[:] = np.nan
      pend_v = np.empty(len(frequency))
      pend_v[:] = np.nan
    else:
      #Ti = workdata['temperature'].dropna('top_100_to_bottom_0').values.copy()
      Ti = np.flipud(workdata['temperature'][0:num[-1]].values.copy())
      Wi = np.flipud(workdata['wetness'][0:num[-1]].values.copy())
      roi = np.flipud(workdata['density'][0:num[-1]].values.copy())
      pci = np.flipud(workdata['correlation_length'][0:num[-1]].values.copy())
      #pci = zeros(len(workdata['pci']))
      #pci[:] = 0.5
      epci = np.flipud(workdata['exp_correlation_length'][0:num[-1]].values.copy())
      di = np.flipud(workdata['thickness'][0:num[-1]].values.copy())
      sal = np.flipud(workdata['salinity'][0:num[-1]].values.copy())
      sitype = np.flipud(workdata['type'][0:num[-1]].values.copy())
      #sitype = zeros(len(workdata['sitype']))
      #sitype[:] = 3.
      si = np.flipud(workdata['snow_ice'][0:num[-1]].values.copy())
      
      frequency, eh, ev, TBH, TBV, TeffH, TeffV, pend_h, pend_v = memls.memls_mod(np.array(num),di,Ti,Wi,roi,pci,sal,sitype,si)
      #print 'pend_h =',pend_h
      #print 'pend_v =',pend_v
    
    fileinput = pd.DataFrame()
    fileinput['Freq'] = frequency
    fileinput['eh'] = eh
    fileinput['ev'] = ev
    fileinput['TBH'] = TBH
    fileinput['TBV'] = TBV
    fileinput['TeffH'] = TeffH
    fileinput['TeffV'] = TeffV
    fileinput['Pen. depth H'] = pend_h
    fileinput['Pen. depth V'] = pend_v
    #fileinput['Pen depth'] = pend
    #ff = outputpath+exp+'/timestep'+str(tt).zfill(4)+'_orig2lay.dat'
    ff = outputpath+'/output_timestep_'+str(tt).zfill(4)+'_'+sens+'.dat'
    np.savetxt(ff, fileinput.values, fmt='%1.2f') 
  print sens+' finished'
print ee