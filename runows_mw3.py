#!/usr/bin/python
"""
Computes maps of optimal values for sea-ice parameters from maps of measured 
brightness temperatures using optimal estimation. Reads and writes netCDF files. 

Usage:
 >>>   runows(<input_netcdf_path>, <output_netcdf_path>)
"""

import numpy as np
from netCDF4 import Dataset
from pyproj import Proj, transform
import optimalest
import sett_input

#%% "Runows"
def runows(nc_f_in, nc_f_out):

    #%% Import netCDF file:
    nc_fid = Dataset(nc_f_in,'r') # Dataset is the class behavior to open the file
                                  # and create an instance of the netCDF4 class        
    # Overview of file content:
    nc_attrs = nc_fid.ncattrs() # Global attributes
    nc_dims = [dim for dim in nc_fid.dimensions] # List of dimensions
    nc_vars = [var for var in nc_fid.variables] # List of variables
    
    #%% Extract data:
    xc = nc_fid.variables['xc'][:] 
    yc = nc_fid.variables['yc'][:]
    time = nc_fid.variables['time'][:]
    Tb6h0 = nc_fid.variables['tb06h'][:] 
    Tb6v0 = nc_fid.variables['tb06v'][:]
    Tb10h0 = nc_fid.variables['tb10h'][:]
    Tb10v0 = nc_fid.variables['tb10v'][:]
    Tb19h0 = nc_fid.variables['tb19h'][:]
    Tb19v0 = nc_fid.variables['tb19v'][:]
    Tb23h0 = nc_fid.variables['tb22v'][:] 
    Tb23v0 = nc_fid.variables['tb22v'][:] # replaced by 22v (not 22h)
    Tb37h0 = nc_fid.variables['tb37h'][:]
    Tb37v0 = nc_fid.variables['tb37v'][:]
    Tb85h0 = nc_fid.variables['tb90h'][:] 
    Tb85v0 = nc_fid.variables['tb90v'][:] 

    # These are masked arrays (same mask):
    mask = Tb6h0.mask
    
    #%% Run loops:
    # Initialize masked arrays
    P_array = np.ma.zeros((len(time),len(xc),len(yc),6))
    P_array.mask = np.tile(mask,6)
    S_array = np.ma.zeros((len(time),len(xc),len(yc),6))
    S_array.mask = np.tile(mask,6)
    cost_array = np.ma.zeros((len(time),len(xc),len(yc)))
    cost_array.mask = mask

    for i in xrange(0,len(xc)):  
        for j in xrange(0,len(yc)):
            if mask[:,i,j]: continue
            
            # Extract brightness temperatures for current pixel
            Tb6h = Tb6h0.data[:,i,j] #6h
            Tb6v = Tb6v0.data[:,i,j] #6v
            Tb10v = Tb10v0.data[:,i,j] #10v
            Tb10h = Tb10h0.data[:,i,j] #10h
            Tb19v = Tb19v0.data[:,i,j] #19v
            Tb19h = Tb19h0.data[:,i,j] #19v
            Tb23v = Tb23v0.data[:,i,j] #23v
            Tb23h = Tb23h0.data[:,i,j] #23v
            Tb37v = Tb37v0.data[:,i,j] #37v
            Tb37h = Tb37h0.data[:,i,j] #37v
            Tb89v = Tb85v0.data[:,i,j] #37v!
            Tb89h = Tb85h0.data[:,i,j] #37n!
            # Brightness temperature vector at current location
            Tb = np.array([Tb6v,Tb6h,Tb10v,Tb10h,Tb19v,Tb19h,Tb23v,Tb23h,Tb37v,Tb37h,Tb89v,Tb89h], dtype='float32')
            
            # Set input values for parameters:
            inputpar=sett_input.setinputpar(Tb)
            
            # Optimize physical parameters using optimal estimation:
            P_final, S_diag, cost = optimalest.optimize(Tb,inputpar)
            
            # Add to array:
            P_array[0,i,j,:]=P_final
            S_array[0,i,j,:]=S_diag
            cost_array[0,i,j]=cost            
 
    #%% Save as netcdf file:
    w_nc_fid = Dataset(nc_f_out,'w',format='NETCDF4')
    w_nc_fid.description = "Calculated values of P, S and cost"

    data = {}
    for dim in nc_dims:
        # Create dimensions: Same size as input file
        w_nc_fid.createDimension(dim,nc_fid.variables[dim].size) 
        # Create the variable (datatype and dimensions): Same as input
        data[dim] = w_nc_fid.createVariable(dim,nc_fid.variables[dim].dtype,(dim,))
        # Set attributes (metadata)
        for ncattr in nc_fid.variables[dim].ncattrs():  
            data[dim].setncattr(ncattr,nc_fid.variables[dim].getncattr(ncattr))
   
    #%% Assign the dimension data to the new NetCDF file.
    w_nc_fid.variables['time'][:] = time
    w_nc_fid.variables['xc'][:] = xc
    w_nc_fid.variables['yc'][:] = yc
            
    #%% Create coordinate system (a variable; 2D array) 
    # Transform xc, yc to lat, lons:
    in_proj = Proj('+datum=WGS84 +ellps=WGS84 +lat_0=90.0 +lon_0=0 +proj=laea') #Proj(init='epsg:3575') # Laea, with northpole and greenwich
    out_proj = Proj(init='epsg:4326') # WGS - geographical coordinates
    k = np.tile(xc,(len(yc),1))*10**3
    m = np.array(np.tile(np.matrix(yc).T,(1,len(xc))))*10**3
    lon,lat = transform(in_proj, out_proj, k,m)
    
    nc_grid_lat = w_nc_fid.createVariable('lat','f8',('xc','yc'))
    nc_grid_lon = w_nc_fid.createVariable('lon','f8',('xc','yc'))
    # Add attributes for compliancy:
    nc_grid_lat.setncatts({'long_name': u"latitude",'units': u"degrees_north", 'standard_name':u'latitude'}) 
    nc_grid_lon.setncatts({'long_name': u"longitude",'units': u"degrees_east", 'standard_name':u'longitude'}) 
       
    #nc_grid_lat.long_name = 'latitude coordinate'
    #nc_grid_lat.standard_name = 'latitude'
    #nc_grid_lat.units = 'degrees_north'
        
    w_nc_fid.variables['lat'][:]= lat
    w_nc_fid.variables['lon'][:]= lon
    
    #%% Create departure variables:
    w_nc_var0 = w_nc_fid.createVariable('Ts','f8',('time','xc','yc')) #64-bit floating
    # Add attributes:
    w_nc_var0.setncatts({'long_name': u"Surface temperature [K]"}) 
    # Set coordinates:
    w_nc_var0.coordinates = 'lat lon'

    w_nc_var1 = w_nc_fid.createVariable('V','f8',('time','xc','yc')) 
    w_nc_var1.setncatts({'long_name': u"Water vapor"}) 
    w_nc_var1.coordinates = 'lat lon'

    w_nc_var2 = w_nc_fid.createVariable('W','f8',('time','xc','yc')) 
    w_nc_var2.setncatts({'long_name': u"Wind speed"}) 
    w_nc_var2.coordinates = 'lat lon'
    w_nc_var3 = w_nc_fid.createVariable('L','f8',('time','xc','yc')) 
    w_nc_var3.setncatts({'long_name': u"Integrated liquid water content"}) 
    w_nc_var3.coordinates = 'lat lon'
    w_nc_var4 = w_nc_fid.createVariable('c_ice','f8',('time','xc','yc')) 
    w_nc_var4.setncatts({'long_name': u"Ice concentration"}) 
    w_nc_var4.coordinates = 'lat lon'
    w_nc_var5 = w_nc_fid.createVariable('sd','f8',('time','xc','yc')) 
    w_nc_var5.setncatts({'long_name': u"Snow depth"}) 
    w_nc_var5.coordinates = 'lat lon'

    w_nc_var6 = w_nc_fid.createVariable('cost','f8',('time','xc','yc')) #64-bit floating
    w_nc_var6.setncatts({'long_name': u"cost"})
    w_nc_var6.coordinates = 'lat lon'
    
    w_nc_varS0 = w_nc_fid.createVariable('S_Ts','f8',('time','xc','yc')) #64-bit floating
    w_nc_varS0.setncatts({'long_name': u"S(Ts)"})
    w_nc_varS0.coordinates = 'lat lon'
    w_nc_varS1 = w_nc_fid.createVariable('S_V','f8',('time','xc','yc')) 
    w_nc_varS1.setncatts({'long_name': u"S(V)"})
    w_nc_varS1.coordinates = 'lat lon'
    w_nc_varS2 = w_nc_fid.createVariable('S_W','f8',('time','xc','yc')) 
    w_nc_varS2.setncatts({'long_name': u"S(W)"})
    w_nc_varS2.coordinates = 'lat lon'
    w_nc_varS3 = w_nc_fid.createVariable('S_L','f8',('time','xc','yc')) 
    w_nc_varS3.setncatts({'long_name': u"S(L)"})
    w_nc_varS3.coordinates = 'lat lon'
    w_nc_varS4 = w_nc_fid.createVariable('S_cice','f8',('time','xc','yc')) 
    w_nc_varS4.setncatts({'long_name': u"S(c_ice)"})
    w_nc_varS4.coordinates = 'lat lon'
    w_nc_varS5 = w_nc_fid.createVariable('S_sd','f8',('time','xc','yc')) 
    w_nc_varS5.setncatts({'long_name': u"S(sd)"})
    w_nc_varS5.coordinates = 'lat lon'

    w_nc_fid.variables['cost'][:]= cost_array[0,:,:]
    w_nc_fid.variables['Ts'][:]= P_array[0,:,:,0]
    w_nc_fid.variables['V'][:]= P_array[0,:,:,1]
    w_nc_fid.variables['W'][:]= P_array[0,:,:,2]
    w_nc_fid.variables['L'][:]= P_array[0,:,:,3]
    w_nc_fid.variables['c_ice'][:]= P_array[0,:,:,4]
    w_nc_fid.variables['sd'][:]= P_array[0,:,:,5]

    w_nc_fid.variables['S_Ts'][:]= S_array[0,:,:,0]
    w_nc_fid.variables['S_V'][:]= S_array[0,:,:,1]
    w_nc_fid.variables['S_W'][:]= S_array[0,:,:,2]
    w_nc_fid.variables['S_L'][:]= S_array[0,:,:,3]
    w_nc_fid.variables['S_cice'][:]= S_array[0,:,:,4]
    w_nc_fid.variables['S_sd'][:]= S_array[0,:,:,5]
   
    #%% Close NetCDF file:
    w_nc_fid.close()
    
#%% Run "runows" for a series of netCDF files: 
year = 2013
month = 3
res = 50

for day in range(21,22): #28):
    if day < 10:
        datename = str(year) + '0' + str(month) + '0' + str(day)
    else:
        datename = str(year) + '0' + str(month) + str(day)
    
    filename_in = './Data/tc_amsr2-gw1_nh_ease2-' + str(res) + '0_' + datename + '120000.nc'
    filename_out = './Out/output_' + str(res) + '_' + datename + '_test.nc'
    
    # Run for map of Tb obtained this day:
    runows(filename_in,filename_out)
    

## Previous version:
#year = 2014
#month = 4
#for day in range(1,2):
#    if day < 10:
#        datename = str(year) + '0' + str(month) + '0' + str(day)
#    else:
#        datename = str(year) + '0' + str(month) + str(day)
#    
#    filename_in = './Data/tc_amsr2_nh_ease2-250_' + datename + '120000.nc'
#    filename_out = './Out/output_' + datename + '.nc'
#    
#    # Run for this day:
#    runows(filename_in,filename_out)