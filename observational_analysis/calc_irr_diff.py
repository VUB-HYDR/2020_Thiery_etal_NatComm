"""calc_irr_diff.py

author: Auke Visser
date: 25.10.2016

This script calculates the change in irrigated fraction between two 30-year time periods from the (regridded) HID data set


"""

import netCDF4 as nc
import numpy as np
import scipy
import os
import matplotlib
matplotlib.rcParams['backend'] = "Qt4Agg"
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import time


#Define the function
def calc_irr_diff(datasource,response,yr_start1,yr_end1,yr_start2,yr_end2):
  time_list = [1900, 1910, 1920, 1930, 1940, 1950, 1960, 1970, 1980, 1985, 1990, 1995, 2000, 2005]
  if datasource == "CRU":
    #Access variables from the HID file
    HIDfile = nc.Dataset('/net/exo/landclim/wthiery/observational_analysis/Data/HID_regridded/hid_v1.0.f_irr_hyde_final_ir.0.5deg.10y.1900-2005.nc','r')
    f_irr = HIDfile.variables['aei_hyde_final_ir'][:]
    lat = HIDfile.variables['lat'][:]
    lon = HIDfile.variables['lon'][:]
  elif datasource == "E-OBS":
    HIDfile = nc.Dataset('/net/exo/landclim/wthiery/observational_analysis/Data/HID_regridded/hid_v1.0.f_irr_hyde_final_ir.0.25deg.10y.1900-2005.nc','r')
    f_irr = HIDfile.variables['aei_hyde_final_ir'][:]
    lat = HIDfile.variables['latitude'][:]
    lon = HIDfile.variables['longitude'][:]
  elif datasource == "HadEX2":
    HIDfile = nc.Dataset('/net/exo/landclim/wthiery/observational_analysis/Data/HID_regridded/hid_v1.0.f_irr_hyde_final_ir.2.5x3.75deg.10y.1900-2005.nc','r')
    f_irr = HIDfile.variables['aei_hyde_final_ir'][:]
    lat = HIDfile.variables['lat'][:]
    lon = HIDfile.variables['lon'][:]
#  elif datasource == "CESM":
#    HIDfile = nc.Dataset('/net/exo/landclim/wthiery/observational_analysis/Data/HID_regridded/CESM/f.e122.F1850PDC5.f09_g16.irrigation-io192.001_constants.nc','r')
#    f_irr = HIDfile.variables['PCT_IRRIG'][:]
#    lat = HIDfile.variables['lat'][:]
#    lon = HIDfile.variables['lon'][:]
  elif datasource in ["CRU_CESM","CRUv4_CESM","Berkeley_CESM","CESM"]:
    HIDfile = nc.Dataset('/net/exo/landclim/wthiery/observational_analysis/Data/HID_regridded/hid_v1.0.f_irr_hyde_final_ir.0.9375x1.25deg.10y.1900-2005.nc','r')
    f_irr = HIDfile.variables['aei_hyde_final_ir'][:]
    lat = HIDfile.variables['lat'][:]
    lon = HIDfile.variables['lon'][:]

  # Constructed a stacked array 3D that contains the AEI for every year from 
  # 1900-2010 (i.e. 5 years for 1900-1904, 10 years for 1905-1914, ... 8 years for
  # 1975-1982, 5 years from 1983-1987, ..., 8 years from 2003-2010).
  if(datasource in ['CRU','E-OBS','HadEX2','CRU_CESM','CRUv4_CESM','Berkeley_CESM','CESM'] and response in ['PC/PD', 'PD']):
    f_irr_stack = np.concatenate((	np.tile(f_irr[0,:,:],(5,1,1)),    \
                                    np.tile(f_irr[1,:,:],(10,1,1)),   \
                                    np.tile(f_irr[2,:,:],(10,1,1)),   \
                                    np.tile(f_irr[3,:,:],(10,1,1)),   \
                                    np.tile(f_irr[4,:,:],(10,1,1)),   \
                                    np.tile(f_irr[5,:,:],(10,1,1)),   \
                                    np.tile(f_irr[6,:,:],(10,1,1)),   \
                                    np.tile(f_irr[7,:,:],(10,1,1)),   \
                                    np.tile(f_irr[8,:,:],(8,1,1)),    \
                                    np.tile(f_irr[9,:,:],(5,1,1)),   \
                                    np.tile(f_irr[10,:,:],(5,1,1)),  \
                                    np.tile(f_irr[11,:,:],(5,1,1)),  \
                                    np.tile(f_irr[12,:,:],(5,1,1)),  \
                                    np.tile(f_irr[13,:,:],(8,1,1))),axis=0)
    
    # Extract two 30-year periods: the reference period (PC_irr), and
    # the present-day period (PD_irr)
    PC_irr = f_irr_stack[yr_start1-1900:yr_end1-1900+1,:,:]
    PD_irr = f_irr_stack[yr_start2-1900:yr_end2-1900+1,:,:]
    del f_irr_stack
      
    
    #Get rid of -9999 values. 
    PC_irr[PC_irr == -9999.0] = np.nan
    PD_irr[PD_irr == -9999.0] = np.nan
    
    #Calculate the mean f_irr for the PC and PD periods
    PC_irr_timavg = np.nanmean(PC_irr,axis=0)
    PD_irr_timavg = np.nanmean(PD_irr,axis=0)
    
    #Store them in a [2 x n_lat x n_lon] array
    f_irr_diff = np.zeros((2,lat.shape[0],lon.shape[0]))
    f_irr_diff[0,:,:] = PC_irr_timavg.data #the addition .data is to omit the data mask from the masked array
    f_irr_diff[1,:,:] = PD_irr_timavg.data
    
    #CESM and HadEX2 are on a different grid, so the hemispheres must be swapped.
    if datasource in ["CRU_CESM", "CRUv4_CESM", "Berkeley_CESM","CESM"]:
      WH = f_irr_diff[:,:,144:288]
      EH = f_irr_diff[:,:,0:144]
      f_irr_diff = np.zeros((2,192,288))
      f_irr_diff[:,:,0:144] = WH
      f_irr_diff[:,:,144:288] = EH
    
    if datasource == "HadEX2":
      WH = f_irr_diff[:,:,48:96]
      EH = f_irr_diff[:,:,0:48]
      f_irr_diff = np.zeros((2,73,96))
      f_irr_diff[:,:,0:48] = WH
      f_irr_diff[:,:,48:96] = EH
    
    #Delete the variables that are no longer needed
    del PC_irr_timavg,PD_irr_timavg
  
  #For the model runs, the PC period (simulation '20cc') assumes no irrigation,
  #so f_irr needs to be set to 0.
  elif datasource in ['CESM']:
    f_irr_diff = np.zeros((2,lat.shape[0],lon.shape[0]))
    f_irr_diff[1,:,0:144] = f_irr[:,144:288]/100.
    f_irr_diff[1,:,144:288] = f_irr[:,0:144]/100.
    lon -= 180.
  
  return lat,lon,f_irr_diff

#These lines can be used to test the method. For regular execution of the
#algorithms, they should be commented out.

#lat,lon,f_irr = calc_irr_diff('CRU','PC/PD',1901,1930,1981,2010)
#lats,lons,f_irr = calc_irr_diff('E-OBS','PC/PD',1901,1910,2001,2010)
#lats,lons,f_irr = calc_irr_diff('HadEX2','PC/PD',1956,1980,1985,2009)
#lats,lons,f_irr = calc_irr_diff('CRU_CESM','PD',1901,1930,1981,2010)
