"""extract_T_irr.py

author: Auke Visser
date: 24.10.2016

This script extracts a user-specified temperature product (mean T, TXx, DTR) from a user-specified
data source (CRU, E-OBS, etc.), and calculates either a transient temperature response or a 
PC/PD (previous century -> present-day) response

First goal: implement PC/PD response for CRU monthly mean temperature

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

def extract_T_irr(datasource,temp_product,response,t_res,yr_start1,yr_end1,yr_start2,yr_end2):
  
  #Define time range
  PC_tr = yr_end1 - yr_start1 + 1
  PD_tr = yr_end2 - yr_start2 + 1
  
  if(datasource in ['CRU','CRUv4_CESM','CRU_CESM'] and temp_product in ['tmp_max', 'dtr_max', 'tmx_max', 'tmn_min','tmn_max']):
    
    #Load netCDF files
    if datasource == 'CRU':
      filepath = '/net/exo/landclim/wthiery/observational_analysis/Data/cru_ts3.22.1901.2013.%s.dat.nc'%temp_product
    elif datasource == 'CRU_CESM':
      filepath = '/net/exo/landclim/wthiery/observational_analysis/Data/cru_ts3.22.1901.2013.0.9375x1.25.%s.dat.nc'%temp_product
    elif datasource == 'CRUv4_CESM':
      filepath = '/home/vissera/Documents/Data/cru_ts4.02.1901.2017.0.9375x1.25.%s.dat.nc'%temp_product
    Tfile = nc.Dataset(filepath,'r')
    n_lat = Tfile.variables['lat'][:].shape[0]
    n_lon = Tfile.variables['lon'][:].shape[0]
    
    #Define start and end indices
    i_s1 = (yr_start1 - 1901)
    i_e1 = (yr_end1 - 1901 + 1)
    i_s2 = (yr_start2 - 1901)
    i_e2 = (yr_end2 - 1901 + 1)
    
    #Define array in which 30-year averaged temperature values are stored
    #for reference and present-day period
    T_output = np.zeros((2,n_lat,n_lon))
    
    #Extract the variable from the input netCDF file
    if temp_product == 'tmp_max':
      T_prod = Tfile.variables['tmp'][:]
    elif temp_product == 'dtr_max':
      T_prod = Tfile.variables['dtr'][:]
    elif temp_product == 'tmx_max':
      T_prod = Tfile.variables['tmx'][:]
    elif temp_product in ['tmn_min','tmn_max']:
      T_prod = Tfile.variables['tmn'][:]
    
    #Calculate average values for PC and PD periods and 
    #store them in the output array
    T_output[0,:,:] = np.nanmean(T_prod[i_s1:i_e1,:,:],axis=0)
    T_output[1,:,:] = np.nanmean(T_prod[i_s2:i_e2,:,:],axis=0)
    
    del T_prod
    
    #Extract HID file to obtain land mask
    if datasource == 'CRU':
      HIDfile = nc.Dataset('/net/exo/landclim/wthiery/observational_analysis/Data/HID_regridded/hid_v1.0.f_irr_hyde_final_ir.0.5deg.10y.1900-2005.nc','r')
    elif datasource in ['CRU_CESM','CRUv4_CESM']:
      HIDfile = nc.Dataset('/net/exo/landclim/wthiery/observational_analysis/Data/HID_regridded/hid_v1.0.f_irr_hyde_final_ir.0.9375x1.25deg.10y.1900-2005.nc','r')
    
    aei = HIDfile.variables['aei_hyde_final_ir'][:]
    
    for i in range(0,2):
      T_output[i,:,:][aei.mask[0,:,:] == True] = np.nan
    
    #Swap the hemispheres for CRU regridded CESM resolution
    if datasource in ['CRU_CESM','CRUv4_CESM']:
      WH = T_output[:,:,144:288]
      EH = T_output[:,:,0:144]
      T_output = np.zeros((2,n_lat,n_lon))
      T_output[:,:,0:144] = WH
      T_output[:,:,144:288] = EH
    
    del aei
  
  elif(datasource == 'CRU' and temp_product in ['tmp', 'dtr', 'tmx', 'tmn']):
    filepath = '/net/exo/landclim/data/dataset/CRUTS/v3.22/0.5deg_lat-lon_1m/original/cru_ts3.22.1901.2013.%s.dat.nc'%temp_product
    Tfile = nc.Dataset(filepath,'r')
    n_lat = Tfile.variables['lat'][:].shape[0]
    n_lon = Tfile.variables['lon'][:].shape[0]

    #Define start and end indices
    i_s1 = (yr_start1 - 1901) * 12
    i_e1 = (yr_end1 - 1901 + 1) * 12
    i_s2 = (yr_start2 - 1901) * 12
    i_e2 = (yr_end2 - 1901 + 1) * 12
    
    if t_res == 'seasonal':
      #Define variables
      T_output = np.zeros((2,4,n_lat,n_lon))
      
      DJF_PC = np.zeros((PC_tr*3,n_lat,n_lon))
      MAM_PC = np.zeros((PC_tr*3,n_lat,n_lon))
      JJA_PC = np.zeros((PC_tr*3,n_lat,n_lon))
      SON_PC = np.zeros((PC_tr*3,n_lat,n_lon))
      
      DJF_PD = np.zeros((PD_tr*3,n_lat,n_lon))
      MAM_PD = np.zeros((PD_tr*3,n_lat,n_lon))
      JJA_PD = np.zeros((PD_tr*3,n_lat,n_lon))
      SON_PD = np.zeros((PD_tr*3,n_lat,n_lon))
      
      T_prod = Tfile.variables['%s'%temp_product][:]
      
      T_PC = T_prod[i_s1:i_e1,:,:]
      T_PD = T_prod[i_s2:i_e2,:,:]
      del T_prod
      
      #Access PC monthly temperature and group by season
      DJF_PC[0:PC_tr*3:3,:,:] = T_PC[0:i_e1-i_s1:12,:,:]
      DJF_PC[1:PC_tr*3:3,:,:] = T_PC[1:i_e1-i_s1:12,:,:]
      DJF_PC[2:PC_tr*3:3,:,:] = T_PC[11:i_e1-i_s1:12,:,:]
      
      MAM_PC[0:PC_tr*3:3,:,:] = T_PC[2:i_e1-i_s1:12,:,:]
      MAM_PC[1:PC_tr*3:3,:,:] = T_PC[3:i_e1-i_s1:12,:,:]
      MAM_PC[2:PC_tr*3:3,:,:] = T_PC[4:i_e1-i_s1:12,:,:]

      JJA_PC[0:PC_tr*3:3,:,:] = T_PC[5:i_e1-i_s1:12,:,:]
      JJA_PC[1:PC_tr*3:3,:,:] = T_PC[6:i_e1-i_s1:12,:,:]
      JJA_PC[2:PC_tr*3:3,:,:] = T_PC[7:i_e1-i_s1:12,:,:]

      SON_PC[0:PC_tr*3:3,:,:] = T_PC[8:i_e1-i_s1:12,:,:]
      SON_PC[1:PC_tr*3:3,:,:] = T_PC[9:i_e1-i_s1:12,:,:]
      SON_PC[2:PC_tr*3:3,:,:] = T_PC[10:i_e1-i_s1:12,:,:]
      
      #Repeat this procedure for present-day
      DJF_PD[0:PC_tr*3:3,:,:] = T_PD[0:i_e2-i_s2:12,:,:]
      DJF_PD[1:PC_tr*3:3,:,:] = T_PD[1:i_e2-i_s2:12,:,:]
      DJF_PD[2:PC_tr*3:3,:,:] = T_PD[11:i_e2-i_s2:12,:,:]
      
      MAM_PD[0:PC_tr*3:3,:,:] = T_PD[2:i_e2-i_s2:12,:,:]
      MAM_PD[1:PC_tr*3:3,:,:] = T_PD[3:i_e2-i_s2:12,:,:]
      MAM_PD[2:PC_tr*3:3,:,:] = T_PD[4:i_e2-i_s2:12,:,:]

      JJA_PD[0:PC_tr*3:3,:,:] = T_PD[5:i_e2-i_s2:12,:,:]
      JJA_PD[1:PC_tr*3:3,:,:] = T_PD[6:i_e2-i_s2:12,:,:]
      JJA_PD[2:PC_tr*3:3,:,:] = T_PD[7:i_e2-i_s2:12,:,:]

      SON_PD[0:PC_tr*3:3,:,:] = T_PD[8:i_e2-i_s2:12,:,:]
      SON_PD[1:PC_tr*3:3,:,:] = T_PD[9:i_e2-i_s2:12,:,:]
      SON_PD[2:PC_tr*3:3,:,:] = T_PD[10:i_e2-i_s2:12,:,:]
      
      #Calculate PC seasonal mean temperature      
      T_output[0,0,:,:] = np.nanmean(DJF_PC,axis=0)
      T_output[0,1,:,:] = np.nanmean(MAM_PC,axis=0)
      T_output[0,2,:,:] = np.nanmean(JJA_PC,axis=0)
      T_output[0,3,:,:] = np.nanmean(SON_PC,axis=0)
      
      #Calculate PD seasonal mean temperature
      T_output[1,0,:,:] = np.nanmean(DJF_PD,axis=0)
      T_output[1,1,:,:] = np.nanmean(MAM_PD,axis=0)
      T_output[1,2,:,:] = np.nanmean(JJA_PD,axis=0)
      T_output[1,3,:,:] = np.nanmean(SON_PD,axis=0)
      
      #Extract HID file to obtain land mask
      HIDfile = nc.Dataset('/net/exo/landclim/wthiery/observational_analysis/Data/HID_regridded/hid_v1.0.f_irr_hyde_final_ir.0.5deg.10y.1900-2005.nc','r')
      aei = HIDfile.variables['aei_hyde_final_ir'][:]
      
      for i in range(0,2):
        for j in range(0,4):
          T_output[i,j,:,:][aei.mask[0,:,:] == True] = np.nan
      
      del DJF_PC,MAM_PC,JJA_PC,SON_PC
      del DJF_PD,MAM_PD,JJA_PD,SON_PD
      del aei
      
    elif t_res == 'monthly':
      #Define variables
      T_output = np.zeros((2,12,n_lat,n_lon))
      
      T_prod = Tfile.variables['%s'%temp_product][:]
      PC = T_prod[i_s1:i_e1,:,:]
      PD = T_prod[i_s2:i_e2,:,:]
      
      del T_prod
      
      #Calculate 30-year average monthly temperature for PC and PD
      for m in range(0,12):
        T_output[0,m,:,:] = np.nanmean(PC[m::12,:,:],axis=0)
        T_output[1,m,:,:] = np.nanmean(PD[m::12,:,:],axis=0)
	
      #Extract HID file to obtain land mask
      HIDfile = nc.Dataset('/net/exo/landclim/wthiery/observational_analysis/Data/HID_regridded/hid_v1.0.f_irr_hyde_final_ir.0.5deg.10y.1900-2005.nc','r')
      aei = HIDfile.variables['aei_hyde_final_ir'][:]
      
      for i in range(0,2):
        for j in range(0,12):
          T_output[i,j,:,:][aei.mask[0,:,:] == True] = np.nan
      
      del PC,PD,aei
      
  elif datasource == 'Berkeley_CESM':
        
    filepath = '/home/vissera/Documents/Data/Berkeley_Complete_tmx_max_AV_LatLong0.9375x1.25.nc'
    Tfile = nc.Dataset(filepath,'r')
    n_lat = Tfile.variables['lat'].shape[0]
    n_lon = Tfile.variables['lon'].shape[0]
    
    i_s1 = yr_start1 - 1850
    i_e1 = yr_end1 + 1 - 1850
    i_s2 = yr_start2 - 1850
    i_e2 = yr_end2 + 1 - 1850
    
    T_output = np.zeros((2,n_lat,n_lon))
    T_output[:,:,:] = np.nan
    
    T_prod = Tfile.variables['%s'%temp_product][:]
    T_vals = T_prod.data
    
    T_output[0,:,:] = np.nanmean(T_vals[i_s1:i_e1,:,:],axis=0)
    T_output[1,:,:] = np.nanmean(T_vals[i_s2:i_e2,:,:],axis=0)
    
    #Swap the hemispheres for Berkeley regridded to CESM resolution
    WH = T_output[:,:,144:288]
    EH = T_output[:,:,0:144]
    T_output = np.zeros((2,n_lat,n_lon))
    T_output[:,:,0:144] = WH
    T_output[:,:,144:288] = EH
          
  elif datasource == 'E-OBS':
    filepath = '/net/exo/landclim/wthiery/observational_analysis/Data/E-OBS_TXx_0.25deg_reg_v14.0.nc'
    Tfile = nc.Dataset(filepath,'r')
    n_lat = Tfile.variables['latitude'][:].shape[0]
    n_lon = Tfile.variables['longitude'][:].shape[0]
    
    T_output = np.zeros((2,n_lat,n_lon))
    T_output[:,:,:] = np.nan
    
    #Extract TXx data, remove nan values
    T_prod = Tfile.variables['%s'%temp_product][:]
    T_vals = T_prod.data
    T_vals[T_vals == -9999] = np.nan
    
    i_s1 = yr_start1 - 1950
    i_e1 = yr_end1 - 1950 + 1
    i_s2 = yr_start2 - 1950
    i_e2 = yr_end2 - 1950 + 1
    
    #Fill the T_output array, ensuring that there are no data gaps
    for i in range(0,n_lat):
      for j in range(0,n_lon):
        if( np.count_nonzero(~np.isnan(T_vals[i_s1:i_e1,i,j])) == (yr_end1 - yr_start1 + 1) and
            np.count_nonzero(~np.isnan(T_vals[i_s2:i_e2,i,j])) == (yr_end2 - yr_start2 + 1)
          ):
          T_output[0,i,j] = np.nanmean(T_vals[i_s1:i_e1,i,j],axis=0)
          T_output[1,i,j] = np.nanmean(T_vals[i_s2:i_e2,i,j],axis=0)
    
    del T_prod,T_vals
      
  elif datasource == 'HadEX2':
    filepath = '/net/exo/landclim/data/dataset/HadEX2/20150106/2.5x3.75deg_lat-lon_1y/original/H2_TXx_1901-2010_RegularGrid_global_3.75x2.5deg_LSmask.nc'
    Tfile = nc.Dataset(filepath,'r')
    n_lat = Tfile.variables['lat'][:].shape[0]
    n_lon = Tfile.variables['lon'][:].shape[0]
    
    T_output = np.zeros((2,n_lat,n_lon))
    T_output[:,:,:] = np.nan
    
    #Extract TXx data, remove nan values
    T_prod = Tfile.variables['%s'%temp_product][:]
    T_vals = np.zeros((T_prod.shape[0],T_prod.shape[1],T_prod.shape[2]))
    T_vals[:,:,0:48] = T_prod[:,:,48:96]
    T_vals[:,:,48:96] = T_prod[:,:,0:48]
    T_vals[T_vals < -99.8] = np.nan
    
    #Define indices for the two time periods
    i_s1 = yr_start1 - 1901
    i_e1 = yr_end1 - 1901 + 1
    i_s2 = yr_start2 - 1901
    i_e2 = yr_end2 - 1901 + 1
    
    #Fill the T_output array with the means over the time period,
    #ensuring that there are no data gaps
    for i in range(0,n_lat):
      for j in range(0,n_lon):
        if(	np.count_nonzero(~np.isnan(T_vals[i_s1:i_e1,i,j])) >= np.ceil(0.8*(yr_end1 - yr_start1 + 1)) and
          np.count_nonzero(~np.isnan(T_vals[i_s2:i_e2,i,j])) >= np.ceil(0.8*(yr_end1 - yr_start1 + 1))
          ):
          T_output[0,i,j] = np.nanmean(T_vals[i_s1:i_e1,i,j],axis=0)
          T_output[1,i,j] = np.nanmean(T_vals[i_s2:i_e2,i,j],axis=0)
      
    del T_prod,T_vals
    
  elif datasource == 'CESM':
    filepath_PD = '/net/exo/landclim/wthiery/observational_analysis/Data/CESM_output/CESM.WT.1981-2010.%s_irr_ensmonmean_yearmax.nc'%temp_product
#    filepath_PC = '/net/exo/landclim/wthiery/observational_analysis/Data/CESM_output/CESM.WT.1981-2010.%s_20cc_ensmonmean_yearmax.nc'%temp_product
    filepath_PC = '/net/cfc/landclim1/wthiery/cesm_output/ensmean/atm/postprocessed/f.e122.F1850PDC5.f09_g16.20cirr-io192.ensmean_atm_h1_%s_monmean_yearmax.nc'%temp_product
    
    #NOTE: there still seems to be an error in the post-processing of there
    #pic simulation (31/1). Until this is fixed, these data will be used 
    #(with caution) in this analysis. 
    PD_file = nc.Dataset(filepath_PD,'r')
    PC_file = nc.Dataset(filepath_PC,'r')
    n_lat = PD_file.variables['lat'][:]
    n_lon = PD_file.variables['lon'][:]
    
    T_output = np.zeros((2,n_lat.shape[0],n_lon.shape[0]))
    T_output[:,:,:] = np.nan
    
    #Extract annual maximum monthly mean daily maximum T
    PD_prod = PD_file.variables['%s'%temp_product][:]
    if filepath_PC == '/net/exo/landclim/wthiery/observational_analysis/Data/CESM_output/CESM.WT.1981-2010.%s_pic_ensmonmean_yearmax.nc'%temp_product:
      PC_prod = PC_file.variables['%s'%temp_product][:][1:31]
    else:
      PC_prod = PC_file.variables['%s'%temp_product][:]
    
    PD_vals = np.zeros((PD_prod.shape[0],PD_prod.shape[1],PD_prod.shape[2]))
    PC_vals = np.zeros((PC_prod.shape[0],PC_prod.shape[1],PC_prod.shape[2]))
    
    #Swap the hemispheres
    PD_vals[:,:,0:144] = PD_prod[:,:,144:288]
    PD_vals[:,:,144:288] = PD_prod[:,:,0:144]
    PC_vals[:,:,0:144] = PC_prod[:,:,144:288]
    PC_vals[:,:,144:288] = PC_prod[:,:,0:144]
    
    T_output[0,:,:] = np.nanmean(PC_vals,axis=0)
    T_output[1,:,:] = np.nanmean(PD_vals,axis=0)
    
    #Convert from Kelvin to degree Celsius
    T_output -= 273.16

  return T_output

#T_out = extract_T_irr("CRU","tmp","PC/PD","seasonal",1901,1930,1981,2010)
#T_out = extract_T_irr("E-OBS","tx","PC/PD","seasonal",1951,1980,1981,2010)
#T_out = extract_T_irr("CRU","tmx_max","PC/PD","seasonal",1901,1930,1981,2010)
#T_out = extract_T_irr("CRU","tmp","PC/PD","monthly",1901,1930,1981,2010)
#T_out = extract_T_irr("HadEX2","Ann","PC/PD","seasonal",1956,1980,1985,2009)
#T_out = extract_T_irr('CESM','TREFHT','PD','seasonal',1901,1930,1981,2010)
