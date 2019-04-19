"""T_mean_seas.py

author: Auke Visser
date: 10.10.2016

This script calculates the seasonal T mean for pre-industrial temperature (PI)
and present-day (PD) climate (i.e. the 1st and last 30 year block of the gridded
temperature dataset)

Output structure: array of dimensions [num_sim,times,seas,lat,lon]
    where num_sim is the number of ensemble simulations, times = 2 (PI & PD),
    seas = 4 (DJF, MAM, JJA, SON), and lat,lon are model-dependent (96 for now)
    
This code is inspired by Quentin Lejeune's NCL version of the
Kumar algorithm.
    
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

sim_list = ['r1i1p1','r2i1p1','r3i1p1','r4i1p1','r5i1p1','r6i1p1']

def T_mean_seas(model,yr_start1,yr_end1,yr_start2,yr_end2):
  lf_file = nc.Dataset('/net/atmos/data/cmip5/historical/fx/sftlf/%s/r0i0p0/sftlf_fx_%s_historical_r0i0p0.nc'%(model,model),'r')
  landfrac = lf_file.variables['sftlf'][:]
  n_lat = landfrac.shape[0]
  n_lon = landfrac.shape[1]
  num_sim = len(os.walk('/net/firebolt/data/vissera/Deforestation_scripts_Quentin/CMIP5_data_preprocessed/historical/Amon/tas/%s/'%model).next()[1])
  T_mean_output = np.zeros((num_sim,2,4,n_lat,n_lon))
  
  PI_tr = yr_end1 - yr_start1 + 1
  PD_tr = yr_end2 - yr_start2 + 1
  
  for i in range(0,num_sim):
    if os.path.isdir('/net/firebolt/data/vissera/Deforestation_scripts_Quentin/CMIP5_data_preprocessed/historical/Amon/tas/%s/%s/'%(model,sim_list[i])) == False:
      continue
    else:
      #Declare path and other variables
      path = '/net/firebolt/data/vissera/Deforestation_scripts_Quentin/CMIP5_data_preprocessed/historical/Amon/tas/%s/%s/tas_Amon_%s_historical_%s_186101-200412.nc'%(model,sim_list[i],model,sim_list[i])
      ds = nc.Dataset(path,'r')
      tas = ds.variables['tas'][:]

      PI_DJF = np.zeros((PI_tr*3,n_lat,n_lon))
      PI_MAM = np.zeros((PI_tr*3,n_lat,n_lon))
      PI_JJA = np.zeros((PI_tr*3,n_lat,n_lon))
      PI_SON = np.zeros((PI_tr*3,n_lat,n_lon))

      PD_DJF = np.zeros((PD_tr*3,n_lat,n_lon))
      PD_MAM = np.zeros((PD_tr*3,n_lat,n_lon))
      PD_JJA = np.zeros((PD_tr*3,n_lat,n_lon))
      PD_SON = np.zeros((PD_tr*3,n_lat,n_lon))

      i_e = (yr_end1 - 1861 + 1) * 12
      i_s = (yr_start1 - 1861) * 12

      #Access PI monthly temperature data and group in seasonal arrays
      PI_DJF[0:PI_tr*3:3,:,:] = tas[i_s:i_e:12,:,:]
      PI_DJF[1:PI_tr*3:3,:,:] = tas[i_s+1:i_e:12,:,:]
      PI_DJF[2:PI_tr*3:3,:,:] = tas[i_s+11:i_e:12,:,:]

      PI_MAM[0:PI_tr*3:3,:,:] = tas[i_s+2:i_e:12,:,:]
      PI_MAM[1:PI_tr*3:3,:,:] = tas[i_s+3:i_e:12,:,:]
      PI_MAM[2:PI_tr*3:3,:,:] = tas[i_s+4:i_e:12,:,:]

      PI_JJA[0:PI_tr*3:3,:,:] = tas[i_s+5:i_e:12,:,:]
      PI_JJA[1:PI_tr*3:3,:,:] = tas[i_s+6:i_e:12,:,:]
      PI_JJA[2:PI_tr*3:3,:,:] = tas[i_s+7:i_e:12,:,:]

      PI_SON[0:PI_tr*3:3,:,:] = tas[i_s+8:i_e:12,:,:]
      PI_SON[1:PI_tr*3:3,:,:] = tas[i_s+9:i_e:12,:,:]
      PI_SON[2:PI_tr*3:3,:,:] = tas[i_s+10:i_e:12,:,:]

      #Calculate PI seasonal mean temperature
      T_mean_output[i,0,0,:,:] = np.nanmean(PI_DJF,axis=0)
      T_mean_output[i,0,1,:,:] = np.nanmean(PI_MAM,axis=0)
      T_mean_output[i,0,2,:,:] = np.nanmean(PI_JJA,axis=0)
      T_mean_output[i,0,3,:,:] = np.nanmean(PI_SON,axis=0)
      
      #Repeat this procedure for PD
      #Access PD monthly temperature data and group in seasonal arrays
      i_e = (yr_end2 - 1861 + 1) * 12
      i_s = (yr_start2 - 1861) * 12
      
      PD_DJF[0:PD_tr*3:3,:,:] = tas[i_s:i_e:12,:,:]
      PD_DJF[1:PD_tr*3:3,:,:] = tas[i_s+1:i_e:12,:,:]
      PD_DJF[2:PD_tr*3:3,:,:] = tas[i_s+11:i_e:12,:,:]

      PD_MAM[0:PD_tr*3:3,:,:] = tas[i_s+2:i_e:12,:,:]
      PD_MAM[1:PD_tr*3:3,:,:] = tas[i_s+3:i_e:12,:,:]
      PD_MAM[2:PD_tr*3:3,:,:] = tas[i_s+4:i_e:12,:,:]

      PD_JJA[0:PD_tr*3:3,:,:] = tas[i_s+5:i_e:12,:,:]
      PD_JJA[1:PD_tr*3:3,:,:] = tas[i_s+6:i_e:12,:,:]
      PD_JJA[2:PD_tr*3:3,:,:] = tas[i_s+7:i_e:12,:,:]

      PD_SON[0:PD_tr*3:3,:,:] = tas[i_s+8:i_e:12,:,:]
      PD_SON[1:PD_tr*3:3,:,:] = tas[i_s+9:i_e:12,:,:]
      PD_SON[2:PD_tr*3:3,:,:] = tas[i_s+10:i_e:12,:,:]

      #Calculate PD seasonal mean temperature
      T_mean_output[i,1,0,:,:] = np.nanmean(PD_DJF,axis=0)
      T_mean_output[i,1,1,:,:] = np.nanmean(PD_MAM,axis=0)
      T_mean_output[i,1,2,:,:] = np.nanmean(PD_JJA,axis=0)
      T_mean_output[i,1,3,:,:] = np.nanmean(PD_SON,axis=0)
      
      for j in range(0,2):
	for k in range(0,4):
	  T_mean_output[i,j,k,:,:][landfrac < 50.] = np.nan #Get rid of non-land cells
	  T_mean_output[i,j,k,0:int(n_lat/6.),:]= np.nan #Get rid of Antarctica
      
      del PI_DJF,PI_MAM,PI_JJA,PI_SON
      del PD_DJF,PD_MAM,PD_JJA,PD_SON
      
  return T_mean_output

#y = T_mean_seas('CanESM2',1861,1870,1995,2004) #Does not work
#y = T_mean_seas('CCSM4',1861,1870,1995,2004) #Does not work
#y = T_mean_seas('HadGEM2-ES',1861,1890,1975,2004)
#y = T_mean_seas('IPSL-CM5A-LR',1861,1870,1995,2004)
#y = T_mean_seas('NorESM1-M',1861,1870,1995,2004) #Does not work
#y = T_mean_seas('GFDL-CM3',1861,1870,1995,2004)
#y = T_mean_seas('IPSL-CM5A-MR',1861,1870,1995,2004)
y = T_mean_seas('MPI-ESM-LR',1861,1870,1995,2004)
#y = T_mean_seas('MPI-ESM-MR',1861,1890,1975,2004)
