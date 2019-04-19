"""treefrac_diff.py

author: Auke Visser
date: 11.10.2016

This script calculates the difference in tree fraction between 
two user-specified periods for a user-specified model.

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

def treefrac_diff(model,yr_start1,yr_end1,yr_start2,yr_end2):
  #Access land fraction files, compute number of simulations and define variables
  lf_file = nc.Dataset('/net/atmos/data/cmip5/historical/fx/sftlf/%s/r0i0p0/sftlf_fx_%s_historical_r0i0p0.nc'%(model,model),'r')
  landfrac = lf_file.variables['sftlf'][:]
  num_sim = len(os.walk('/net/firebolt/data/vissera/Deforestation_scripts_Quentin/CMIP5_data_preprocessed/historical/Amon/tas/%s/'%model).next()[1])
  tf_mean_output = np.zeros((num_sim,landfrac.shape[0],landfrac.shape[1]))
  
  #Define time range for periods
  PI_tr = yr_end1 - yr_start1 + 1
  PD_tr = yr_end2 - yr_start2 + 1
  
  for i in range(0,num_sim):
    if os.path.isdir('/net/firebolt/data/vissera/Deforestation_scripts_Quentin/CMIP5_data_preprocessed/historical/Lmon/treeFrac/%s/%s/'%(model,sim_list[i])) == False:
      continue
    else:
      #Access tree fraction from the simulation of interest
      if(model == 'CanESM2' or model == 'CCSM4'):
	path = '/net/firebolt/data/vissera/Deforestation_scripts_Quentin/CMIP5_data_preprocessed/historical/Lmon/treeFrac/%s/%s/treeFrac_yr_%s_historical_%s_1861-2004.nc'%(model,sim_list[i],model,sim_list[i])
      else:
	path = '/net/firebolt/data/vissera/Deforestation_scripts_Quentin/CMIP5_data_preprocessed/historical/Lmon/treeFrac/%s/%s/treeFrac_Lmon_%s_historical_%s_186101-200412.nc'%(model,sim_list[i],model,sim_list[i])
      ds = nc.Dataset(path,'r')
      
      if(model == 'CanESM2'):
	treeFrac = ds.variables['FRAC'][:]
      elif(model == 'CCSM4'):
	treeFrac = ds.variables['PCT_PFT'][:]
      else:
	treeFrac = ds.variables['treeFrac'][:]
      
      i_s1 = (yr_start1 - 1861) * 12
      i_e1 = (yr_end1 - 1861 + 1) * 12
      i_s2 = (yr_start2 - 1861) * 12
      i_e2 = (yr_end2 - 1861 + 1) * 12
      
      #Calculate mean tree fraction for PI and PD
      tf_PI = np.nanmean(treeFrac[i_s1:i_e1,:,:],axis=0)
      tf_PD = np.nanmean(treeFrac[i_s2:i_e2,:,:],axis=0)
      
      tf_mean_output[i,:,:] = tf_PD - tf_PI
      tf_mean_output[i,:,:][landfrac < 50.] = np.nan #Get rid of non-land cells
      tf_mean_output[i,0:int(landfrac.shape[0]/6.),:] = np.nan #Get rid of Antarctica
  
  return tf_mean_output

#tf_out = treefrac_diff('CanESM2',1861,1890,1975,2004) #Does not work (11.10,2:30pm)
#tf_out = treefrac_diff('CCSM4',1861,1890,1975,2004) #Does not work (11.10,2:30pm)
#tf_out = treefrac_diff('HadGEM2-ES',1861,1890,1975,2004)
#tf_out = treefrac_diff('IPSL-CM5A-LR',1861,1890,1975,2004)
#tf_out = treefrac_diff('NORESM1-M',1861,1890,1975,2004) #Does not work (11.10,2:30pm)
#tf_out = treefrac_diff('GFDL-CM3',1861,1890,1975,2004)
#tf_out = treefrac_diff('IPSL-CM5A-MR',1861,1890,1975,2004)
#tf_out = treefrac_diff('MPI-ESM-LR',1861,1890,1975,2004)
#tf_out = treefrac_diff('MPI-ESM-MR',1861,1890,1975,2004)
