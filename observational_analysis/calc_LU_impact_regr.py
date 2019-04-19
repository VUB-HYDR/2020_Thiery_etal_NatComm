"""calc_LU_impact_regr.py

author: Auke Visser
date: 19.10.2016

This script calculates the land use impact on temperature
following the algorithm by Kumar et al. (2013) and 
Lejeune et al. (2016, in rev.)

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
from sklearn import linear_model

execfile('T_max_seas.py')
execfile('treefrac_diff.py')

sim_list = ['r1i1p1','r2i1p1','r3i1p1','r4i1p1','r5i1p1','r6i1p1']
seas_list = ['DJF','MAM','JJA','SON']
thres_tf = -15.

def calc_LU_impact_regr(model,exp,seas,percentile,yr_start1,yr_end1,yr_start2,yr_end2):
  #Run other programs needed in the calculation
  Tmean_out,Tmax_out = T_max_seas(model,percentile,yr_start1,yr_end1,yr_start2,yr_end2)
  tf_out = treefrac_diff(model,yr_start1,yr_end1,yr_start2,yr_end2)
  
  #Find the indices of the season and experiment of interest, and access the variables of interest
  seas_ind = seas_list.index(seas)
  exp_ind = sim_list.index(exp)
  Tmean_diff = Tmean_out[exp_ind,1,seas_ind,:,:] - Tmean_out[exp_ind,0,seas_ind,:,:]
  Tmax_diff = Tmax_out[exp_ind,1,seas_ind,:,:] - Tmax_out[exp_ind,0,seas_ind,:,:]
  tf_diff = tf_out[exp_ind,:,:]
  
  #Initialize and access variables
  LU_impact_mean = np.zeros((tf_diff.shape[0],tf_diff.shape[1]))
  LU_impact_max = np.zeros((tf_diff.shape[0],tf_diff.shape[1]))
  
  #Access orography and make lat/lon files
  or_file = nc.Dataset('/net/atmos/data/cmip5/historical/fx/orog/%s/r0i0p0/orog_fx_%s_historical_r0i0p0.nc'%(model,model),'r')
  orog = or_file.variables['orog'][:]
  lat = or_file.variables['lat'][:]
  lon = or_file.variables['lon'][:]
  lats = np.transpose(np.tile(lat,(lon.shape[0],1)))
  lons = np.tile(lon,(lat.shape[0],1))
  
  LU_impact_mean[:,:] = np.nan
  LU_impact_max[:,:] = np.nan
  
  for i in range(2,tf_diff.shape[0]-2):
    for j in range(2,tf_diff.shape[1]-2):
      
      if(np.count_nonzero(~np.isnan(tf_diff[i-2:i+3,j-2:j+3])) >= 15 and 
	 np.count_nonzero(tf_diff[i-2:i+2,j-2:j+2] != 0) >= 2 and
	 tf_diff[i,j] <= thres_tf # Added to reproduce Fig. 1 from QL16b 
	 ):
	
	# Find the defrate, lat, lon, and elev. and reshape them into 1D arrays, using
	# only land-based cells (non-land cells are nans in tf_diff)
	tf_diff_box = tf_diff[i-2:i+3,j-2:j+3]
	tf_box = tf_diff_box.ravel()[~np.isnan(tf_diff_box.ravel())]
	dTmean_box = Tmean_diff[i-2:i+3,j-2:j+3].ravel()[~np.isnan(tf_diff_box.ravel())]
	dTmax_box = Tmax_diff[i-2:i+3,j-2:j+3].ravel()[~np.isnan(tf_diff_box.ravel())]
	
	arr_eff = np.zeros((tf_box.shape[0],4))
	arr_eff[:,0] = tf_box
	arr_eff[:,1] = orog[i-2:i+3,j-2:j+3].ravel()[~np.isnan(tf_diff_box.ravel())]
	arr_eff[:,2] = lats[i-2:i+3,j-2:j+3].ravel()[~np.isnan(tf_diff_box.ravel())]
	arr_eff[:,3] = lons[i-2:i+3,j-2:j+3].ravel()[~np.isnan(tf_diff_box.ravel())]
	
	# Define the multi-linear regression model
	regmean = linear_model.LinearRegression()
	regmean.fit(arr_eff,dTmean_box)
	
	regmax = linear_model.LinearRegression()
	regmax.fit(arr_eff,dTmax_box)
	
	# Extract the regression coefficient for deforestation and multiply it by the 
	# defrate in the cell of interest
	LU_impact_mean[i,j] = regmean.coef_[0] * tf_diff[i,j]
	LU_impact_max[i,j] = regmax.coef_[0] * tf_diff[i,j]
  
  return(LU_impact_mean,LU_impact_max)

modelname = 'MPI-ESM-MR'

if modelname == 'IPSL-CM5A-LR':
  dTmean_def = np.zeros((6,96,96))
  dTmax_def = np.zeros((6,96,96))
elif modelname == 'IPSL-CM5A-MR':
  dTmean_def = np.zeros((3,143,144))
  dTmax_def = np.zeros((3,143,144))
elif(modelname == 'MPI-ESM-LR' or modelname == 'MPI-ESM-MR'):
  dTmean_def = np.zeros((3,96,192))
  dTmax_def = np.zeros((3,96,192))

for k in range(0,dTmean_def.shape[0]):
  dTmean_def[k,:,:], dTmax_def[k,:,:] = calc_LU_impact_regr(modelname,sim_list[k],seas_list[2],99,1861,1890,1975,2004)