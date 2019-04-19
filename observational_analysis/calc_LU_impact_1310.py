"""calc_LU_impact_1310.py

author: Auke Visser
date: 13.10.2016
last edits: 26.10.2016

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

execfile('T_mean_seas.py')
execfile('treefrac_diff.py')

sim_list = ['r1i1p1','r2i1p1','r3i1p1','r4i1p1','r5i1p1','r6i1p1']
seas_list = ['DJF','MAM','JJA','SON']

def calc_LU_impact(model,exp,seas,yr_start1,yr_end1,yr_start2,yr_end2):
  #Run other programs needed in the calculation
  T_out = T_mean_seas(model,yr_start1,yr_end1,yr_start2,yr_end2)
  tf_out = treefrac_diff(model,yr_start1,yr_end1,yr_start2,yr_end2)
  
  #Find the indices of the season and experiment of interest, and access the variables of interest
  seas_ind = seas_list.index(seas)
  exp_ind = sim_list.index(exp)
  T_diff = T_out[exp_ind,1,seas_ind,:,:] - T_out[exp_ind,0,seas_ind,:,:]
  tf_diff = tf_out[exp_ind,:,:]
    
  #Initialize and access variables
  LU_impact = np.zeros((tf_diff.shape[0],tf_diff.shape[1]))
  CtSt = np.zeros((tf_diff.shape[0],tf_diff.shape[1]))
  thres_tf = -15.
  BB = 0
  count=0
  
  LU_impact[:,:] = np.nan
  CtSt[:,:] = np.nan
  
  for i in range(3,tf_diff.shape[0]-3):
    for j in range(4,tf_diff.shape[1]-4):
      
      #Start the analysis if a cell has experienced a tree cover decrease exceeding the threshold
      if(~np.isnan(T_diff[i,j]) and ~np.isnan(tf_diff[i,j]) and tf_diff[i,j] <= thres_tf):
	LUC_window_1 = tf_diff[i-2:i+3,j-2:j+3]
	LUC_window_3 = tf_diff[i-3:i+4,j-4:j+5]
	
	CtSt1 = float(np.count_nonzero(LUC_window_1 <= thres_tf)) / float(np.count_nonzero(~np.isnan(LUC_window_1)))
	CtSt3 = float(np.count_nonzero(LUC_window_3 <= thres_tf)) / float(np.count_nonzero(~np.isnan(LUC_window_3)))

	# If the conditions are met in the smallest box
	if(	0.35 <= CtSt1 <= 0.65 and np.count_nonzero(~np.isnan(LUC_window_1)) >= 8 and
		np.count_nonzero(LUC_window_1 <= thres_tf) >= 3 and
		np.count_nonzero(LUC_window_1 >= thres_tf) >= 3
	  ): 
	  
	  BB = 1
	    
	else:
	  LUC_window_2 = tf_diff[i-3:i+4,j-3:j+4]
	  CtSt2 = float(np.count_nonzero(LUC_window_2 <= thres_tf)) / float(np.count_nonzero(~np.isnan(LUC_window_2)))
	  
	  if(	0.35 <= CtSt2 <= 0.65 and np.count_nonzero(~np.isnan(LUC_window_2)) >= 8 and
		np.count_nonzero(LUC_window_2 <= thres_tf) >= 3 and
		np.count_nonzero(LUC_window_2 >= thres_tf) >= 3
	    ):
	    
	    BB = 2
	  
	  else:
	    LUC_window_3 = tf_diff[i-3:i+4,j-4:j+5]
	    CtSt3 = float(np.count_nonzero(LUC_window_3 <= thres_tf)) / float(np.count_nonzero(~np.isnan(LUC_window_3)))
	    
	    if(CtSt3 < 0.35 or CtSt3 > 0.65
	      ):
	      
	      best_CtSt = np.nanmin([abs(CtSt1-0.5),abs(CtSt2-0.5),abs(CtSt3-0.5)])
	      
	      if(	abs(CtSt1-0.5) == best_CtSt and np.count_nonzero(~np.isnan(LUC_window_1)) >= 8 and
			np.count_nonzero(LUC_window_1 <= thres_tf) >= 3 and
			np.count_nonzero(LUC_window_1 >= thres_tf) >= 3
		):
		
		BB = 1
		
	      elif(	abs(CtSt2-0.5) == best_CtSt and np.count_nonzero(~np.isnan(LUC_window_2)) >= 8 and
			np.count_nonzero(LUC_window_2 <= thres_tf) >= 3 and
			np.count_nonzero(LUC_window_2 >= thres_tf) >= 3
		):
		
		BB = 2
	      
	      elif(	abs(CtSt2-0.5) == best_CtSt and np.count_nonzero(~np.isnan(LUC_window_3)) >= 8 and
			np.count_nonzero(LUC_window_3 <= thres_tf) >= 3 and
			np.count_nonzero(LUC_window_3 >= thres_tf) >= 3
		):
		
		BB = 2
	    
	    elif(	0.35 <= CtSt3 <= 0.65 and np.count_nonzero(~np.isnan(LUC_window_3)) >= 8 and
			np.count_nonzero(LUC_window_3 <= thres_tf) >= 3 and
			np.count_nonzero(LUC_window_3 >= thres_tf) >= 3
		):
	      
	      BB = 3

      # Calculate the land use impact as:
      # The mean of the cells that have experience a tree cover decrease minus
      # The mean of the cells that have not experienced a tree cover decrease
      if BB == 3:
	T_window = T_diff[i-3:i+4,j-4:j+5]
	#print T_window.shape
	CtSt[i,j] = CtSt3
	LU_impact[i,j] = np.nanmean(T_window[LUC_window_3 <= thres_tf]) - np.nanmean(T_window[LUC_window_3 > thres_tf])
	count += 1
	#print BB,CtSt[i,j],LU_impact[i,j]
	
      elif BB == 2:
	T_window = T_diff[i-3:i+4,j-3:j+4]
	#print T_window.shape
	CtSt[i,j] = CtSt2
	LU_impact[i,j] = np.nanmean(T_window[LUC_window_2 <= thres_tf]) - np.nanmean(T_window[LUC_window_2 > thres_tf])
	count += 1
	#print BB,CtSt[i,j],LU_impact[i,j]
      
      elif BB == 1:
	T_window = T_diff[i-2:i+3,j-2:j+3]
	#print T_window.shape
	CtSt[i,j] = CtSt1
	LU_impact[i,j] = np.nanmean(T_window[LUC_window_1 <= thres_tf]) - np.nanmean(T_window[LUC_window_1 > thres_tf])
	count += 1
	#print BB,CtSt[i,j],LU_impact[i,j]
      
      elif BB == 0:
	LU_impact[i,j] = np.nan
	CtSt[i,j] = np.nan
		
      #Reinitialization
      BB = 0 
      
  print('COUNT:'+str(count))    
  count = 0
  return(LU_impact,CtSt)

modelname = 'MPI-ESM-MR'

if modelname == 'GFDL-CM3':
  dT_def = np.zeros((5,4,90,144))
  CtSt = np.zeros((5,4,90,144))
elif modelname == 'HadGEM2-ES':
  dT_def = np.zeros((4,4,145,192))
  CtSt = np.zeros((4,4,145,192))
elif modelname == 'IPSL-CM5A-LR':
  dT_def = np.zeros((6,4,96,96))
  CtSt = np.zeros((6,4,96,96))
elif modelname == 'IPSL-CM5A-MR':
  dT_def = np.zeros((3,4,143,144))
  CtSt = np.zeros((3,4,143,144))
elif(modelname == 'MPI-ESM-LR' or modelname == 'MPI-ESM-MR'):
  dT_def = np.zeros((3,4,96,192))
  CtSt = np.zeros((3,4,96,192))

for k in range(0,dT_def.shape[0]):
  for l in range(0,4):
    dT_def[k,l,:,:],CtSt[k,l,:,:] = calc_LU_impact(modelname,sim_list[k],seas_list[l],1861,1890,1975,2004)