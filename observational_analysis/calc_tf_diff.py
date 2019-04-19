"""calc_tf_diff.py

author: Auke Visser
date: 25.10.2016

This script calculates the change in tree fraction between two 30-year time periods from the historical 
LUC dataset from Meiyappan and Jain (2012). 


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

def calc_tf_diff(datasource,response,yr_start1,yr_end1,yr_start2,yr_end2):
  
  if datasource == "CRU":
    PC_tf = np.zeros((30,360,720))
    PD_tf = np.zeros((30,360,720))
    
    #Access variables from the Meiyappan and Jain (2012) Historical LUC data set
    for i in range(0,30):
      PC_filename = '/net/exo/landclim/wthiery/observational_analysis/Data/LUC/HYDE_AREAVEG/land-cover_hyde_landcover_yr%i.nc'%(i+yr_start1)
      PD_filename = '/net/exo/landclim/wthiery/observational_analysis/Data/LUC/HYDE_AREAVEG/land-cover_hyde_landcover_yr%i.nc'%(i+yr_start2)
      
      PC_tffile = nc.Dataset(PC_filename,'r')
      PD_tffile = nc.Dataset(PD_filename,'r')
      
      PC_tf[i,:,:] = PC_tffile.variables['TrpEBF'][:][0]+PC_tffile.variables['TrpDBF'][:][0]+PC_tffile.variables['TmpEBF'][:][0]+PC_tffile.variables['TmpENF'][:][0]+PC_tffile.variables['TmpDBF'][:][0]+PC_tffile.variables['BorENF'][:][0]+PC_tffile.variables['BorDNF'][:][0]
      PD_tf[i,:,:] = PD_tffile.variables['TrpEBF'][:][0]+PD_tffile.variables['TrpDBF'][:][0]+PD_tffile.variables['TmpEBF'][:][0]+PD_tffile.variables['TmpENF'][:][0]+PD_tffile.variables['TmpDBF'][:][0]+PD_tffile.variables['BorENF'][:][0]+PD_tffile.variables['BorDNF'][:][0]
      
      PC_tf[i,:,:][PC_tf[i,:,:] == -9900.0] = np.nan
      PD_tf[i,:,:][PD_tf[i,:,:] == -9900.0] = np.nan

    tf_diff = (np.nanmean(PD_tf,axis=0) - np.nanmean(PC_tf,axis=0)) * 0.01
    
    #Shift center of dataset to (0,0) instead of (0,180)
    tf_diff_shift = np.zeros((360,720))
    tf_diff_shift[:,0:360] = tf_diff[:,360:720]
    tf_diff_shift[:,360:720] = tf_diff[:,0:360]
      
  return tf_diff_shift
    
#tf_diff = calc_tf_diff("CRU","PC/PD",1901,1930,1981,2010)
