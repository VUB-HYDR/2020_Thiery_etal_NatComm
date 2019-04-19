"""dT_irr_SREX.py

author: Auke Visser
date: 27.10.2016

This script calls a routine to calculate the irrigation impact on 
temperature with either the threshold- or the regression-based
approach. The user has the option to print or visualize output.

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
import scipy.stats as stats
import matplotlib.patches as mpatches

execfile('calc_irr_diff.py')

###############################
# User-specified options
###############################

#Choose region, season and time period
SREX_region = 'GIL' #'GIL' (global irrigated land), 'CNA', 'EAS', 'MED', 'SAS', 'SEA', 'WAS', 'WNA'
EOBS_region = 'MED' #'EUR', 'MED' (Only works for E-OBS)
yr_start1 = 1901
yr_end1 = 1930
yr_start2 = 1981
yr_end2 = 2010

#Select method, data product and response

datasource = 'CRU' #'CRU', 'E-OBS'
response = 'PC/PD'

binwidth=0.025
figshow = True
figsave = True
figformat = 'png' #'pdf', 'png'
###############################
# End of user-specified options
###############################

#perform calculations
lat,lon,df_irr = calc_irr_diff(datasource,response,yr_start1,yr_end1,yr_start2,yr_end2)
HIDfile = nc.Dataset('/net/exo/landclim/wthiery/observational_analysis/Data/HID_regridded/hid_v1.0.f_irr_hyde_final_ir.0.5deg.10y.1900-2005.nc','r')
f_irr = HIDfile.variables['aei_hyde_final_ir'][:][12,:,:].data
f_irr[f_irr == -9999] = np.nan

if datasource == 'CRU':
  if SREX_region == 'GIL':
    df_array = df_irr.ravel()[~np.isnan(df_irr.ravel())]
    f_array = f_irr.ravel()[~np.isnan(f_irr.ravel())]
  elif SREX_region == 'CNA':
    df_array = df_irr[237:280,150:190].ravel()[~np.isnan(df_irr[237:280,150:190].ravel())]
    f_array = f_irr[237:280,150:190].ravel()[~np.isnan(f_irr[237:280,150:190].ravel())]
  elif SREX_region == 'EAS':
    df_array = df_irr[220:280,560:650].ravel()[~np.isnan(df_irr[220:280,560:650].ravel())]
    f_array = f_irr[220:280,560:650].ravel()[~np.isnan(f_irr[220:280,560:650].ravel())]
  elif SREX_region == 'MED':
    df_array = df_irr[240:270,340:440].ravel()[~np.isnan(df_irr[240:270,340:440].ravel())]
    f_array = f_irr[240:270,340:440].ravel()[~np.isnan(f_irr[240:270,340:440].ravel())]
  elif SREX_region == 'SAS':
    SAS1 = df_irr[190:240,480:550].ravel()[~np.isnan(df_irr[190:240,480:550].ravel())]
    SAS2 = df_irr[220:240,550:560].ravel()[~np.isnan(df_irr[220:240,550:560].ravel())]
    df_array = np.append(SAS1,SAS2)
    
    S1 = f_irr[190:240,480:550].ravel()[~np.isnan(f_irr[190:240,480:550].ravel())]
    S2 = f_irr[220:240,550:560].ravel()[~np.isnan(f_irr[220:240,550:560].ravel())]
    f_array = np.append(S1,S2)
  elif SREX_region == 'SEA':
    df_array = df_irr[160:220,550:670].ravel()[~np.isnan(df_irr[160:220,550:670].ravel())]
    f_array = f_irr[160:220,550:670].ravel()[~np.isnan(f_irr[160:220,550:670].ravel())]
  elif SREX_region == 'WAS':
    df_array = df_irr[210:280,440:480].ravel()[~np.isnan(df_irr[210:280,440:480].ravel())]
    f_array = f_irr[210:280,440:480].ravel()[~np.isnan(f_irr[210:280,440:480].ravel())]
  elif SREX_region == 'WNA':
    df_array = df_irr[237:300,100:150].ravel()[~np.isnan(df_irr[237:300,100:150].ravel())]
    f_array = f_irr[237:300,100:150].ravel()[~np.isnan(f_irr[237:300,100:150].ravel())]
  
  print np.percentile(df_array[f_array >= 0.05],[25,50,75])
  print np.nanmean(df_array[f_array >= 0.05])
  
  if figshow == True:
    bins = np.linspace(0,1,1/binwidth+1)
    plt.figure(figsize=(12,8))
    plt.hist(df_array,bins=bins,color='green',alpha=0.5,label='%s'%SREX_region)
    plt.grid(True)
    plt.xlabel('$\Delta f_{irr}$ bins',size=15)
    plt.ylabel('Count [-]',size=15)
    plt.yscale('log', nonposy='clip')
    plt.legend()
    plt.title('region: %s, count: %i, ref: %i-%i'%(SREX_region,df_array.shape[0],yr_start1,yr_end1),size=18,weight='bold')
    plt.show(block=False)
    
  if figsave == True:
    figpath = '/net/exo/landclim/wthiery/observational_analysis/Figures/HID/df_irr.hist.%1.3f.%s.CRU.%s-%s.%s-%s.%s'%(binwidth,SREX_region,yr_start1,yr_end1,yr_start2,yr_end2,figformat)
    if os.path.exists(figpath):
      print('Figure already exists at path: ' + figpath)
    else:  
      print('Saving figure as: ' + figpath)
      plt.savefig(figpath,bbox_inches='tight')
    
    
