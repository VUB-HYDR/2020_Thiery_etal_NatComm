"""SREX_scatter_binned.py

author: Auke Visser
date: 12.12.2016

This script makes scatter plots of df_irr and dT based to figure
out the type of response between df_irr and dT.

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
from scipy import stats
import matplotlib.lines as mlines

execfile('extract_T_irr.py')
execfile('calc_irr_diff.py')

#################################
#User-specified options
#################################
SREX_region = 'SAS'

datasource = "CRU"
temp_product = "tmx_max"
response = "PD"
t_res = "seasonal"

seas_ind = 2
month_ind = 5

yr_start1 = 1901
yr_end1 = 1930
yr_start2 = 1981
yr_end2 = 2010

thres_irr_PD = 0.02

p_value = 0.01
figsave = False
figformat = 'png'
#################################

#Perform calculations
_,_, f_irr = calc_irr_diff(datasource,response,yr_start1,yr_end1,yr_start2,yr_end2)
T_out = extract_T_irr(datasource,temp_product,response,t_res,yr_start1,yr_end1,yr_start2,yr_end2)
print('Done importing changes in temperature and irrigated fraction')

#df_irr is defined differently based on the user-selected response
if response == 'PC/PD':
  df_irr = f_irr[1,:,:] - f_irr[0,:,:]
  response = 'PC-PD'
elif response == 'PD':
  df_irr = f_irr[1,:,:]
df_irr[f_irr[1,:,:] < thres_irr_PD] = np.nan

if(temp_product in ['tmp','dtr','tmn','tmx'] and t_res == 'seasonal'):
  dT = T_out[1,seas_ind,:,:] - T_out[0,seas_ind,:,:]
elif(temp_product in ['tmp','dtr','tmn','tmx'] and t_res == 'monthly'):
  dT = T_out[1,month_ind,:,:] - T_out[0,month_ind,:,:]
elif(temp_product in ['tmx_max','tmp_max','tmn_max']):
  dT = T_out[1,:,:] - T_out[0,:,:]

#df_plot and dT_plot are defined as empty lists, to which the values for every region are appended
#The final list has 8 entries of 1D-arrays. 
df_plot = []
dT_plot = []

#Select pixels for all regions and append them to the list, ensuring that neither the dT_irr-value nor 
#the df_irr value is NaN
df_plot.append(df_irr.ravel()[(~np.isnan(df_irr.ravel())) & (~np.isnan(dT.ravel()))])
#u = dT
#u[190:250,480:560] = np.nan
dT_plot.append(dT.ravel()[(~np.isnan(df_irr.ravel())) & (~np.isnan(dT.ravel()))])
df_plot.append(df_irr[237:280,150:190].ravel()[(~np.isnan(df_irr[237:280,150:190].ravel())) & (~np.isnan(dT[237:280,150:190].ravel()))])
dT_plot.append(dT[237:280,150:190].ravel()[(~np.isnan(df_irr[237:280,150:190].ravel())) & (~np.isnan(dT[237:280,150:190].ravel()))])
df_plot.append(df_irr[220:280,560:650].ravel()[(~np.isnan(df_irr[220:280,560:650].ravel())) & (~np.isnan(dT[220:280,560:650].ravel()))])
dT_plot.append(dT[220:280,560:650].ravel()[(~np.isnan(df_irr[220:280,560:650].ravel())) & (~np.isnan(dT[220:280,560:650].ravel()))])
df_plot.append(df_irr[240:270,340:440].ravel()[(~np.isnan(df_irr[240:270,340:440].ravel())) & (~np.isnan(dT[240:270,340:440].ravel()))])
dT_plot.append(dT[240:270,340:440].ravel()[(~np.isnan(df_irr[240:270,340:440].ravel())) & (~np.isnan(dT[240:270,340:440].ravel()))])
df_SAS1 = df_irr[190:250,480:550].ravel()[(~np.isnan(df_irr[190:250,480:550].ravel())) & (~np.isnan(dT[190:250,480:550].ravel()))]
df_SAS2 = df_irr[220:250,550:560].ravel()[(~np.isnan(df_irr[220:250,550:560].ravel())) & (~np.isnan(dT[220:250,550:560].ravel()))]
df_plot.append(np.append(df_SAS1,df_SAS2))
dT_SAS1 = dT[190:250,480:550].ravel()[(~np.isnan(df_irr[190:250,480:550].ravel())) & (~np.isnan(dT[190:250,480:550].ravel()))]
dT_SAS2 = dT[220:250,550:560].ravel()[(~np.isnan(df_irr[220:250,550:560].ravel())) & (~np.isnan(dT[220:250,550:560].ravel()))]
dT_plot.append(np.append(dT_SAS1,dT_SAS2))
del df_SAS1,df_SAS2,dT_SAS1,dT_SAS2
df_plot.append(df_irr[160:220,550:670].ravel()[(~np.isnan(df_irr[160:220,550:670].ravel())) & (~np.isnan(dT[160:220,550:670].ravel()))])
dT_plot.append(dT[160:220,550:670].ravel()[(~np.isnan(df_irr[160:220,550:670].ravel())) & (~np.isnan(dT[160:220,550:670].ravel()))])
df_plot.append(df_irr[210:280,440:480].ravel()[(~np.isnan(df_irr[210:280,440:480].ravel())) & (~np.isnan(dT[210:280,440:480].ravel()))])
dT_plot.append(dT[210:280,440:480].ravel()[(~np.isnan(df_irr[210:280,440:480].ravel())) & (~np.isnan(dT[210:280,440:480].ravel()))])
df_plot.append(df_irr[237:300,100:150].ravel()[(~np.isnan(df_irr[237:300,100:150].ravel())) & (~np.isnan(dT[237:300,100:150].ravel()))])
dT_plot.append(dT[237:300,100:150].ravel()[(~np.isnan(df_irr[237:300,100:150].ravel())) & (~np.isnan(dT[237:300,100:150].ravel()))])

#Define a figure with 4x2 subplots
f, axs = plt.subplots(2,4,sharex='col',sharey='row',figsize=(20,9))
SREX_list = ['Global land','Central North America','East Asia','Mediterranean','South Asia','Southeast Asia','West Asia','West North America']

#Define the lower and upper plotting bound based on the values in df.
b_low = round(round(min([x.min() for x in df_plot])/0.05)*0.05, -int(np.floor(np.log10(0.05))))
b_high = round(round(max([x.max() for x in df_plot])/0.05)*0.05, -int(np.ceil(np.log10(0.05))))
bins = np.round(np.linspace(b_low,b_high,(abs(b_low)+abs(b_high))/0.05+1),2)
bins_plot = bins[:-1]+0.025

#dT_plot = [x/8. for x in dT_plot]

#Fill the subplots in a for loop
for ii in range(0,2):
  for jj in range(0,4):
    #Select the df and dT values of interest
    data_dT = dT_plot[4*ii+jj]
    data_df = df_plot[4*ii+jj]
    
    #Define lists that contain the binned medians, percentiles and counts for the subplot
    medians = []
    pctls = []
    counts = []
    
    #Add the median, Q25 and Q75 if there is enough data, otherwise add
    #NaNs to ensure the variable is not plotted
    for x in range(0,len(bins)-1):
      if data_dT[(data_df >= bins[x]) & (data_df < bins[x+1])].shape[0] >= 10:
        counts = data_dT[(data_df >= bins[x]) & (data_df < bins[x+1])].shape[0]
        medians.append(np.median(data_dT[(data_df >= bins[x]) & (data_df < bins[x+1])]))
        pctls.append(np.percentile(data_dT[(data_df >= bins[x]) & (data_df < bins[x+1])],[25,75]))
      else:
        counts = data_dT[(data_df >= bins[x]) & (data_df < bins[x+1])].shape[0]
        medians.append(np.nan)
        pctls.append(np.array([np.nan,np.nan]))
    
    pctls_l = [x[0] for x in pctls]
    pctls_h = [x[1] for x in pctls]
    
    #Some extra plotting options
    axs[ii,jj].plot(bins_plot,medians,'ks',markersize=4)
    ax2 = axs[ii,jj].twinx()
    ax2.hist(data_df,bins,color='white',alpha=0.4)
    ax2.set_yscale('log')
    
    if 4*ii+jj == 0:
      ax2.set_ylim([10,1e9])
      ax2.set_yticks([10,1e3,1e5])
    elif jj == 3:
      ax2.set_ylim([1,1e8])
      ax2.set_yticks([1,1e2,1e4])
    else:
      ax2.set_ylim([1,1e8])
      ax2.set_yticks([])
    
    for u in range(0,len(pctls_l)):
      axs[ii,jj].plot((bins_plot[u],bins_plot[u]), (pctls_l[u],pctls_h[u]), 'k-',linewidth=1.5)
      
    axs[ii,jj].set_xlim([0,0.9])
    axs[ii,jj].set_xticks(np.arange(0,1.2,0.3))
    axs[ii,jj].set_yticks(np.arange(-2,3,1))
    axs[ii,jj].annotate('n=%i'%(data_df.shape[0]),xy=(0.88,1.9),ha='right',va='top')
    axs[ii,jj].grid(True)
    axs[ii,jj].fill_between([-2,2],0,3,color='salmon',alpha=0.5)
    axs[ii,jj].fill_between([-2,2],-3,0,color='lightblue',alpha=0.7)
    axs[ii,jj].set_title(SREX_list[4*ii+jj])
    axs[ii,jj].set_ylim([-2,2])

#Plot the axis labels
f.text(0.47, 0.02, '$f_{irr}$ [-]', ha='center', va='center',fontsize=16)
f.text(0.01, 0.5, '$\Delta T$ [$^\circ$C]', ha='center', va='center', rotation='vertical',fontsize=16)
f.text(0.95, 0.5, 'Count [-]', ha='center', va='center', rotation=270,fontsize=16)
f.text(0.47, 0.96, '$t_{ref}$: %i-%i, $t_{PD}$: %i-%i'%(yr_start1,yr_end1,yr_start2,yr_end2),ha='center',va='center',fontsize=14)
plt.subplots_adjust(left=0.03, right=0.91, top=0.9, bottom=0.1)
plt.show(block=False)
  
#Save results  
if figsave == True:
  figpath = '/net/exo/landclim/wthiery/observational_analysis/Figures/%s/%s/binnedboxplots/binnedboxplots_hist.%s.%s.%s.all.%1.2f.%i-%i.%i-%i.%s'%(datasource,temp_product,datasource,temp_product,response,thres_irr_PD,yr_start1,yr_end1,yr_start2,yr_end2,figformat)
  
  if os.path.exists(figpath):
    print('Figure already exists at path: ' + figpath)
  else:
    print('Saving figure as:' + figpath)
    plt.savefig(figpath,bbox_inches='tight')
