"""SREX_irr_boxplots.py

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
execfile('calc_irr_impact_regr.py')
execfile('calc_irr_impact_thres.py')

#################################
#User-specified options
#################################
SREX_region = 'SAS'

datasource = "CRU_CESM"
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
figsave = True
figformat = 'png'
#################################

#Perform calculations
_,_, f_irr = calc_irr_diff(datasource,response,yr_start1,yr_end1,yr_start2,yr_end2)
T_out = extract_T_irr(datasource,temp_product,response,t_res,yr_start1,yr_end1,yr_start2,yr_end2)
print('Done importing changes in temperature and irrigated fraction')
dT_irr,_,_,_,_ = calc_irr_impact_regr(datasource,temp_product,response,t_res,thres_irr_PD,False,yr_start1,yr_end1,yr_start2,yr_end2)
#dT_irr, _,_ = calc_irr_impact_thres(datasource,temp_product,response,t_res,thres_irr_PD,0.5,yr_start1,yr_end1,yr_start2,yr_end2)

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

df_plot = []
dT_plot = []

# Arrange data in list format for plotting in subplots
df_plot.append(df_irr.ravel()[(~np.isnan(df_irr.ravel())) & (~np.isnan(dT.ravel()))])
#u = dT
#u[190:250,480:560] = np.nan
dT_plot.append(dT.ravel()[(~np.isnan(df_irr.ravel())) & (~np.isnan(dT.ravel()))])

df_plot.append(df_irr.ravel()[(~np.isnan(df_irr.ravel())) & (~np.isnan(dT_irr.ravel()))])
dT_plot.append(dT_irr.ravel()[(~np.isnan(df_irr.ravel())) & (~np.isnan(dT_irr.ravel()))])

if SREX_region == 'CNA':
  df_plot.append(df_irr[237:280,150:190].ravel()[(~np.isnan(df_irr[237:280,150:190].ravel())) & (~np.isnan(dT[237:280,150:190].ravel()))])
  dT_plot.append(dT[237:280,150:190].ravel()[(~np.isnan(df_irr[237:280,150:190].ravel())) & (~np.isnan(dT[237:280,150:190].ravel()))])
elif SREX_region == 'EAS':
  df_plot.append(df_irr[220:280,560:650].ravel()[(~np.isnan(df_irr[220:280,560:650].ravel())) & (~np.isnan(dT[220:280,560:650].ravel()))])
  dT_plot.append(dT[220:280,560:650].ravel()[(~np.isnan(df_irr[220:280,560:650].ravel())) & (~np.isnan(dT[220:280,560:650].ravel()))])
elif SREX_region == 'MED':
  df_plot.append(df_irr[240:270,340:440].ravel()[(~np.isnan(df_irr[240:270,340:440].ravel())) & (~np.isnan(dT[240:270,340:440].ravel()))])
  dT_plot.append(dT[240:270,340:440].ravel()[(~np.isnan(df_irr[240:270,340:440].ravel())) & (~np.isnan(dT[240:270,340:440].ravel()))])
elif SREX_region == 'SAS':
  df_SAS1 = df_irr[190:250,480:550].ravel()[(~np.isnan(df_irr[190:250,480:550].ravel())) & (~np.isnan(dT[190:250,480:550].ravel()))]
  df_SAS2 = df_irr[220:250,550:560].ravel()[(~np.isnan(df_irr[220:250,550:560].ravel())) & (~np.isnan(dT[220:250,550:560].ravel()))]
  df_plot.append(np.append(df_SAS1,df_SAS2))
  dT_SAS1 = dT[190:250,480:550].ravel()[(~np.isnan(df_irr[190:250,480:550].ravel())) & (~np.isnan(dT[190:250,480:550].ravel()))]
  dT_SAS2 = dT[220:250,550:560].ravel()[(~np.isnan(df_irr[220:250,550:560].ravel())) & (~np.isnan(dT[220:250,550:560].ravel()))]
  dT_plot.append(np.append(dT_SAS1,dT_SAS2))
  del df_SAS1,df_SAS2,dT_SAS1,dT_SAS2
elif SREX_region == 'SEA':
  df_plot.append(df_irr[160:220,550:670].ravel()[(~np.isnan(df_irr[160:220,550:670].ravel())) & (~np.isnan(dT[160:220,550:670].ravel()))])
  dT_plot.append(dT[160:220,550:670].ravel()[(~np.isnan(df_irr[160:220,550:670].ravel())) & (~np.isnan(dT[160:220,550:670].ravel()))])
elif SREX_region == 'WAS':
  df_plot.append(df_irr[210:280,440:480].ravel()[(~np.isnan(df_irr[210:280,440:480].ravel())) & (~np.isnan(dT[210:280,440:480].ravel()))])
  dT_plot.append(dT[210:280,440:480].ravel()[(~np.isnan(df_irr[210:280,440:480].ravel())) & (~np.isnan(dT[210:280,440:480].ravel()))])
elif SREX_region == 'WNA':
  df_plot.append(df_irr[237:300,100:150].ravel()[(~np.isnan(df_irr[237:300,100:150].ravel())) & (~np.isnan(dT[237:300,100:150].ravel()))])
  dT_plot.append(dT[237:300,100:150].ravel()[(~np.isnan(df_irr[237:300,100:150].ravel())) & (~np.isnan(dT[237:300,100:150].ravel()))])
  
if SREX_region == 'CNA':
  df_plot.append(df_irr[237:280,150:190].ravel()[(~np.isnan(df_irr[237:280,150:190].ravel())) & (~np.isnan(dT_irr[237:280,150:190].ravel()))])
  dT_plot.append(dT_irr[237:280,150:190].ravel()[(~np.isnan(df_irr[237:280,150:190].ravel())) & (~np.isnan(dT_irr[237:280,150:190].ravel()))])
elif SREX_region == 'EAS':
  df_plot.append(df_irr[220:280,560:650].ravel()[(~np.isnan(df_irr[220:280,560:650].ravel())) & (~np.isnan(dT_irr[220:280,560:650].ravel()))])
  dT_plot.append(dT_irr[220:280,560:650].ravel()[(~np.isnan(df_irr[220:280,560:650].ravel())) & (~np.isnan(dT_irr[220:280,560:650].ravel()))])
elif SREX_region == 'MED':
  df_plot.append(df_irr[240:270,340:440].ravel()[(~np.isnan(df_irr[240:270,340:440].ravel())) & (~np.isnan(dT_irr[240:270,340:440].ravel()))])
  dT_plot.append(dT_irr[240:270,340:440].ravel()[(~np.isnan(df_irr[240:270,340:440].ravel())) & (~np.isnan(dT_irr[240:270,340:440].ravel()))])
elif SREX_region == 'SAS':
  df_irr_SAS1 = df_irr[190:250,480:550].ravel()[(~np.isnan(df_irr[190:250,480:550].ravel())) & (~np.isnan(dT_irr[190:250,480:550].ravel()))]
  df_irr_SAS2 = df_irr[220:250,550:560].ravel()[(~np.isnan(df_irr[220:250,550:560].ravel())) & (~np.isnan(dT_irr[220:250,550:560].ravel()))]
  df_plot.append(np.append(df_irr_SAS1,df_irr_SAS2))
  dT_irr_SAS1 = dT_irr[190:250,480:550].ravel()[(~np.isnan(df_irr[190:250,480:550].ravel())) & (~np.isnan(dT_irr[190:250,480:550].ravel()))]
  dT_irr_SAS2 = dT_irr[220:250,550:560].ravel()[(~np.isnan(df_irr[220:250,550:560].ravel())) & (~np.isnan(dT_irr[220:250,550:560].ravel()))]
  dT_plot.append(np.append(dT_irr_SAS1,dT_irr_SAS2))
  del df_irr_SAS1,df_irr_SAS2,dT_irr_SAS1,dT_irr_SAS2
elif SREX_region == 'SEA':
  df_plot.append(df_irr[160:220,550:670].ravel()[(~np.isnan(df_irr[160:220,550:670].ravel())) & (~np.isnan(dT_irr[160:220,550:670].ravel()))])
  dT_plot.append(dT_irr[160:220,550:670].ravel()[(~np.isnan(df_irr[160:220,550:670].ravel())) & (~np.isnan(dT_irr[160:220,550:670].ravel()))])
elif SREX_region == 'WAS':
  df_plot.append(df_irr[210:280,440:480].ravel()[(~np.isnan(df_irr[210:280,440:480].ravel())) & (~np.isnan(dT_irr[210:280,440:480].ravel()))])
  dT_plot.append(dT_irr[210:280,440:480].ravel()[(~np.isnan(df_irr[210:280,440:480].ravel())) & (~np.isnan(dT_irr[210:280,440:480].ravel()))])
elif SREX_region == 'WNA':
  df_plot.append(df_irr[237:300,100:150].ravel()[(~np.isnan(df_irr[237:300,100:150].ravel())) & (~np.isnan(dT_irr[237:300,100:150].ravel()))])
  dT_plot.append(dT_irr[237:300,100:150].ravel()[(~np.isnan(df_irr[237:300,100:150].ravel())) & (~np.isnan(dT_irr[237:300,100:150].ravel()))])

f, axs = plt.subplots(2,2,sharex='col',sharey='row',figsize=(12,10))
SREX_list = ['Global land','Global irrigated land','%s'%SREX_region,'%s irrigated'%SREX_region]

b_low = round(round(min([x.min() for x in df_plot])/0.05)*0.05, -int(np.floor(np.log10(0.05))))
b_high = round(round(max([x.max() for x in df_plot])/0.05)*0.05, -int(np.ceil(np.log10(0.05))))
bins = np.round(np.linspace(b_low,b_high,(abs(b_low)+abs(b_high))/0.05+1),2)
bins_plot = bins[:-1]+0.025

#dT_plot = [x/8. for x in dT_plot]

for ii in range(0,2):
  for jj in range(0,2):
    data_dT = dT_plot[2*ii+jj]
    data_df = df_plot[2*ii+jj]
    
    means = []
    pctls = []
    counts = []
    for x in range(0,len(bins)-1):
      if data_dT[(data_df >= bins[x]) & (data_df < bins[x+1])].shape[0] >= 10:
        counts = data_dT[(data_df >= bins[x]) & (data_df < bins[x+1])].shape[0]
        #means.append(data_dT[(data_df >= bins[x]) & (data_df < bins[x+1])].mean())
        means.append(np.median(data_dT[(data_df >= bins[x]) & (data_df < bins[x+1])]))
        pctls.append(np.percentile(data_dT[(data_df >= bins[x]) & (data_df < bins[x+1])],[25,75]))
      else:
        counts = data_dT[(data_df >= bins[x]) & (data_df < bins[x+1])].shape[0]
        means.append(np.nan)
        pctls.append(np.array([np.nan,np.nan]))
    
    pctls_l = [x[0] for x in pctls]
    pctls_h = [x[1] for x in pctls]
    
    axs[ii,jj].plot(bins_plot,means,'ks',markersize=4)
    ax2 = axs[ii,jj].twinx()
    ax2.hist(data_df,bins,color='white',alpha=0.5)
    
    if ii == 0:
      ax2.set_ylim([0,8000])
      if jj == 1:
        ax2.set_yticks([0,2000,4000])
      else:
        ax2.set_yticks([])
        #axs[ii,jj].axhline(y=2,c='k',linewidth=3)
      #axs[ii,jj].axhline(y=-2,c='k',linewidth=3)
      #axs[ii,jj].axvline(x=-0.2,c='k',linewidth=3)
      #axs[ii,jj].axvline(x=0.85,c='k',linewidth=3)
    elif ii == 1:
      ax2.set_ylim([0,800])
      if jj == 1:
        ax2.set_yticks([0,200,400])
      else:
        ax2.set_yticks([])
    
    for u in range(0,len(pctls_l)):
      axs[ii,jj].plot((bins_plot[u],bins_plot[u]), (pctls_l[u],pctls_h[u]), 'k-',linewidth=1.5)
      
    #axs[ii,jj].set_xlim([b_low+0.05,b_high-0.05])
    axs[ii,jj].set_xlim([0,0.9])
    axs[ii,jj].set_xticks(np.arange(0,1.2,0.3))
    axs[ii,jj].set_yticks(np.arange(-2,3,1))
    axs[ii,jj].annotate('n=%i'%(data_df.shape[0]),xy=(0.88,1.9),ha='right',va='top')
    axs[ii,jj].grid(True)
    axs[ii,jj].fill_between([-2,2],0,3,color='salmon',alpha=0.5)
    axs[ii,jj].fill_between([-2,2],-3,0,color='lightblue',alpha=0.7)
    axs[ii,jj].set_title(SREX_list[2*ii+jj])
    axs[ii,jj].set_ylim([-2,2])
    
#axs[0,3].set_ylim([-0.25,0.25])
#axs[1,3].set_ylim([-0.25,0.25])
f.text(0.47, 0.02, '$f_{irr}$ [-]', ha='center', va='center',fontsize=16)
f.text(0.001, 0.5, '$\Delta T$ [$^\circ$C]', ha='center', va='center', rotation='vertical',fontsize=16)
f.text(0.49, 0.5, '$\Delta T_{irr}$ [$^\circ$C]', ha='center', va='center', rotation='vertical',fontsize=16)
f.text(0.98, 0.5, 'Count [-]', ha='center', va='center', rotation=270,fontsize=16)
f.text(0.47, 0.96, '$t_{ref}$: %i-%i, $t_{PD}$: %i-%i'%(yr_start1,yr_end1,yr_start2,yr_end2),ha='center',va='center',fontsize=14)
plt.subplots_adjust(left=0.03, right=0.91, top=0.9, bottom=0.1)
#plt.show(block=False)
  
#Save results  
if figsave == True:
  figpath = '/net/exo/landclim/wthiery/observational_analysis/Figures/%s/%s/binnedboxplots/binnedboxplots_hist.%s.%s.%s.%s.all.%1.2f.%i-%i.%i-%i.TEST.%s'%(datasource,temp_product,datasource,temp_product,SREX_region,response,thres_irr_PD,yr_start1,yr_end1,yr_start2,yr_end2,figformat)
  
  if os.path.exists(figpath):
    print('Figure already exists at path: ' + figpath)
  else:
    print('Saving figure as:' + figpath)
    plt.savefig(figpath,bbox_inches='tight')
