"""SREX_scatter.py

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
response = "PC/PD"
t_res = "seasonal"

seas_ind = 2
month_ind = 5

yr_start1 = 1901
yr_end1 = 1930
yr_start2 = 1981
yr_end2 = 2010

plot_opt = 'all' #'region', 'all'
thres_irr_PD = 0.02

p_value = 0.01
figsave = True
figformat = 'pdf'
#################################

#Perform calculations
_,_, f_irr = calc_irr_diff(datasource,response,yr_start1,yr_end1,yr_start2,yr_end2)
T_out = extract_T_irr(datasource,temp_product,response,t_res,yr_start1,yr_end1,yr_start2,yr_end2)
print('Done importing changes in temperature and irrigated fraction')

df_irr = f_irr[1,:,:] - f_irr[0,:,:]
df_irr[f_irr[1,:,:] < thres_irr_PD] = np.nan

if(temp_product in ['tmp','dtr','tmn','tmx'] and t_res == 'seasonal'):
  dT = T_out[1,seas_ind,:,:] - T_out[0,seas_ind,:,:]
elif(temp_product in ['tmp','dtr','tmn','tmx'] and t_res == 'monthly'):
  dT = T_out[1,month_ind,:,:] - T_out[0,month_ind,:,:]
elif(temp_product == 'tmx_max'):
  dT = T_out[1,:,:] - T_out[0,:,:]

if plot_opt == 'region':
  if SREX_region == 'GIL':
    df_plot = df_irr.ravel()[(~np.isnan(df_irr.ravel())) & (~np.isnan(dT.ravel()))]
    dT_plot = dT.ravel()[(~np.isnan(df_irr.ravel())) & (~np.isnan(dT.ravel()))]
  elif SREX_region == 'CNA':
    df_plot = df_irr[237:280,150:190].ravel()[(~np.isnan(df_irr[237:280,150:190].ravel())) & (~np.isnan(dT[237:280,150:190].ravel()))]
    dT_plot = dT[237:280,150:190].ravel()[(~np.isnan(df_irr[237:280,150:190].ravel())) & (~np.isnan(dT[237:280,150:190].ravel()))]
  elif SREX_region == 'EAS':
    df_plot = df_irr[220:280,560:650].ravel()[(~np.isnan(df_irr[220:280,560:650].ravel())) & (~np.isnan(dT[220:280,560:650].ravel()))]
    dT_plot = dT[220:280,560:650].ravel()[(~np.isnan(df_irr[220:280,560:650].ravel())) & (~np.isnan(dT[220:280,560:650].ravel()))]
  elif SREX_region == 'MED':
    df_plot = df_irr[240:270,340:440].ravel()[(~np.isnan(df_irr[240:270,340:440].ravel())) & (~np.isnan(dT[240:270,340:440].ravel()))]
    dT_plot = dT[240:270,340:440].ravel()[(~np.isnan(df_irr[240:270,340:440].ravel())) & (~np.isnan(dT[240:270,340:440].ravel()))]
  elif SREX_region == 'SAS':
    df_plot1 = df_irr[190:250,480:550].ravel()[(~np.isnan(df_irr[190:250,480:550].ravel())) & (~np.isnan(dT[190:250,480:550].ravel()))]
    df_plot2 = df_irr[220:250,550:560].ravel()[(~np.isnan(df_irr[220:250,550:560].ravel())) & (~np.isnan(dT[220:250,550:560].ravel()))]
    df_plot = np.append(df_plot1,df_plot2)
    
    dT_plot1 = dT[190:250,480:550].ravel()[(~np.isnan(df_irr[190:250,480:550].ravel())) & (~np.isnan(dT[190:250,480:550].ravel()))]
    dT_plot2 = dT[220:250,550:560].ravel()[(~np.isnan(df_irr[220:250,550:560].ravel())) & (~np.isnan(dT[220:250,550:560].ravel()))]
    dT_plot = np.append(dT_plot1,dT_plot2)
  elif SREX_region == 'SEA':
    df_plot = df_irr[160:220,550:670].ravel()[(~np.isnan(df_irr[160:220,550:670].ravel())) & (~np.isnan(dT[160:220,550:670].ravel()))]
    dT_plot = dT[160:220,550:670].ravel()[(~np.isnan(df_irr[160:220,550:670].ravel())) & (~np.isnan(dT[160:220,550:670].ravel()))]
  elif SREX_region == 'WAS':
    df_plot = df_irr[210:280,440:480].ravel()[(~np.isnan(df_irr[210:280,440:480].ravel())) & (~np.isnan(dT[210:280,440:480].ravel()))]
    dT_plot = dT[210:280,440:480].ravel()[(~np.isnan(df_irr[210:280,440:480].ravel())) & (~np.isnan(dT[210:280,440:480].ravel()))]
  elif SREX_region == 'WNA':
    df_plot = df_irr[237:300,100:150].ravel()[(~np.isnan(df_irr[237:300,100:150].ravel())) & (~np.isnan(dT[237:300,100:150].ravel()))]
    dT_plot = dT[237:300,100:150].ravel()[(~np.isnan(df_irr[237:300,100:150].ravel())) & (~np.isnan(dT[237:300,100:150].ravel()))]
    
  #Calculate statistics
  #Experiment with non-parametric Theil-Sen slope estimate
  slope,_,lb,ub = stats.theilslopes(df_plot,dT_plot,alpha=1-p_value)
  S = stats.linregress(df_plot,dT_plot)
  rmse = np.linalg.norm(np.poly1d(np.polyfit(df_plot,dT_plot,1))(df_plot[df_plot>0]) - dT_plot[df_plot>0]) / np.sqrt(df_plot[df_plot>0].shape[0])

  plt.figure(figsize=(12,8))
  plt.scatter(df_plot,dT_plot,color='black',s=2,label='data points')
  if ub < 0:
    plt.plot(np.unique(df_plot), np.poly1d(np.polyfit(df_plot,dT_plot,1))(np.unique(df_plot)),'g-',linewidth=2,label='1st order fit')
    plt.plot(np.unique(df_plot), np.poly1d(np.polyfit(df_plot,dT_plot,2))(np.unique(df_plot)),'g--',linewidth=2,label='2nd order fit')
  else:
    plt.plot(np.unique(df_plot), np.poly1d(np.polyfit(df_plot,dT_plot,1))(np.unique(df_plot)),'b-',label='1st order fit')
    plt.plot(np.unique(df_plot), np.poly1d(np.polyfit(df_plot,dT_plot,2))(np.unique(df_plot)),'b--',label='2nd order fit')

  plt.legend(loc='upper right')
  plt.xlabel('$\Delta f_{irr}$ [-]')
  plt.ylabel('$\Delta T$ [$^\circ$C]')
  plt.text(0.5,-1.7,'$r^2 = %1.3f$'%(S[2]**2),ha='left',va='top')
  plt.text(0.5,-1.85,'$RMSE = %1.3f$'%rmse,ha='left',va='top')
  plt.xlim([-0.1,0.7])
  plt.ylim([-2,2])
  plt.title('%s %s %s, %i-%i -> %i-%i'%(datasource,temp_product,SREX_region,yr_start1,yr_end1,yr_start2,yr_end2))
  plt.show(block=False)
elif plot_opt == 'all':
  df_plot = []
  dT_plot = []
  
  df_plot.append(df_irr.ravel()[(~np.isnan(df_irr.ravel())) & (~np.isnan(dT.ravel()))])
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
  
  slope_ub = []
  slope_lb = []
  slopes = [stats.linregress(x,y)[3] for x,y in zip(df_plot,dT_plot)]
  cod = [stats.linregress(x,y)[2]**2 for x,y in zip(df_plot,dT_plot)]
  
  for i in range(1,8):
    slope_ub.append(stats.theilslopes(df_plot[i],dT_plot[i],alpha=1-p_value)[3])
    slope_lb.append(stats.theilslopes(df_plot[i],dT_plot[i],alpha=1-p_value)[2])
    
  f, axs = plt.subplots(2,4,sharex='col',sharey='row',figsize=(21,6))
  
  SREX_list = ['Global land','CNA','EAS','MED','SAS+','SEA','WAS','WNA']
  for i in range(0,2):
    for j in range(0,4):
      scat = axs[i,j].scatter(df_plot[4*i+j],dT_plot[4*i+j],color='black',s=2,label='data points')
      if slopes[4*i+j] < p_value:
        axs[i,j].plot(np.unique(df_plot[4*i+j]), np.poly1d(np.polyfit(df_plot[4*i+j],dT_plot[4*i+j],1))(np.unique(df_plot[4*i+j])),'g-',linewidth=2,label='1st order fit')
        #axs[i,j].plot(np.unique(df_plot[4*i+j]), np.poly1d(np.polyfit(df_plot[4*i+j],dT_plot[4*i+j],2))(np.unique(df_plot[4*i+j])),'g--',linewidth=2,label='2nd order fit')
      else:                                                                                                               
        axs[i,j].plot(np.unique(df_plot[4*i+j]), np.poly1d(np.polyfit(df_plot[4*i+j],dT_plot[4*i+j],1))(np.unique(df_plot[4*i+j])),'b-',label='1st order fit')
        #axs[i,j].plot(np.unique(df_plot[4*i+j]), np.poly1d(np.polyfit(df_plot[4*i+j],dT_plot[4*i+j],2))(np.unique(df_plot[4*i+j])),'b--',label='2nd order fit')
      axs[i,j].set_title(SREX_list[4*i+j])
      if yr_start1 < 1950:
        axs[i,j].set_xlim([-0.4,0.8])
        axs[i,j].set_ylim([-3,3])
        axs[i,j].annotate('n=%i, $r^2$=%1.2f'%(df_plot[4*i+j].shape[0],cod[4*i+j]),xy=(0.78,-2.8),ha='right',va='bottom')
      else:
        axs[i,j].set_xlim([-0.2,0.4])
        axs[i,j].set_ylim([-3,3])
        axs[i,j].annotate('n=%i, $r^2$=%1.2f'%(df_plot[4*i+j].shape[0],cod[4*i+j]),xy=(0.39,-2.8),ha='right',va='bottom')
      axs[i,j].grid(True)
      
  #Plot axis labels
  f.text(0.5125, 0.04, '$\Delta f_{irr}$ [-]', ha='center', va='center',fontsize=12)
  f.text(0.1, 0.5, '$\Delta T$ [$^\circ$C]', ha='center', va='center', rotation='vertical',fontsize=12)
  f.text(0.5, 0.96, '%s %s, %i-%i -> %i-%i, PD thres: %1.2f'%(datasource,temp_product,yr_start1,yr_end1,yr_start2,yr_end2,thres_irr_PD),ha='center',va='center',fontsize=14,weight='bold')
  
  #Plot legends
  linline = mlines.Line2D([], [], color='blue',label='linear fit')
  quadline = mlines.Line2D([], [], color='blue',linestyle='--',label='quadratic fit')
  f.legend((scat,linline,quadline),('data points','linear fit','quadratic fit'),loc=7)
  
  plt.show(block=False)
  
#Save results  
if figsave == True:
  if plot_opt == 'region':
    figpath = '/net/exo/landclim/wthiery/observational_analysis/Figures/%s/%s/scatterplots/scatter.%s.%s.%s.%i-%i.%i.%i.%s'%(datasource,temp_product,datasource,temp_product,SREX_region,yr_start1,yr_end1,yr_start2,yr_end2,figformat)
  elif plot_opt == 'all':
    figpath = '/net/exo/landclim/wthiery/observational_analysis/Figures/%s/%s/scatterplots/scatter.%s.%s.all.%1.2f.%i-%i.%i-%i.%s'%(datasource,temp_product,datasource,temp_product,thres_irr_PD,yr_start1,yr_end1,yr_start2,yr_end2,figformat)
  
  if os.path.exists(figpath):
    print('Figure already exists at path: ' + figpath)
  else:
    print('Saving figure as:' + figpath)
    plt.savefig(figpath,bbox_inches='tight')
    
  
  
