"""dT_asfuncof_thres.py

author: Auke Visser
date: 23.11.2016

This script calls a routine to calculate the irrigation impact on the sign 
of the T change with either the threshold- or the regression-based approach
The user has the option to print or visualize output.

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
import matplotlib.colors as colors

execfile('calc_irr_impact_regr.py')
execfile('calc_irr_impact_thres.py')

#################################
#User-specified options
#################################
datasource = "CRU"
temp_product = "tmx_max"
response = "PD"
t_res = "seasonal"
method = "regression"
scaling = False

yr_start1 = 1901
yr_end1 = 1930
yr_start2 = 1981
yr_end2 = 2010

u = 0.5

p_value = 0.01
figsave = False
figformat = 'png'
#################################

#Perform calculations for different threshold values
#The method is executed multiple times for different threshold values
_,_,f_irr = calc_irr_diff(datasource,response,yr_start1,yr_end1,yr_start2,yr_end2)

if datasource in ["CRU","CESM","CRU_CESM"]:
  if method == "regression":
    
    t1, _,_,_,_ = calc_irr_impact_regr(datasource,temp_product,response,t_res,0.1,False,yr_start1,yr_end1,yr_start2,yr_end2)
    t2, _,_,_,_ = calc_irr_impact_regr(datasource,temp_product,response,t_res,0.2,False,yr_start1,yr_end1,yr_start2,yr_end2)
    t3, _,_,_,_ = calc_irr_impact_regr(datasource,temp_product,response,t_res,0.3,False,yr_start1,yr_end1,yr_start2,yr_end2)
    t4, _,_,_,_ = calc_irr_impact_regr(datasource,temp_product,response,t_res,0.4,False,yr_start1,yr_end1,yr_start2,yr_end2)
    t5, _,_,_,_ = calc_irr_impact_regr(datasource,temp_product,response,t_res,0.5,False,yr_start1,yr_end1,yr_start2,yr_end2)
    t6, _,_,_,_ = calc_irr_impact_regr(datasource,temp_product,response,t_res,0.6,False,yr_start1,yr_end1,yr_start2,yr_end2)
    t7, _,_,_,_ = calc_irr_impact_regr(datasource,temp_product,response,t_res,0.7,False,yr_start1,yr_end1,yr_start2,yr_end2)
    t8, _,_,_,_ = calc_irr_impact_regr(datasource,temp_product,response,t_res,0.8,False,yr_start1,yr_end1,yr_start2,yr_end2)
    
  elif method == "threshold":
    if scaling == True:
      _,_, t1 = calc_irr_impact_thres(datasource,temp_product,response,t_res,0.1,u,yr_start1,yr_end1,yr_start2,yr_end2)
      _,_, t2 = calc_irr_impact_thres(datasource,temp_product,response,t_res,0.2,u,yr_start1,yr_end1,yr_start2,yr_end2)
      _,_, t3 = calc_irr_impact_thres(datasource,temp_product,response,t_res,0.3,u,yr_start1,yr_end1,yr_start2,yr_end2)
      _,_, t4 = calc_irr_impact_thres(datasource,temp_product,response,t_res,0.4,u,yr_start1,yr_end1,yr_start2,yr_end2)
      _,_, t5 = calc_irr_impact_thres(datasource,temp_product,response,t_res,0.5,u,yr_start1,yr_end1,yr_start2,yr_end2)
      _,_, t6 = calc_irr_impact_thres(datasource,temp_product,response,t_res,0.6,u,yr_start1,yr_end1,yr_start2,yr_end2)
      _,_, t7 = calc_irr_impact_thres(datasource,temp_product,response,t_res,0.7,u,yr_start1,yr_end1,yr_start2,yr_end2)
      _,_, t8 = calc_irr_impact_thres(datasource,temp_product,response,t_res,0.8,u,yr_start1,yr_end1,yr_start2,yr_end2)
    elif scaling == False:
      t1, _,_ = calc_irr_impact_thres(datasource,temp_product,response,t_res,0.1,u,yr_start1,yr_end1,yr_start2,yr_end2)
      t2, _,_ = calc_irr_impact_thres(datasource,temp_product,response,t_res,0.2,u,yr_start1,yr_end1,yr_start2,yr_end2)
      t3, _,_ = calc_irr_impact_thres(datasource,temp_product,response,t_res,0.3,u,yr_start1,yr_end1,yr_start2,yr_end2)
      t4, _,_ = calc_irr_impact_thres(datasource,temp_product,response,t_res,0.4,u,yr_start1,yr_end1,yr_start2,yr_end2)
      t5, _,_ = calc_irr_impact_thres(datasource,temp_product,response,t_res,0.5,u,yr_start1,yr_end1,yr_start2,yr_end2)
      t6, _,_ = calc_irr_impact_thres(datasource,temp_product,response,t_res,0.6,u,yr_start1,yr_end1,yr_start2,yr_end2)
      t7, _,_ = calc_irr_impact_thres(datasource,temp_product,response,t_res,0.7,u,yr_start1,yr_end1,yr_start2,yr_end2)
      t8, _,_ = calc_irr_impact_thres(datasource,temp_product,response,t_res,0.8,u,yr_start1,yr_end1,yr_start2,yr_end2)
elif datasource == "HadEX2":
  if method == "regression":
    output, _,_,_,_ = calc_irr_impact_regr(datasource,temp_product,response,t_res,0.0,False,yr_start1,yr_end1,yr_start2,yr_end2)
    
    t1 = output[(f_irr[1,:,:] >= 0.05) & (~np.isnan(output))]
    t2 = output[(f_irr[1,:,:] >= 0.1) & (~np.isnan(output))]
    t3 = output[(f_irr[1,:,:] >= 0.15) & (~np.isnan(output))]
    t4 = output[(f_irr[1,:,:] >= 0.2) & (~np.isnan(output))]
    t5 = output[(f_irr[1,:,:] >= 0.25) & (~np.isnan(output))]
    t6 = output[(f_irr[1,:,:] >= 0.3) & (~np.isnan(output))]
    t7 = output[(f_irr[1,:,:] >= 0.35) & (~np.isnan(output))]
    t8 = output[(f_irr[1,:,:] >= 0.4) & (~np.isnan(output))]
  elif method == "threshold":
    if scaling == True:
      _,_, t1 = calc_irr_impact_thres(datasource,temp_product,response,t_res,u,0.05,yr_start1,yr_end1,yr_start2,yr_end2)
      _,_, t2 = calc_irr_impact_thres(datasource,temp_product,response,t_res,u,0.1,yr_start1,yr_end1,yr_start2,yr_end2)
      _,_, t3 = calc_irr_impact_thres(datasource,temp_product,response,t_res,u,0.15,yr_start1,yr_end1,yr_start2,yr_end2)
      _,_, t4 = calc_irr_impact_thres(datasource,temp_product,response,t_res,u,0.2,yr_start1,yr_end1,yr_start2,yr_end2)
      _,_, t5 = calc_irr_impact_thres(datasource,temp_product,response,t_res,u,0.25,yr_start1,yr_end1,yr_start2,yr_end2)
      _,_, t6 = calc_irr_impact_thres(datasource,temp_product,response,t_res,u,0.3,yr_start1,yr_end1,yr_start2,yr_end2)
      _,_, t7 = calc_irr_impact_thres(datasource,temp_product,response,t_res,u,0.35,yr_start1,yr_end1,yr_start2,yr_end2)
      _,_, t8 = calc_irr_impact_thres(datasource,temp_product,response,t_res,u,0.4,yr_start1,yr_end1,yr_start2,yr_end2)
    elif scaling == False:
      t1, _,_ = calc_irr_impact_thres(datasource,temp_product,response,t_res,u,0.05,yr_start1,yr_end1,yr_start2,yr_end2)
      t2, _,_ = calc_irr_impact_thres(datasource,temp_product,response,t_res,u,0.1,yr_start1,yr_end1,yr_start2,yr_end2)
      t3, _,_ = calc_irr_impact_thres(datasource,temp_product,response,t_res,u,0.15,yr_start1,yr_end1,yr_start2,yr_end2)
      t4, _,_ = calc_irr_impact_thres(datasource,temp_product,response,t_res,u,0.2,yr_start1,yr_end1,yr_start2,yr_end2)
      t5, _,_ = calc_irr_impact_thres(datasource,temp_product,response,t_res,u,0.25,yr_start1,yr_end1,yr_start2,yr_end2)
      t6, _,_ = calc_irr_impact_thres(datasource,temp_product,response,t_res,u,0.3,yr_start1,yr_end1,yr_start2,yr_end2)
      t7, _,_ = calc_irr_impact_thres(datasource,temp_product,response,t_res,u,0.35,yr_start1,yr_end1,yr_start2,yr_end2)
      t8, _,_ = calc_irr_impact_thres(datasource,temp_product,response,t_res,u,0.4,yr_start1,yr_end1,yr_start2,yr_end2)

#Change to 1D, remove NaNs
out1 = t1.ravel()[~np.isnan(t1.ravel())]
out2 = t2.ravel()[~np.isnan(t2.ravel())]
out3 = t3.ravel()[~np.isnan(t3.ravel())]
out4 = t4.ravel()[~np.isnan(t4.ravel())]
out5 = t5.ravel()[~np.isnan(t5.ravel())]
out6 = t6.ravel()[~np.isnan(t6.ravel())]
out7 = t7.ravel()[~np.isnan(t7.ravel())]
out8 = t8.ravel()[~np.isnan(t8.ravel())]

#Change into format readable by boxplot, removing entries without data if needed
data = [out1, out2, out3, out4, out5, out6, out7, out8]
dataplot = [x for x in data if x.shape[0] > 0]
lenlist = [x.shape[0] for x in dataplot]

#Plot results
plt.figure(figsize=(12,8))
bp = plt.boxplot(dataplot,patch_artist=True)

#Perform statistical tests on the mean and for normality
for a in range(0,len(dataplot)):
  p_ttest = stats.ttest_1samp(dataplot[a],0)[1]
  if dataplot[a].shape[0] >= 3:
    p_shapiro = stats.shapiro(dataplot[a])[1]
  if p_ttest <= p_value:
    bp['boxes'][a].set(edgecolor='green',linewidth=3)
  if p_shapiro <= p_value:
    bp['boxes'][a].set(hatch='//')
  bp['boxes'][a].set(facecolor='white')

#Plot the sample size at the bottom of the plot
if datasource in ["CRU","CESM","CRU_CESM"]:
  plt.xticks(np.arange(1,len(dataplot)+1),np.arange(0.1,len(dataplot)/10.+0.1,0.1))
elif datasource == "HadEX2":
  plt.xticks(np.arange(1,len(dataplot)+1),np.arange(0.05,len(dataplot)/20.+0.05,0.05))
  
#Some further plotting commands
plt.xlabel('$f_{irr,PD}$ [-]')
plt.ylabel(r'$\Delta T_{irr}$ [$^\circ$C]')
plt.ylim([-3,3])
for i in range(0,len(dataplot)):
  plt.text(i+1,-2.9,str(lenlist[i]),horizontalalignment='center',verticalalignment='center')
meanpatch = mpatches.Patch(edgecolor='green',facecolor='white',linewidth=3,label=r'reject $H_0: \mu = 0$')
normalpatch = mpatches.Patch(facecolor='white',linewidth=0,hatch='//',label=r'reject $H_0: X \sim N(\mu,\sigma^2)$')
plt.legend(handles=[meanpatch,normalpatch],loc='upper right')
plt.grid(True)
if method == 'threshold':
  plt.title('method: %s %s, dataset: %s %s, region: GIL, p: %s, ref: %i-%i, u: %1.2f'%(method, response, datasource,temp_product,p_value,yr_start1,yr_end1,u))
elif method == 'regression':
  plt.title('method: %s %s, dataset: %s %s, region: GIL, p: %s, ref: %i-%i'%(method, response, datasource,temp_product,p_value,yr_start1,yr_end1))
plt.show(block=False)

if figsave == True:
  if response == 'PC/PD':
    resp = 'PC-PD'
  elif response == 'PD':
    resp = 'PD'
  
  if method == 'threshold':
    figpath = '/net/exo/landclim/wthiery/observational_analysis/Figures/%s/boxplots/%s/dT_irr.thresvals.%s.%s.%s.%s.%s.%1.2f.%i-%i.%i-%i.%s'%(datasource,temp_product,datasource,temp_product,method,scaling,resp,u,yr_start1,yr_end1,yr_start2,yr_end2,figformat)
  elif method == 'regression':
    figpath = '/net/exo/landclim/wthiery/observational_analysis/Figures/%s/boxplots/%s/dT_irr.thresvals.%s.%s.%s.%s.%s.%i-%i.%i-%i.%s'%(datasource,temp_product,datasource,temp_product,method,scaling,resp,yr_start1,yr_end1,yr_start2,yr_end2,figformat)
  
  if os.path.exists(figpath):
    print('Figure already exists at path: ' + figpath)
  else:
    print('Saving figure as: ' + figpath)
    plt.savefig(figpath,bbox_inches='tight')
    
#del t1,t2,t3,t4,t5,t6,t7,t8,out1,out2,out3,out4,out5,out6,out7,out8,data
