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
import matplotlib.lines as mlines
import matplotlib.colors as colors

execfile('calc_irr_impact_regr.py')
execfile('calc_irr_impact_thres.py')

#################################
#User-specified options
#################################
datasource = "CRU_CESM"
temp_product_CRU = "tmp_max"
temp_product_CESM = "TREFHT"
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
figformat = 'pdf'
#################################

#Perform calculations for different threshold values
_,_,f_irr = calc_irr_diff('CRU_CESM',response,yr_start1,yr_end1,yr_start2,yr_end2)

df_irr = (f_irr[1] - f_irr[0]) * 100.
df_irr[df_irr == 0] = np.nan
df_irr[f_irr[1] < 0.1] = np.nan

CESM_irrfile = nc.Dataset('/home/wthiery/documents/research/cesm_present/postprocessing/constants/f.e122.F1850PDC5.f09_g16.irrigation-io192.001_constants.nc','r')
f_irr_CESM = CESM_irrfile.variables['PCT_IRRIG'][:]
f_irr_CESM[f_irr_CESM == 0] = np.nan
f_irr_CESM[f_irr_CESM < 0.1] = np.nan

#Perform calculations for CRU and CESM for different thresholds
CRU_t1, _,_,_,_ = calc_irr_impact_regr('CRU_CESM',temp_product_CRU,response,t_res,0.1,False,yr_start1,yr_end1,yr_start2,yr_end2)
CRU_t2, _,_,_,_ = calc_irr_impact_regr('CRU_CESM',temp_product_CRU,response,t_res,0.2,False,yr_start1,yr_end1,yr_start2,yr_end2)
CRU_t3, _,_,_,_ = calc_irr_impact_regr('CRU_CESM',temp_product_CRU,response,t_res,0.3,False,yr_start1,yr_end1,yr_start2,yr_end2)
CRU_t4, _,_,_,_ = calc_irr_impact_regr('CRU_CESM',temp_product_CRU,response,t_res,0.4,False,yr_start1,yr_end1,yr_start2,yr_end2)
CRU_t5, _,_,_,_ = calc_irr_impact_regr('CRU_CESM',temp_product_CRU,response,t_res,0.5,False,yr_start1,yr_end1,yr_start2,yr_end2)
CRU_t6, _,_,_,_ = calc_irr_impact_regr('CRU_CESM',temp_product_CRU,response,t_res,0.6,False,yr_start1,yr_end1,yr_start2,yr_end2)
CRU_t7, _,_,_,_ = calc_irr_impact_regr('CRU_CESM',temp_product_CRU,response,t_res,0.7,False,yr_start1,yr_end1,yr_start2,yr_end2)
CRU_t8, _,_,_,_ = calc_irr_impact_regr('CRU_CESM',temp_product_CRU,response,t_res,0.8,False,yr_start1,yr_end1,yr_start2,yr_end2)

CESM_t1, _,_,_,_ = calc_irr_impact_regr('CESM',temp_product_CESM,response,t_res,0.1,False,yr_start1,yr_end1,yr_start2,yr_end2)
CESM_t2, _,_,_,_ = calc_irr_impact_regr('CESM',temp_product_CESM,response,t_res,0.2,False,yr_start1,yr_end1,yr_start2,yr_end2)
CESM_t3, _,_,_,_ = calc_irr_impact_regr('CESM',temp_product_CESM,response,t_res,0.3,False,yr_start1,yr_end1,yr_start2,yr_end2)
CESM_t4, _,_,_,_ = calc_irr_impact_regr('CESM',temp_product_CESM,response,t_res,0.4,False,yr_start1,yr_end1,yr_start2,yr_end2)
CESM_t5, _,_,_,_ = calc_irr_impact_regr('CESM',temp_product_CESM,response,t_res,0.5,False,yr_start1,yr_end1,yr_start2,yr_end2)
CESM_t6, _,_,_,_ = calc_irr_impact_regr('CESM',temp_product_CESM,response,t_res,0.6,False,yr_start1,yr_end1,yr_start2,yr_end2)
CESM_t7, _,_,_,_ = calc_irr_impact_regr('CESM',temp_product_CESM,response,t_res,0.7,False,yr_start1,yr_end1,yr_start2,yr_end2)
CESM_t8, _,_,_,_ = calc_irr_impact_regr('CESM',temp_product_CESM,response,t_res,0.8,False,yr_start1,yr_end1,yr_start2,yr_end2)

#Change to 1D, remove NaNs
CRU1 = CRU_t1.ravel()[~np.isnan(CRU_t1.ravel())]
CRU2 = CRU_t2.ravel()[~np.isnan(CRU_t2.ravel())]
CRU3 = CRU_t3.ravel()[~np.isnan(CRU_t3.ravel())]
CRU4 = CRU_t4.ravel()[~np.isnan(CRU_t4.ravel())]
CRU5 = CRU_t5.ravel()[~np.isnan(CRU_t5.ravel())]
CRU6 = CRU_t6.ravel()[~np.isnan(CRU_t6.ravel())]
CRU7 = CRU_t7.ravel()[~np.isnan(CRU_t7.ravel())]
CRU8 = CRU_t8.ravel()[~np.isnan(CRU_t8.ravel())]

CESM1 = CESM_t1.ravel()[~np.isnan(CESM_t1.ravel())]
CESM2 = CESM_t2.ravel()[~np.isnan(CESM_t2.ravel())]
CESM3 = CESM_t3.ravel()[~np.isnan(CESM_t3.ravel())]
CESM4 = CESM_t4.ravel()[~np.isnan(CESM_t4.ravel())]
CESM5 = CESM_t5.ravel()[~np.isnan(CESM_t5.ravel())]
CESM6 = CESM_t6.ravel()[~np.isnan(CESM_t6.ravel())]
CESM7 = CESM_t7.ravel()[~np.isnan(CESM_t7.ravel())]
CESM8 = CESM_t8.ravel()[~np.isnan(CESM_t8.ravel())]

#Change into list format readable by boxplot, removing entries without data if needed
CRU_data = [CRU1, CRU2, CRU3, CRU4, CRU5, CRU6, CRU7, CRU8]
CRU_dataplot = [x for x in CRU_data if x.shape[0] > 0]
CRU_lenlist = [x.shape[0] for x in CRU_dataplot]

CESM_data = [CESM1, CESM2, CESM3, CESM4, CESM5, CESM6, CESM7, CESM8]
CESM_dataplot = [x[(x > -999) & (x < 999)] for x in CESM_data if x.shape[0] > 0]
CESM_lenlist = [x.shape[0] for x in CESM_dataplot]

#Plot results
plt.figure(figsize=(12,8))
pos = np.arange(1,9)

bp1 = plt.boxplot(CRU_dataplot,'',positions=pos-0.2,widths=0.4,patch_artist=True)
bp2 = plt.boxplot(CESM_dataplot,'',positions=pos+0.2,widths=0.4,patch_artist=True)

for a in range(0,len(CRU_dataplot)):
  p_ttest_CRU = stats.ttest_1samp(CRU_dataplot[a],0)[1]
  p_ttest_CESM = stats.ttest_1samp(CESM_dataplot[a],0)[1]
  
  #Perform test for normality
  if CRU_dataplot[a].shape[0] >= 3:
    p_shapiro_CRU = stats.shapiro(CRU_dataplot[a])[1]
  if CESM_dataplot[a].shape[0] >= 3:
    p_shapiro_CESM = stats.shapiro(CESM_dataplot[a])[1]
  
  #Test if the mean significantly differs from 0
  if p_ttest_CRU <= p_value:
    bp1['boxes'][a].set(edgecolor='black',linewidth=2)
  if p_shapiro_CRU <= p_value:
    bp1['boxes'][a].set(hatch='//')
  if p_ttest_CESM <= p_value:
    bp2['boxes'][a].set(edgecolor='black',linewidth=2)
  if p_shapiro_CESM <= p_value:
    bp2['boxes'][a].set(hatch='//')
  
  bp1['boxes'][a].set(facecolor='tan')
  bp2['boxes'][a].set(facecolor='lightblue')

  plt.text(a+0.8,-2.9,str(CRU_lenlist[a]),horizontalalignment='center',verticalalignment='center')
  plt.text(a+1.2,-2.9,str(CESM_lenlist[a]),horizontalalignment='center',verticalalignment='center')

#Some additional plotting options
plt.xticks(np.arange(1,len(CRU_dataplot)+1),np.arange(0.1,len(CRU_dataplot)/10.+0.1,0.1))
plt.xlabel('$f_{irr,PD}$ [-]')
plt.ylabel(r'$\Delta T_{irr} \left[^\circ C\right]$')
plt.xlim([0.5,len(CRU_dataplot)+.5])
plt.ylim([-3,1])
plt.title('CRU: %s - CESM: %s'%(temp_product_CRU,temp_product_CESM))
CRUpatch = mpatches.Patch(edgecolor='black',facecolor='tan',linewidth=1,label=r'CRU')
CESMpatch = mpatches.Patch(edgecolor='black',facecolor='lightblue',linewidth=1,label=r'CESM')
meanpatch = mpatches.Patch(edgecolor='black',facecolor='white',linewidth=2,label=r'reject $H_0: \mu = 0$')
normalpatch = mpatches.Patch(facecolor='white',linewidth=0,hatch='//',label=r'reject $H_0: X \sim N(\mu,\sigma^2)$')
plt.legend(handles=[CRUpatch,CESMpatch,meanpatch,normalpatch],loc='upper right')
plt.grid(True)

#Save figure
if figsave == True:
  if response == 'PC/PD':
    resp = 'PC-PD'
  elif response == 'PD':
    resp = 'PD'
  
  figpath = '/net/exo/landclim/wthiery/observational_analysis/Figures/CRU_CESM_comp/dT_irr.thresvals.%s.%s.%s.%i-%i.%i-%i.%s'%(temp_product_CRU,method,resp,yr_start1,yr_end1,yr_start2,yr_end2,figformat)
  
  if os.path.exists(figpath):
    print('Figure already exists at path: ' + figpath)
  else:
    print('Saving figure as: ' + figpath)
    plt.savefig(figpath,bbox_inches='tight')
    
#del t1,t2,t3,t4,t5,t6,t7,t8,out1,out2,out3,out4,out5,out6,out7,out8,data
