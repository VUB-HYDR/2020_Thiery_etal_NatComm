"""monthlyboxplots.py

author: Auke Visser
date: 29.11.2016

This script calls a routine to calculate the irrigation impact on 
temperature with either the threshold- or the regression-based
approach for HadEX2 and plots monthly boxplots

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

execfile('calc_irr_impact_thres.py')
execfile('calc_irr_impact_regr.py')

####################################
#######User-specified options#######
datasource = 'HadEX2'
method = 'regression'
response = 'PD'
thres_irr = 0.2
u = 0.5
t_res = 'monthly'

yr_start1 = 1951
yr_end1 = 1980
yr_start2 = 1981
yr_end2 = 2010

p_value = 0.01
def_regr = False #Not implemented for HadEX2

figsave = True
figformat = 'png'

####################################

monthlist = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec', 'Ann']
data = []

_,_,f_irr = calc_irr_diff(datasource,response,yr_start1,yr_end1,yr_start2,yr_end2)
df_irr = f_irr[1,:,:] - f_irr[0,:,:]

#Calculate dT_irr for every month using the TB- or RB-approach
for i in range(0,len(monthlist)):
  if method == 'regression':
    dT_irr, _,_,_,_ = calc_irr_impact_regr(datasource,monthlist[i],response,t_res,thres_irr,def_regr,yr_start1,yr_end1,yr_start2,yr_end2)
    data.append(dT_irr[(~np.isnan(dT_irr)) & (df_irr >= thres_irr)])
  elif method == 'threshold':
    dT_irr, _,_ = calc_irr_impact_thres(datasource,monthlist[i],response,t_res,thres_irr,u,yr_start1,yr_end1,yr_start2,yr_end2)
    data.append(dT_irr[~np.isnan(dT_irr)])
  
#Plot the results
plt.figure(figsize=(12,8))
bp = plt.boxplot(data,patch_artist=True)

#Perform statistical tests
for a in range(0,len(data)):
  p_ttest = stats.ttest_1samp(data[a],0)[1]
  if data[a].shape[0] >= 3:
    p_shapiro = stats.shapiro(data[a])[1]
    if p_shapiro <= p_value:
      bp['boxes'][a].set(hatch='//')
  if p_ttest <= p_value:
    bp['boxes'][a].set(edgecolor='green',linewidth=3)
  bp['boxes'][a].set(facecolor='white')

#Some plotting options
plt.xticks(np.arange(1,14),monthlist)
plt.axvline(x=12.5, color = 'black')
plt.xlabel('$f_{irr}$ [-]')
plt.ylabel('$\Delta T_{irr}$ [$^\circ$C]')
meanpatch = mpatches.Patch(edgecolor='green',facecolor='white',linewidth=3,label=r'reject $H_0: \mu = 0$')
normalpatch = mpatches.Patch(facecolor='white',linewidth=0,hatch='//',label=r'reject $H_0: X \sim N(\mu,\sigma^2)$')
plt.legend(handles=[meanpatch,normalpatch],loc='upper right')
plt.grid(True)
plt.title('method: %s %s, dataset: %s, threshold: %1.2f, region: GIL, p: %s, ref: %i-%i'%(method,response,datasource,thres_irr,p_value,yr_start1,yr_end1))
plt.show(block=False)

#Save result
if figsave == True:
  if response == 'PC/PD':
    resp = 'PC-PD'
  elif response == 'PD':
    resp = 'PD'
  figpath = '/net/exo/landclim/wthiery/observational_analysis/Figures/%s/boxplots/dT_irr.monthly.%s.%s.%s.%1.2f.%i-%i.%i-%i.%s'%(datasource,datasource,method,resp,thres_irr,yr_start1,yr_end1,yr_start2,yr_end2,figformat)
  if os.path.exists(figpath):
    print('Figure already exists at path: ' + figpath)
  else:
    print('Saving figure as: ' + figpath)
    plt.savefig(figpath,bbox_inches='tight')
