"""calc_coolfrac.py

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

execfile('calc_irr_impact_thres.py')
execfile('calc_irr_impact_regr.py')

#################################
#User-specified options
#################################
datasource = "CRU"
temp_product = "tmx_max"
response = "PD"
t_res = "seasonal"
u = 0.5

yr_start1 = 1901
yr_end1 = 1930
yr_start2 = 1981
yr_end2 = 2010

thres_irr_max = 0.8

figsave = True
figformat = 'pdf'
#################################

#fcool_regdf = np.zeros((10))
fcool_regndf = np.zeros((int(thres_irr_max/0.05)))
fcool_thres = np.zeros((int(thres_irr_max/0.05)))
for i in range(0,int(thres_irr_max/0.05)):
  print('Currently processing threshold: ' + str((i+1)/20.))
  
  #dT_irr_regdf, _,_,_,_ = calc_irr_impact_regr(datasource,temp_product,response,t_res,(i+1)/20.,True,yr_start1,yr_end1,yr_start2,yr_end2)
  dT_irr_reg, _,_,_,_ = calc_irr_impact_regr(datasource,temp_product,response,t_res,(i+1)/20.,False,yr_start1,yr_end1,yr_start2,yr_end2)
  dT_irr_thres, _,_ = calc_irr_impact_thres(datasource,temp_product,response,t_res,(i+1)/20.,u,yr_start1,yr_end1,yr_start2,yr_end2)
  
  #Calculate fraction of cooled cells
  #fcool_regdf[i] = np.count_nonzero(dT_irr_regdf < 0.0) / np.float(np.count_nonzero(~np.isnan(dT_irr_regdf)))
  fcool_regndf[i] = np.count_nonzero(dT_irr_reg < 0.0) / np.float(np.count_nonzero(~np.isnan(dT_irr_reg)))
  fcool_thres[i] = np.count_nonzero(dT_irr_thres < 0.0) / np.float(np.count_nonzero(~np.isnan(dT_irr_thres)))
  
  #del dT_irr_regdf,dT_irr_reg,dT_irr_thres
  del dT_irr_reg,dT_irr_thres
  
plt.figure(figsize=(12,8))
plt.plot(np.arange(0.05,thres_irr_max+0.05,0.05),fcool_thres,'ro-',label='TB')
plt.plot(np.arange(0.05,thres_irr_max+0.05,0.05),fcool_regndf,'bo-',label='RB')
plt.xlabel('$f_{irr,PD}$ [-]')
plt.ylabel('Fraction of cooled cells [-]')
plt.xlim([0,thres_irr_max+0.05])
#plt.ylim([0.5,1])
plt.grid(True)
plt.legend(loc='lower right')
plt.title('Data source: %s %s, ref: %s-%s'%(datasource,temp_product,yr_start1,yr_end1))
plt.show(block=False)

if figsave == True:
  figpath = '/net/exo/landclim/wthiery/observational_analysis/Figures/%s/%s/GIL/coolfrac.%s.%s.%s.%i-%i.%i-%i.%s'%(datasource,temp_product,datasource,temp_product,response,yr_start1,yr_end1,yr_start2,yr_end2,figformat)
  if os.path.exists(figpath):
    print('Figure already exists at path: ' + figpath)
  else:
    print('Saving figure as: ' + figpath)
    plt.savefig(figpath,bbox_inches='tight')
