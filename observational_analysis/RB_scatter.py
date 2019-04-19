"""RB_scatter.py

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
import matplotlib.colors as colors

execfile('calc_irr_impact_thres.py')
execfile('calc_irr_impact_regr.py')

#################################
#User-specified options
#################################
datasource = "CRU"
temp_product = "tmx"
response = "PD"
t_res = "seasonal"
method = "threshold"
thres_irr = 0.25

yr_start1 = 1901
yr_end1 = 1930
yr_start2 = 1981
yr_end2 = 2010

#figsave = True
#figformat = 'png'
#################################

_, indvals, codvals = calc_irr_impact_regr(datasource,temp_product,response,t_res,thres_irr,False,yr_start1,yr_end1,yr_start2,yr_end2)

irrvals = indvals[0][0]
Tvals = indvals[0][1]

for i in range(1,len(indvals)):
  irrvals = np.append(irrvals,indvals[i][0],axis=0)
  Tvals = np.append(Tvals,indvals[i][1],axis=0)
  
s = stats.linregress(irrvals,Tvals).slope
i = stats.linregress(irrvals,Tvals).intercept

plt.figure(figsize=(12,8))
plt.scatter(irrvals,Tvals)
plt.plot([-0.2,1],[-0.2*s+i,s+i],'k--')
plt.xlabel('PD $f_{irr}$ [-]')
plt.ylabel('$\Delta T_{irr}$ [$^\circ$C]')
plt.title('Dataset: %s %s %s %1.2f'%(datasource,temp_product,response,thres_irr))

cc = np.corrcoef(irrvals,Tvals)[1,0]
plt.text(0.8,0,'$r^2 = %1.4f$'%cc)

plt.show(block=False)