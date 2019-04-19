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
SREX_region = "GIL"
temp_product = "tmx_max"
response = "PD"
t_res = "seasonal"
method = "regression"
scaling = False
thres_irr = 0.25

u = 0.5

p_value = 0.01
figsave = True
figformat = 'pdf'
#################################

#Perform calculations for different threshold values

if datasource in ["CRU","CESM","CRU_CESM"]:
  if method == "regression":
    t1, _,_,_,_ = calc_irr_impact_regr(datasource,temp_product,response,t_res,thres_irr,False,1901,1930,1981,2010)
    t2, _,_,_,_ = calc_irr_impact_regr(datasource,temp_product,response,t_res,thres_irr,False,1911,1940,1981,2010)
    t3, _,_,_,_ = calc_irr_impact_regr(datasource,temp_product,response,t_res,thres_irr,False,1921,1950,1981,2010)
    t4, _,_,_,_ = calc_irr_impact_regr(datasource,temp_product,response,t_res,thres_irr,False,1931,1960,1981,2010)
    t5, _,_,_,_ = calc_irr_impact_regr(datasource,temp_product,response,t_res,thres_irr,False,1941,1970,1981,2010)
    t6, _,_,_,_ = calc_irr_impact_regr(datasource,temp_product,response,t_res,thres_irr,False,1951,1980,1981,2010)
  elif method == "threshold":
    if scaling == True:
      _,_, t1 = calc_irr_impact_thres(datasource,temp_product,response,t_res,thres_irr,u,1901,1930,1981,2010)
      _,_, t2 = calc_irr_impact_thres(datasource,temp_product,response,t_res,thres_irr,u,1911,1940,1981,2010)
      _,_, t3 = calc_irr_impact_thres(datasource,temp_product,response,t_res,thres_irr,u,1921,1950,1981,2010)
      _,_, t4 = calc_irr_impact_thres(datasource,temp_product,response,t_res,thres_irr,u,1931,1960,1981,2010)
      _,_, t5 = calc_irr_impact_thres(datasource,temp_product,response,t_res,thres_irr,u,1941,1970,1981,2010)
      _,_, t6 = calc_irr_impact_thres(datasource,temp_product,response,t_res,thres_irr,u,1951,1980,1981,2010)
    elif scaling == False:
      t1, _,_ = calc_irr_impact_thres(datasource,temp_product,response,t_res,thres_irr,u,1901,1930,1981,2010)
      t2, _,_ = calc_irr_impact_thres(datasource,temp_product,response,t_res,thres_irr,u,1911,1940,1981,2010)
      t3, _,_ = calc_irr_impact_thres(datasource,temp_product,response,t_res,thres_irr,u,1921,1950,1981,2010)
      t4, _,_ = calc_irr_impact_thres(datasource,temp_product,response,t_res,thres_irr,u,1931,1960,1981,2010)
      t5, _,_ = calc_irr_impact_thres(datasource,temp_product,response,t_res,thres_irr,u,1941,1970,1981,2010)
      t6, _,_ = calc_irr_impact_thres(datasource,temp_product,response,t_res,thres_irr,u,1951,1980,1981,2010)

#Change to 1D, remove NaNs
if datasource == 'CRU':
  if SREX_region == 'GIL':
    out1 = t1.ravel()[~np.isnan(t1.ravel())]
    out2 = t2.ravel()[~np.isnan(t2.ravel())]
    out3 = t3.ravel()[~np.isnan(t3.ravel())]
    out4 = t4.ravel()[~np.isnan(t4.ravel())]
    out5 = t5.ravel()[~np.isnan(t5.ravel())]
    out6 = t6.ravel()[~np.isnan(t6.ravel())]
  elif SREX_region == 'CNA':
    out1 = t1[237:280,150:190].ravel()[~np.isnan(t1[237:280,150:190].ravel())]
    out2 = t2[237:280,150:190].ravel()[~np.isnan(t2[237:280,150:190].ravel())]
    out3 = t3[237:280,150:190].ravel()[~np.isnan(t3[237:280,150:190].ravel())]
    out4 = t4[237:280,150:190].ravel()[~np.isnan(t4[237:280,150:190].ravel())]
    out5 = t5[237:280,150:190].ravel()[~np.isnan(t5[237:280,150:190].ravel())]
    out6 = t6[237:280,150:190].ravel()[~np.isnan(t6[237:280,150:190].ravel())]
  elif SREX_region == 'EAS':
    out1 = t1[220:280,560:650].ravel()[~np.isnan(t1[220:280,560:650].ravel())]
    out2 = t2[220:280,560:650].ravel()[~np.isnan(t2[220:280,560:650].ravel())]
    out3 = t3[220:280,560:650].ravel()[~np.isnan(t3[220:280,560:650].ravel())]
    out4 = t4[220:280,560:650].ravel()[~np.isnan(t4[220:280,560:650].ravel())]
    out5 = t5[220:280,560:650].ravel()[~np.isnan(t5[220:280,560:650].ravel())]
    out6 = t6[220:280,560:650].ravel()[~np.isnan(t6[220:280,560:650].ravel())]
  elif SREX_region == 'MED':
    out1 = t1[240:270,340:440].ravel()[~np.isnan(t1[240:270,340:440].ravel())]
    out2 = t2[240:270,340:440].ravel()[~np.isnan(t2[240:270,340:440].ravel())]
    out3 = t3[240:270,340:440].ravel()[~np.isnan(t3[240:270,340:440].ravel())]
    out4 = t4[240:270,340:440].ravel()[~np.isnan(t4[240:270,340:440].ravel())]
    out5 = t5[240:270,340:440].ravel()[~np.isnan(t5[240:270,340:440].ravel())]
    out6 = t6[240:270,340:440].ravel()[~np.isnan(t6[240:270,340:440].ravel())]
  elif SREX_region == 'SAS':
    out1 = t1[170:270,460:580].ravel()[~np.isnan(t1[170:270,460:580].ravel())]
    out2 = t2[170:270,460:580].ravel()[~np.isnan(t2[170:270,460:580].ravel())]
    out3 = t3[170:270,460:580].ravel()[~np.isnan(t3[170:270,460:580].ravel())]
    out4 = t4[170:270,460:580].ravel()[~np.isnan(t4[170:270,460:580].ravel())]
    out5 = t5[170:270,460:580].ravel()[~np.isnan(t5[170:270,460:580].ravel())]
    out6 = t6[170:270,460:580].ravel()[~np.isnan(t6[170:270,460:580].ravel())]
  elif SREX_region == 'SEA':
    out1 = t1[160:220,550:670].ravel()[~np.isnan(t1[160:220,550:670].ravel())]
    out2 = t2[160:220,550:670].ravel()[~np.isnan(t2[160:220,550:670].ravel())]
    out3 = t3[160:220,550:670].ravel()[~np.isnan(t3[160:220,550:670].ravel())]
    out4 = t4[160:220,550:670].ravel()[~np.isnan(t4[160:220,550:670].ravel())]
    out5 = t5[160:220,550:670].ravel()[~np.isnan(t5[160:220,550:670].ravel())]
    out6 = t6[160:220,550:670].ravel()[~np.isnan(t6[160:220,550:670].ravel())]
  elif SREX_region == 'WAS':
    out1 = t1[210:280,440:480].ravel()[~np.isnan(t1[210:280,440:480].ravel())]
    out2 = t2[210:280,440:480].ravel()[~np.isnan(t2[210:280,440:480].ravel())]
    out3 = t3[210:280,440:480].ravel()[~np.isnan(t3[210:280,440:480].ravel())]
    out4 = t4[210:280,440:480].ravel()[~np.isnan(t4[210:280,440:480].ravel())]
    out5 = t5[210:280,440:480].ravel()[~np.isnan(t5[210:280,440:480].ravel())]
    out6 = t6[210:280,440:480].ravel()[~np.isnan(t6[210:280,440:480].ravel())]
  elif SREX_region == 'WNA':
    out1 = t1[237:300,100:150].ravel()[~np.isnan(t1[237:300,100:150].ravel())]
    out2 = t2[237:300,100:150].ravel()[~np.isnan(t2[237:300,100:150].ravel())]
    out3 = t3[237:300,100:150].ravel()[~np.isnan(t3[237:300,100:150].ravel())]
    out4 = t4[237:300,100:150].ravel()[~np.isnan(t4[237:300,100:150].ravel())]
    out5 = t5[237:300,100:150].ravel()[~np.isnan(t5[237:300,100:150].ravel())]
    out6 = t6[237:300,100:150].ravel()[~np.isnan(t6[237:300,100:150].ravel())]

#Change into format readable by boxplot, removing entries without data if needed
data = [out1, out2, out3, out4, out5, out6]
dataplot = [x for x in data if x.shape[0] > 0]
lenlist = [x.shape[0] for x in dataplot]
s_max = max([x.shape[0] for x in dataplot])

#Plot results
plt.figure(figsize=(12,8))
#bp = plt.boxplot(dataplot,widths = [0.9 * (q / np.float(s_max)) for q in lenlist],patch_artist=True)
bp = plt.boxplot(dataplot,patch_artist=True)
for a in range(0,len(dataplot)):
  p_ttest = stats.ttest_1samp(dataplot[a],0)[1]
  if dataplot[a].shape[0] >= 3:
    p_shapiro = stats.shapiro(dataplot[a])[1]
  if p_ttest <= p_value:
    bp['boxes'][a].set(edgecolor='green',linewidth=3)
  if p_shapiro <= p_value:
    bp['boxes'][a].set(hatch='//')
  bp['boxes'][a].set(facecolor='white')

plt.xticks(np.arange(1,7),['1901-1930','1911-1940','1921-1950','1931-1960','1941-1970','1951-1980'])
plt.xlabel('Reference periods [-]')
plt.ylabel('$\Delta T_{irr}$ [$^\circ$C]')
plt.ylim([-2,2])
for i in range(0,len(dataplot)):
  plt.text(i+1,-1.9,str(lenlist[i]),horizontalalignment='center',verticalalignment='center')
meanpatch = mpatches.Patch(edgecolor='green',facecolor='white',linewidth=3,label=r'reject $H_0: \mu = 0$')
normalpatch = mpatches.Patch(facecolor='white',linewidth=0,hatch='//',label=r'reject $H_0: X \sim N(\mu,\sigma^2)$')
plt.legend(handles=[meanpatch,normalpatch],loc='upper right')
plt.plot(np.arange(1,7),[x.mean() for x in data],'k--')
plt.grid(True)
if method == 'threshold':
  plt.title('method: %s %s, dataset: %s %s, region: %s, p: %s, thres_irr: %1.2f, u: %1.2f'%(method, response, datasource,temp_product,SREX_region,p_value,thres_irr,u))
elif method == 'regression':
  plt.title('method: %s %s, dataset: %s %s, region: %s, p: %s, thres_irr: %1.2f'%(method, response, datasource,temp_product,SREX_region,p_value,thres_irr))
plt.show(block=False)

if figsave == True:
  if response == 'PC/PD':
    resp = 'PC-PD'
  elif response == 'PD':
    resp = 'PD'
  
  if method == 'threshold':
    figpath = '/net/exo/landclim/wthiery/observational_analysis/Figures/%s/boxplots/%s/dT_irr.refperiod.%s.%s.%s.%s.%s.%1.2f.%1.2f.%s'%(datasource,temp_product,datasource,temp_product,method,scaling,resp,thres_irr,figformat)
  elif method == 'regression':
    figpath = '/net/exo/landclim/wthiery/observational_analysis/Figures/%s/boxplots/%s/dT_irr.refperiod.%s.%s.%s.%s.%s.%1.2f.%s'%(datasource,temp_product,datasource,temp_product,method,scaling,resp,thres_irr,figformat)
  
  if os.path.exists(figpath):
    print('Figure already exists at path: ' + figpath)
  else:
    print('Saving figure as: ' + figpath)
    plt.savefig(figpath,bbox_inches='tight')
    
#del t1,t2,t3,t4,t5,t6,t7,t8,out1,out2,out3,out4,out5,out6,out7,out8,data
