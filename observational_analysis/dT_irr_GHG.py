"""dT_irr_GHG.py

author: Auke Visser
date: 5.12.2016



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
method = "threshold"

yr_start1 = 1901
yr_end1 = 1930
yr_start2 = 1981
yr_end2 = 2010

p_value = 0.01
figsave = False
figformat = 'png'
#################################

#Perform calculations for different threshold values
_,_,f_irr = calc_irr_diff(datasource,response,yr_start1,yr_end1,yr_start2,yr_end2)

if datasource == "CRU":
  if method == "regression":
    output, _,_, loc_eff, reg_eff = calc_irr_impact_regr(datasource,temp_product,response,t_res,0.0,False,yr_start1,yr_end1,yr_start2,yr_end2)
    
    t1 = output[(f_irr[1,:,:] > 0.1) & (~np.isnan(output)) & (~np.isnan(loc_eff))]
    t2 = output[(f_irr[1,:,:] > 0.2) & (~np.isnan(output)) & (~np.isnan(loc_eff))]
    t3 = output[(f_irr[1,:,:] > 0.3) & (~np.isnan(output)) & (~np.isnan(loc_eff))]
    t4 = output[(f_irr[1,:,:] > 0.4) & (~np.isnan(output)) & (~np.isnan(loc_eff))]
    t5 = output[(f_irr[1,:,:] > 0.5) & (~np.isnan(output)) & (~np.isnan(loc_eff))]
    t6 = output[(f_irr[1,:,:] > 0.6) & (~np.isnan(output)) & (~np.isnan(loc_eff))]
    t7 = output[(f_irr[1,:,:] > 0.7) & (~np.isnan(output)) & (~np.isnan(loc_eff))]
    t8 = output[(f_irr[1,:,:] > 0.8) & (~np.isnan(output)) & (~np.isnan(loc_eff))]
    
    loc_eff1 = loc_eff[(f_irr[1,:,:] > 0.1) & (~np.isnan(output)) & (~np.isnan(loc_eff))]
    loc_eff2 = loc_eff[(f_irr[1,:,:] > 0.2) & (~np.isnan(output)) & (~np.isnan(loc_eff))]
    loc_eff3 = loc_eff[(f_irr[1,:,:] > 0.3) & (~np.isnan(output)) & (~np.isnan(loc_eff))]
    loc_eff4 = loc_eff[(f_irr[1,:,:] > 0.4) & (~np.isnan(output)) & (~np.isnan(loc_eff))]
    loc_eff5 = loc_eff[(f_irr[1,:,:] > 0.5) & (~np.isnan(output)) & (~np.isnan(loc_eff))]
    loc_eff6 = loc_eff[(f_irr[1,:,:] > 0.6) & (~np.isnan(output)) & (~np.isnan(loc_eff))]
    loc_eff7 = loc_eff[(f_irr[1,:,:] > 0.7) & (~np.isnan(output)) & (~np.isnan(loc_eff))]
    loc_eff8 = loc_eff[(f_irr[1,:,:] > 0.8) & (~np.isnan(output)) & (~np.isnan(loc_eff))]
    
  elif method == "threshold":
    t1, CtSt = calc_irr_impact_thres(datasource,temp_product,response,t_res,0.1,yr_start1,yr_end1,yr_start2,yr_end2)
    t2, CtSt = calc_irr_impact_thres(datasource,temp_product,response,t_res,0.2,yr_start1,yr_end1,yr_start2,yr_end2)
    t3, CtSt = calc_irr_impact_thres(datasource,temp_product,response,t_res,0.3,yr_start1,yr_end1,yr_start2,yr_end2)
    t4, CtSt = calc_irr_impact_thres(datasource,temp_product,response,t_res,0.4,yr_start1,yr_end1,yr_start2,yr_end2)
    t5, CtSt = calc_irr_impact_thres(datasource,temp_product,response,t_res,0.5,yr_start1,yr_end1,yr_start2,yr_end2)
    t6, CtSt = calc_irr_impact_thres(datasource,temp_product,response,t_res,0.6,yr_start1,yr_end1,yr_start2,yr_end2)
    t7, CtSt = calc_irr_impact_thres(datasource,temp_product,response,t_res,0.7,yr_start1,yr_end1,yr_start2,yr_end2)
    t8, CtSt = calc_irr_impact_thres(datasource,temp_product,response,t_res,0.8,yr_start1,yr_end1,yr_start2,yr_end2)
    del CtSt
elif datasource == "HadEX2":
  if method == "regression":
    t1, ind_vals = calc_irr_impact_regr(datasource,temp_product,response,t_res,0.05,False,yr_start1,yr_end1,yr_start2,yr_end2)
    t2, ind_vals = calc_irr_impact_regr(datasource,temp_product,response,t_res,0.1,False,yr_start1,yr_end1,yr_start2,yr_end2)
    t3, ind_vals = calc_irr_impact_regr(datasource,temp_product,response,t_res,0.15,False,yr_start1,yr_end1,yr_start2,yr_end2)
    t4, ind_vals = calc_irr_impact_regr(datasource,temp_product,response,t_res,0.2,False,yr_start1,yr_end1,yr_start2,yr_end2)
    t5, ind_vals = calc_irr_impact_regr(datasource,temp_product,response,t_res,0.25,False,yr_start1,yr_end1,yr_start2,yr_end2)
    t6, ind_vals = calc_irr_impact_regr(datasource,temp_product,response,t_res,0.3,False,yr_start1,yr_end1,yr_start2,yr_end2)
    t7, ind_vals = calc_irr_impact_regr(datasource,temp_product,response,t_res,0.35,False,yr_start1,yr_end1,yr_start2,yr_end2)
    t8, ind_vals = calc_irr_impact_regr(datasource,temp_product,response,t_res,0.4,False,yr_start1,yr_end1,yr_start2,yr_end2)
  elif method == "threshold":
    t1, CtSt = calc_irr_impact_thres(datasource,temp_product,response,t_res,0.05,yr_start1,yr_end1,yr_start2,yr_end2)
    t2, CtSt = calc_irr_impact_thres(datasource,temp_product,response,t_res,0.1,yr_start1,yr_end1,yr_start2,yr_end2)
    t3, CtSt = calc_irr_impact_thres(datasource,temp_product,response,t_res,0.15,yr_start1,yr_end1,yr_start2,yr_end2)
    t4, CtSt = calc_irr_impact_thres(datasource,temp_product,response,t_res,0.2,yr_start1,yr_end1,yr_start2,yr_end2)
    t5, CtSt = calc_irr_impact_thres(datasource,temp_product,response,t_res,0.25,yr_start1,yr_end1,yr_start2,yr_end2)
    t6, CtSt = calc_irr_impact_thres(datasource,temp_product,response,t_res,0.3,yr_start1,yr_end1,yr_start2,yr_end2)
    t7, CtSt = calc_irr_impact_thres(datasource,temp_product,response,t_res,0.35,yr_start1,yr_end1,yr_start2,yr_end2)
    t8, CtSt = calc_irr_impact_thres(datasource,temp_product,response,t_res,0.4,yr_start1,yr_end1,yr_start2,yr_end2)

dT_irr = [t1,t2,t3,t4,t5,t6,t7,t8]
loc_eff = [loc_eff1,loc_eff2,loc_eff3,loc_eff4,loc_eff5,loc_eff6,loc_eff7,loc_eff8]
reg_eff = [reg_eff1,reg_eff2,reg_eff3,reg_eff4,reg_eff5,reg_eff6,reg_eff7,reg_eff8]
dT_otherforcings = [np.nanmean(x1) - np.nanmean(x2) for (x1,x2) in zip(loc_eff,dT_irr)]
n_vals = [np.count_nonzero(~np.isnan(x)) for x in dT_irr]