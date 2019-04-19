"""SREX_irr_boxplots.py

author: Auke Visser
date: 12.12.2016

This script makes binned boxplots of dT_tot and dT_irr.
First, the upper two panels are filled (dT_tot), after which
the window searching algorithm is applied to CRU and CESM
output to calculate dT_irr and fill the bottom subplots

"""

import netCDF4 as nc
import numpy as np
import scipy
import os
import matplotlib
matplotlib.rcParams['backend'] = "Qt4Agg"
matplotlib.rcParams['legend.numpoints'] = 1
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LogNorm
import matplotlib.lines as mlines
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
temp_product_CRU = "tmx_max"
temp_product_CESM = "TREFHTMX"

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
figformat = 'jpg' ##'pdf'
#################################

SREX_region = 'SAS'

#Define a figure with four subplots and the x-locations of the 
#data points on the plot for CRU and CESM
f, axs = plt.subplots(2,2,sharex='col',sharey='row',figsize=(12,10))
bins_plot_CESM = np.arange(-0.01,0.89,0.1)
bins_plot_CRU = np.arange(0.01,0.91,0.1)

#Calculate df_irr
_,_,f_irr_CRU = calc_irr_diff('CRU_CESM',response,yr_start1,yr_end1,yr_start2,yr_end2)
_,_,f_irr_CESM = calc_irr_diff('CESM',response,yr_start1,yr_end1,yr_start2,yr_end2)

################################################
#
#Calculate dT_tot
#Only the temperature signal is needed for this,
#obtained from extract_T_irr
#
#N.B. This fills only the upper two panels of 
#the plot
#
################################################

T_out_CRU = extract_T_irr('CRU_CESM','tmx_max','PD', t_res, yr_start1, yr_end1, yr_start2, yr_end2)
T_out_CESM = extract_T_irr('CESM', 'TREFHTMX', 'PD', t_res, yr_start1, yr_end1, yr_start2, yr_end2)
dT_tot_CRU = T_out_CRU[1] - T_out_CRU[0]
dT_tot_CESM = T_out_CESM[1] - T_out_CESM[0]

#Define empty lists where the total global signal 
#and the total regional temperature signal are stored
dT_plot_CRU = []
dT_plot_CESM = []
df_plot_CRU = []
df_plot_CESM = []

#Append the total signals for CRU and CESM to the list
df_plot_CRU.append(f_irr_CRU[1].ravel()[(~np.isnan(f_irr_CRU[1].ravel())) & (~np.isnan(dT_tot_CRU.ravel()))])
dT_plot_CRU.append(dT_tot_CRU.ravel()[(~np.isnan(f_irr_CRU[1].ravel())) & (~np.isnan(dT_tot_CRU.ravel()))])
df_plot_CESM.append(f_irr_CESM[1].ravel()[(~np.isnan(f_irr_CESM[1].ravel())) & (~np.isnan(dT_tot_CESM.ravel()))])
dT_plot_CESM.append(dT_tot_CESM.ravel()[(~np.isnan(f_irr_CESM[1].ravel())) & (~np.isnan(dT_tot_CESM.ravel()))])

#Append the regional signal for CRU and CESM to the list
if SREX_region == 'SAS':
    df_plot_CRU.append(np.append(f_irr_CRU[1,101:133,192:220].ravel()[(~np.isnan(f_irr_CRU[1,101:133,192:220].ravel())) & (~np.isnan(dT_tot_CRU[101:133,192:220].ravel()))],f_irr_CRU[1,117:133,220:225].ravel()[(~np.isnan(f_irr_CRU[1,117:133,220:225].ravel())) & (~np.isnan(dT_tot_CRU[117:133,220:225].ravel()))]))
    dT_plot_CRU.append(np.append(dT_tot_CRU[101:133,192:220].ravel()[(~np.isnan(f_irr_CRU[1,101:133,192:220].ravel())) & (~np.isnan(dT_tot_CRU[101:133,192:220].ravel()))],dT_tot_CRU[117:133,220:225].ravel()[(~np.isnan(f_irr_CRU[1,117:133,220:225].ravel())) & (~np.isnan(dT_tot_CRU[117:133,220:225].ravel()))]))
    df_plot_CESM.append(np.append(f_irr_CESM[1,101:133,192:220].ravel()[(~np.isnan(f_irr_CESM[1,101:133,192:220].ravel())) & (~np.isnan(dT_tot_CESM[101:133,192:220].ravel()))],f_irr_CESM[1,117:133,220:225].ravel()[(~np.isnan(f_irr_CESM[1,117:133,220:225].ravel())) & (~np.isnan(dT_tot_CESM[117:133,220:225].ravel()))]))
    dT_plot_CESM.append(np.append(dT_tot_CESM[101:133,192:220].ravel()[(~np.isnan(f_irr_CESM[1,101:133,192:220].ravel())) & (~np.isnan(dT_tot_CESM[101:133,192:220].ravel()))],dT_tot_CESM[117:133,220:225].ravel()[(~np.isnan(f_irr_CESM[1,117:133,220:225].ravel())) & (~np.isnan(dT_tot_CESM[117:133,220:225].ravel()))]))

#Calculate the medians and Q25-Q75 and plot the results
bins = np.arange(0,1.,0.1)
for j in range(0,2):
    data_dT_CRU = dT_plot_CRU[j]
    data_dT_CESM = dT_plot_CESM[j]
    data_df_CRU = df_plot_CRU[j]
    data_df_CESM = df_plot_CESM[j]
    
    medians = []
    pctls = []
    counts = []
    
    #If there are enough data points in the bin, calculate the bin size, median and quartiles
    for x in range(0,9):
        if data_dT_CRU[(data_df_CRU >= bins[x]) & (data_df_CRU < bins[x+1])].shape[0] >= 5:
            counts.append(data_dT_CRU[(data_df_CRU >= bins[x]) & (data_df_CRU < bins[x+1])].shape[0])
            medians.append(np.median(data_dT_CRU[(data_df_CRU >= bins[x]) & (data_df_CRU < bins[x+1])]))
            pctls.append(np.percentile(data_dT_CRU[(data_df_CRU >= bins[x]) & (data_df_CRU < bins[x+1])],[25,75]))
        else:
            counts.append(data_dT_CRU[(data_df_CRU >= bins[x]) & (data_df_CRU < bins[x+1])].shape[0])
            medians.append(np.nan)
            pctls.append(np.array([np.nan,np.nan]))
    
    #Save Q25, median and Q75 in the lists
    #these have concise notations for a for-loop in python for lists, called a list comprehension
    counts = [x if x>=5 else np.nan for x in counts] 
    pctls_l = [x[0] for x in pctls] #Q25
    pctls_h = [x[1] for x in pctls] #Q75
    
    #Plot the medians...
    axs[0,j].plot(bins_plot_CRU,medians,'ks',markeredgewidth=0.0,markersize=6)
    
    #... and quartiles (as a vertical line)
    for u in range(0,9):
        axs[0,j].plot((bins_plot_CRU[u],bins_plot_CRU[u]), (pctls_l[u],pctls_h[u]), 'k--',linewidth=2.5)
    
    #Define a secondary axis (log scale) and plot the histograms
    ax2 = axs[0,j].twinx()
    ax2.bar(bins_plot_CRU-0.01,counts,0.02,color='grey',alpha=0.5)
    ax2.set_yscale('log')
    ax2.set_ylim([1,1e10])
    if j == 0:
        ax2.set_yticks([])
    else:
        ax2.set_yticks([1e0,1e1,1e3,1e5])
        ax2.set_ylabel('Cell count [-]')
        ax2.yaxis.set_label_coords(1.08, 0.25)
    
    #Repeat the plotting procedure for CESM
    
    #Define empty lists
    medians = []
    pctls = []
    counts = []
    
    #If there are enough data points in the bin, calculate the bin size, median and quartiles
    for x in range(0,9):
        if data_dT_CESM[(data_df_CESM >= bins[x]) & (data_df_CESM < bins[x+1])].shape[0] >= 5:
            counts.append(data_dT_CESM[(data_df_CESM >= bins[x]) & (data_df_CESM < bins[x+1])].shape[0])
            medians.append(np.median(data_dT_CESM[(data_df_CESM >= bins[x]) & (data_df_CESM < bins[x+1])]))
            pctls.append(np.percentile(data_dT_CESM[(data_df_CESM >= bins[x]) & (data_df_CESM < bins[x+1])],[25,75]))
        else:
            counts.append(data_dT_CESM[(data_df_CESM >= bins[x]) & (data_df_CESM < bins[x+1])].shape[0])
            medians.append(np.nan)
            pctls.append(np.array([np.nan,np.nan]))

    counts = [x if x>=5 else np.nan for x in counts]
    pctls_l = [x[0] for x in pctls]
    pctls_h = [x[1] for x in pctls]
    
    #Plot the medians...
    axs[0,j].plot(bins_plot_CESM,medians,'ks',markeredgewidth=0.0,markersize=6)
    
    #... and the quartiles as a vertical line
    for u in range(0,9):
        axs[0,j].plot((bins_plot_CESM[u],bins_plot_CESM[u]), (pctls_l[u],pctls_h[u]),'k-',linewidth=2.5)
    
    #Plot the CESM cell count histogram on the secondary axis
    ax2.bar(bins_plot_CESM-0.01,counts,0.02,color='white',alpha=0.5)

###################################
#
#Calculate dT_irr
#Execute the algorithm to calculate
#dT_irr globally and regionally.
#
#N.B. this procedure fills the 
#bottom two panels of the plot.
#
###################################

#Define stacked array in which the results are stored
CRU_t = np.zeros((9, 192, 288))
CESM_t = np.zeros((9, 192, 288))

#Execute the algorithm for different threshold values for CRU and CESM
CRU_t[0,:,:], _,_,_,_ = calc_irr_impact_regr('CRU_CESM',temp_product_CRU,response,t_res,0.0,False,yr_start1,yr_end1,yr_start2,yr_end2)
CRU_t[1,:,:], _,_,_,_ = calc_irr_impact_regr('CRU_CESM',temp_product_CRU,response,t_res,0.1,False,yr_start1,yr_end1,yr_start2,yr_end2)
CRU_t[2,:,:], _,_,_,_ = calc_irr_impact_regr('CRU_CESM',temp_product_CRU,response,t_res,0.2,False,yr_start1,yr_end1,yr_start2,yr_end2)
CRU_t[3,:,:], _,_,_,_ = calc_irr_impact_regr('CRU_CESM',temp_product_CRU,response,t_res,0.3,False,yr_start1,yr_end1,yr_start2,yr_end2)
CRU_t[4,:,:], _,_,_,_ = calc_irr_impact_regr('CRU_CESM',temp_product_CRU,response,t_res,0.4,False,yr_start1,yr_end1,yr_start2,yr_end2)
CRU_t[5,:,:], _,_,_,_ = calc_irr_impact_regr('CRU_CESM',temp_product_CRU,response,t_res,0.5,False,yr_start1,yr_end1,yr_start2,yr_end2)
CRU_t[6,:,:], _,_,_,_ = calc_irr_impact_regr('CRU_CESM',temp_product_CRU,response,t_res,0.6,False,yr_start1,yr_end1,yr_start2,yr_end2)
CRU_t[7,:,:], _,_,_,_ = calc_irr_impact_regr('CRU_CESM',temp_product_CRU,response,t_res,0.7,False,yr_start1,yr_end1,yr_start2,yr_end2)
CRU_t[8,:,:], _,_,_,_ = calc_irr_impact_regr('CRU_CESM',temp_product_CRU,response,t_res,0.8,False,yr_start1,yr_end1,yr_start2,yr_end2)

CESM_t[0,:,:], _,_,_,_ = calc_irr_impact_regr('CESM',temp_product_CESM,response,t_res,0.0,False,yr_start1,yr_end1,yr_start2,yr_end2)
CESM_t[1,:,:], _,_,_,_ = calc_irr_impact_regr('CESM',temp_product_CESM,response,t_res,0.1,False,yr_start1,yr_end1,yr_start2,yr_end2)
CESM_t[2,:,:], _,_,_,_ = calc_irr_impact_regr('CESM',temp_product_CESM,response,t_res,0.2,False,yr_start1,yr_end1,yr_start2,yr_end2)
CESM_t[3,:,:], _,_,_,_ = calc_irr_impact_regr('CESM',temp_product_CESM,response,t_res,0.3,False,yr_start1,yr_end1,yr_start2,yr_end2)
CESM_t[4,:,:], _,_,_,_ = calc_irr_impact_regr('CESM',temp_product_CESM,response,t_res,0.4,False,yr_start1,yr_end1,yr_start2,yr_end2)
CESM_t[5,:,:], _,_,_,_ = calc_irr_impact_regr('CESM',temp_product_CESM,response,t_res,0.5,False,yr_start1,yr_end1,yr_start2,yr_end2)
CESM_t[6,:,:], _,_,_,_ = calc_irr_impact_regr('CESM',temp_product_CESM,response,t_res,0.6,False,yr_start1,yr_end1,yr_start2,yr_end2)
CESM_t[7,:,:], _,_,_,_ = calc_irr_impact_regr('CESM',temp_product_CESM,response,t_res,0.7,False,yr_start1,yr_end1,yr_start2,yr_end2)
CESM_t[8,:,:], _,_,_,_ = calc_irr_impact_regr('CESM',temp_product_CESM,response,t_res,0.8,False,yr_start1,yr_end1,yr_start2,yr_end2)

#Change to 1D, remove NaNs
CRU0 = CRU_t[0,:,:].ravel()[~np.isnan(CRU_t[0,:,:].ravel())]
CRU1 = CRU_t[1,:,:].ravel()[~np.isnan(CRU_t[1,:,:].ravel())]
CRU2 = CRU_t[2,:,:].ravel()[~np.isnan(CRU_t[2,:,:].ravel())]
CRU3 = CRU_t[3,:,:].ravel()[~np.isnan(CRU_t[3,:,:].ravel())]
CRU4 = CRU_t[4,:,:].ravel()[~np.isnan(CRU_t[4,:,:].ravel())]
CRU5 = CRU_t[5,:,:].ravel()[~np.isnan(CRU_t[5,:,:].ravel())]
CRU6 = CRU_t[6,:,:].ravel()[~np.isnan(CRU_t[6,:,:].ravel())]
CRU7 = CRU_t[7,:,:].ravel()[~np.isnan(CRU_t[7,:,:].ravel())]
CRU8 = CRU_t[8,:,:].ravel()[~np.isnan(CRU_t[8,:,:].ravel())]

CESM0 = CESM_t[0,:,:].ravel()[~np.isnan(CESM_t[0,:,:].ravel())]
CESM1 = CESM_t[1,:,:].ravel()[~np.isnan(CESM_t[1,:,:].ravel())]
CESM2 = CESM_t[2,:,:].ravel()[~np.isnan(CESM_t[2,:,:].ravel())]
CESM3 = CESM_t[3,:,:].ravel()[~np.isnan(CESM_t[3,:,:].ravel())]
CESM4 = CESM_t[4,:,:].ravel()[~np.isnan(CESM_t[4,:,:].ravel())]
CESM5 = CESM_t[5,:,:].ravel()[~np.isnan(CESM_t[5,:,:].ravel())]
CESM6 = CESM_t[6,:,:].ravel()[~np.isnan(CESM_t[6,:,:].ravel())]
CESM7 = CESM_t[7,:,:].ravel()[~np.isnan(CESM_t[7,:,:].ravel())]
CESM8 = CESM_t[8,:,:].ravel()[~np.isnan(CESM_t[8,:,:].ravel())]

#Plotting results

########################################
#Fill bottom left subfigure (dT_irr GIL)
########################################

#Change into format readable by boxplot, removing entries without data if needed
CRU_data = [CRU0, CRU1, CRU2, CRU3, CRU4, CRU5, CRU6, CRU7, CRU8]
CRU_dataplot = [x for x in CRU_data if x.shape[0] > 0]
CRU_lenlist = [x.shape[0] for x in CRU_dataplot]

CESM_data = [CESM0, CESM1, CESM2, CESM3, CESM4, CESM5, CESM6, CESM7, CESM8]
CESM_dataplot = [x[(x > -999) & (x < 999)] for x in CESM_data if x.shape[0] > 0]
CESM_lenlist = [x.shape[0] for x in CESM_dataplot]

#Calculate medians and percentiles
CRU_medians = [np.median(x) for x in CRU_dataplot]
CESM_medians = [np.median(x) for x in CESM_dataplot]

CRU_pctls = [np.percentile(x,[25,75]) for x in CRU_dataplot]
CESM_pctls = [np.percentile(x,[25,75]) for x in CESM_dataplot]

#Only plot the results if there are enough data points for the f_irr bin
for u in range(0,9):
    if CRU_lenlist[u] < 5:
        CRU_medians[u] = np.nan
        CRU_pctls[u] = np.array([np.nan,np.nan])
    
    if CESM_lenlist[u] < 5:
        CESM_medians[u] = np.nan
        CESM_pctls[u] = np.array([np.nan,np.nan])
    
    #Plot the IQR of dT_irr as a vertical line
    axs[1,0].plot((bins_plot_CRU[u],bins_plot_CRU[u]), (CRU_pctls[u][0],CRU_pctls[u][1]), 'k--',linewidth=2.5)
    axs[1,0].plot((bins_plot_CESM[u],bins_plot_CESM[u]), (CESM_pctls[u][0],CESM_pctls[u][1]),'-',color='k',linewidth=2.5)
    
#Plot the dT_irr medians
axs[1,0].plot(bins_plot_CRU,CRU_medians,'ks',markeredgewidth=0.0,markersize=6)
axs[1,0].plot(bins_plot_CESM,CESM_medians,'ks',markeredgewidth=0.0,markersize=6)

##############################################
#Fill bottom right subfigure (dT_irr regional)
##############################################

#Extract the SAS signal for each threshold for CRU and CESM, removing NaN values 
CRU0_reg = np.append(CRU_t[0,101:133,192:220].ravel()[~np.isnan(CRU_t[0,101:133,192:220].ravel())],CRU_t[0,117:133,220:225].ravel()[~np.isnan(CRU_t[0,117:133,220:225].ravel())])
CRU1_reg = np.append(CRU_t[1,101:133,192:220].ravel()[~np.isnan(CRU_t[1,101:133,192:220].ravel())],CRU_t[1,117:133,220:225].ravel()[~np.isnan(CRU_t[1,117:133,220:225].ravel())])
CRU2_reg = np.append(CRU_t[2,101:133,192:220].ravel()[~np.isnan(CRU_t[2,101:133,192:220].ravel())],CRU_t[2,117:133,220:225].ravel()[~np.isnan(CRU_t[2,117:133,220:225].ravel())])
CRU3_reg = np.append(CRU_t[3,101:133,192:220].ravel()[~np.isnan(CRU_t[3,101:133,192:220].ravel())],CRU_t[3,117:133,220:225].ravel()[~np.isnan(CRU_t[3,117:133,220:225].ravel())])
CRU4_reg = np.append(CRU_t[4,101:133,192:220].ravel()[~np.isnan(CRU_t[4,101:133,192:220].ravel())],CRU_t[4,117:133,220:225].ravel()[~np.isnan(CRU_t[4,117:133,220:225].ravel())])
CRU5_reg = np.append(CRU_t[5,101:133,192:220].ravel()[~np.isnan(CRU_t[5,101:133,192:220].ravel())],CRU_t[5,117:133,220:225].ravel()[~np.isnan(CRU_t[5,117:133,220:225].ravel())])
CRU6_reg = np.append(CRU_t[6,101:133,192:220].ravel()[~np.isnan(CRU_t[6,101:133,192:220].ravel())],CRU_t[6,117:133,220:225].ravel()[~np.isnan(CRU_t[6,117:133,220:225].ravel())])
CRU7_reg = np.append(CRU_t[7,101:133,192:220].ravel()[~np.isnan(CRU_t[7,101:133,192:220].ravel())],CRU_t[7,117:133,220:225].ravel()[~np.isnan(CRU_t[7,117:133,220:225].ravel())])
CRU8_reg = np.append(CRU_t[8,101:133,192:220].ravel()[~np.isnan(CRU_t[8,101:133,192:220].ravel())],CRU_t[8,117:133,220:225].ravel()[~np.isnan(CRU_t[8,117:133,220:225].ravel())])

CESM0_reg = np.append(CESM_t[0,101:133,192:220].ravel()[~np.isnan(CESM_t[0,101:133,192:220].ravel())],CESM_t[0,117:133,220:225].ravel()[~np.isnan(CESM_t[0,117:133,220:225].ravel())])
CESM1_reg = np.append(CESM_t[1,101:133,192:220].ravel()[~np.isnan(CESM_t[1,101:133,192:220].ravel())],CESM_t[1,117:133,220:225].ravel()[~np.isnan(CESM_t[1,117:133,220:225].ravel())])
CESM2_reg = np.append(CESM_t[2,101:133,192:220].ravel()[~np.isnan(CESM_t[2,101:133,192:220].ravel())],CESM_t[2,117:133,220:225].ravel()[~np.isnan(CESM_t[2,117:133,220:225].ravel())])
CESM3_reg = np.append(CESM_t[3,101:133,192:220].ravel()[~np.isnan(CESM_t[3,101:133,192:220].ravel())],CESM_t[3,117:133,220:225].ravel()[~np.isnan(CESM_t[3,117:133,220:225].ravel())])
CESM4_reg = np.append(CESM_t[4,101:133,192:220].ravel()[~np.isnan(CESM_t[4,101:133,192:220].ravel())],CESM_t[4,117:133,220:225].ravel()[~np.isnan(CESM_t[4,117:133,220:225].ravel())])
CESM5_reg = np.append(CESM_t[5,101:133,192:220].ravel()[~np.isnan(CESM_t[5,101:133,192:220].ravel())],CESM_t[5,117:133,220:225].ravel()[~np.isnan(CESM_t[5,117:133,220:225].ravel())])
CESM6_reg = np.append(CESM_t[6,101:133,192:220].ravel()[~np.isnan(CESM_t[6,101:133,192:220].ravel())],CESM_t[6,117:133,220:225].ravel()[~np.isnan(CESM_t[6,117:133,220:225].ravel())])
CESM7_reg = np.append(CESM_t[7,101:133,192:220].ravel()[~np.isnan(CESM_t[7,101:133,192:220].ravel())],CESM_t[7,117:133,220:225].ravel()[~np.isnan(CESM_t[7,117:133,220:225].ravel())])
CESM8_reg = np.append(CESM_t[8,101:133,192:220].ravel()[~np.isnan(CESM_t[8,101:133,192:220].ravel())],CESM_t[8,117:133,220:225].ravel()[~np.isnan(CESM_t[8,117:133,220:225].ravel())])

#Store the results for SAS into a format readable by the boxplot function
CRU_data_reg = [CRU0_reg, CRU1_reg, CRU2_reg, CRU3_reg, CRU4_reg, CRU5_reg, CRU6_reg, CRU7_reg, CRU8_reg]
CRU_dataplot_reg = [x for x in CRU_data_reg if x.shape[0] > 0]
CRU_lenlist_reg = [x.shape[0] for x in CRU_dataplot_reg]

CESM_data_reg = [CESM0_reg, CESM1_reg, CESM2_reg, CESM3_reg, CESM4_reg, CESM5_reg, CESM6_reg, CESM7_reg, CESM8_reg]
CESM_dataplot_reg = [x for x in CESM_data_reg if x.shape[0] > 0]
CESM_lenlist_reg = [x.shape[0] for x in CESM_dataplot_reg]

#Calculate the medians and the percentiles
CRU_reg_medians = [np.median(x) for x in CRU_dataplot_reg]
CESM_reg_medians = [np.median(x) for x in CESM_dataplot_reg]

CRU_reg_pctls = [np.percentile(x,[25,75]) for x in CRU_dataplot_reg]
CESM_reg_pctls = [np.percentile(x,[25,75]) for x in CESM_dataplot_reg]

#Change values into NaNs if there is not enough data in the bin
for u in range(0,9):
    if CRU_lenlist_reg[u] < 5:
        CRU_reg_medians[u] = np.nan
        CRU_reg_pctls[u] = np.array([np.nan,np.nan])
    if CESM_lenlist_reg[u] < 5:
        CESM_reg_medians[u] = np.nan
        CESM_reg_pctls[u] = np.array([np.nan,np.nan])
    
    #Plot the IQR for CRU and CESM
    axs[1,1].plot((bins_plot_CRU[u],bins_plot_CRU[u]), (CRU_reg_pctls[u][0],CRU_reg_pctls[u][1]), 'k--',linewidth=2.5)
    axs[1,1].plot((bins_plot_CESM[u],bins_plot_CESM[u]), (CESM_reg_pctls[u][0],CESM_reg_pctls[u][1]),'k-',linewidth=2.5)

#Plot the medians for CRU and CESM
axs[1,1].plot(bins_plot_CRU,CRU_reg_medians,'ks',markeredgewidth=0.0,markersize=6)
axs[1,1].plot(bins_plot_CESM,CESM_reg_medians,'ks',markeredgewidth=0.0,markersize=6)

##############################
#Some general plotting options
##############################

panel_list = ['a)','b)','c)','d)']
if SREX_region == 'SAS':
    title_list = ['Global land', 'South Asia', 'Global irrigated land', 'South Asia irrigated']
for ii in range(0,2):
    for jj in range(0,2):
        if (ii == 0 and jj == 0):
            axs[ii,jj].set_ylabel(r'$\Delta T_{tot}$ [$^\circ$C]',size=15)
        if (ii == 1 and jj == 0):
            axs[ii,jj].set_ylabel(r'$\Delta T_{irr}$ [$^\circ$C]',size=15)
        
        axs[ii,jj].set_xlim([-0.05,0.75])
        axs[ii,jj].set_ylim([-2.5,2.5])
        axs[ii,jj].set_xticks(np.arange(0,0.9,0.1))
        axs[ii,jj].set_yticks(np.arange(-2,3,1))
        axs[ii,jj].set_title(title_list[2*ii+jj])
        axs[ii,jj].fill_between([-2,2],0,3,color='salmon',alpha=0.5)
        axs[ii,jj].fill_between([-2,2],-3,0,color='lightblue',alpha=0.7)
        axs[ii,jj].grid(True)
        axs[ii,jj].annotate(panel_list[2*ii+jj],xy=(-0.035,2.1),fontsize=14,ha='left',va='bottom')

f.text(0.5, 0.02, '$f_{irr,PD}$ [-]', ha='center', va='center',fontsize=15)
f.text(0.5, 0.96, '$t_{ref}$: %i-%i, $t_{PD}$: %i-%i'%(yr_start1,yr_end1,yr_start2,yr_end2),ha='center',va='center',fontsize=14)
plt.subplots_adjust(left=0.05, right=0.95, top=0.9, bottom=0.08)

#Plot legends
CESMline = mlines.Line2D([], [], color='k',marker='s',linewidth=2,markeredgewidth=0.0,label='CESM')
CRUline = mlines.Line2D([], [], color='k',marker='s',linestyle='--',linewidth=2,markeredgewidth=0.0,label='CRU')
axs[0,0].legend((CESMline,CRUline),('CESM','CRU'))

#plt.show(block=False)
    
#Save results  
if figsave == True:
    figpath = 'dT_tot_irr.%s.%s.all.%i-%i.%i-%i.WT_testwith20cirr.%s'%(SREX_region,response,yr_start1,yr_end1,yr_start2,yr_end2,figformat)
    
    # if os.path.exists(figpath):
    #   print('Figure already exists at path: ' + figpath)
    # else:
    #   print('Saving figure as:' + figpath)
    plt.savefig(figpath, bbox_inches='tight')
