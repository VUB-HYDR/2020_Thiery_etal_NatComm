"""dT_irr_SREX.py

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

#Call functions from other files
execfile('calc_irr_impact_thres.py')
execfile('calc_irr_impact_regr.py')

###############################
# User-specified options
###############################

#Choose region, season and time period
SREX_region = 'GIL' #'GIL' (global irrigated land), 'CNA', 'EAS', 'MED', 'SAS', 'SEA', 'WAS', 'WNA'
EOBS_region = 'MED' #'EUR', 'MED' (Only works for E-OBS)
plot_season = 2 # 0 = DJF, 1 = MAM, 2 = JJA, 3 = SON (only works if plots are generated for CRU seasonal)\
plot_month = 5 #month number (only works if spatial plots are generated for CRU monthly)
yr_start1 = 1901
yr_end1 = 1930
yr_start2 = 1981
yr_end2 = 2010

#Select method, data product and response

#NOTE: thres_irr works differently for the threshold-based and regression-based method.
#      For the TB, it functions as the cutoff threshold used in the calculation, and
#      for the RB, it functions as the threshold for the variable to be plotted
method = 'regression' #'threshold', 'regression'
thres_irr = 0.25
u = 0.5
datasource = 'CRU_CESM' #'CRU', 'E-OBS', 'HadEX2', 'CESM', 'CRU_CESM'
temp_product = 'tmx_max' 
#CRU: 		'tmp', 'tmx', 'tmn', 'dtr', 'tmp_max', 'dtr_max', 'tmx_max', 'tmn_min','tmn_max' 
#E-OBS: 	'tx'
#HadEX2:	'Ann' (TXx)
#CESM: 'TREFHTMX'
#CRU_CESM: 'tmx_max' 'tmp_max'
response = 'PD' #'PC/PD', 'PD'
t_res = 'seasonal' #'seasonal', 'monthly' (only works for CRU standard products!)
def_regr = False #Whether to include the deforestation rate as a regressor for the RB method (only works for CRU products!)
scaling = False #Whether to scale the TB results with df_irr in the grid cell of interest

output_format = 'both' #'print', 'spatial plot', 'both', 'boxplot'
p_value = 0.01
figshow = True #Currently only works for monthly boxplots
figsave = False
figformat = 'pdf' #'pdf', 'png'
###############################
# End of user-specified options
###############################

if output_format == 'print':
  figsave = False

season_list = ['DJF', 'MAM', 'JJA', 'SON']
month_list = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov','Dec']

#Function that draws a rectangular shape along the SREX region of interest
def plot_SREX_shape(bmap, region):
  #bmap: Basemap
  #region: SREX region for which the shape will be plotted
  if region == 'CNA':
    lot = [(-105,28.566),(-105,50),(-85,50),(-85,28.566)]
  elif region == 'EAS':
    lot = [(100,20),(100,50),(145,50),(145,20)]
  elif region == 'MED':
    lot = [(-10,30),(-10,45),(40,45),(40,30)]
  elif region == 'SAS':
    lot = [(60,5),(60,35),(100,35),(100,20),(95,20),(95,5)]
  elif region == 'SEA':
    lot = [(95,-10),(95,20),(155,20),(155,-10)]
  elif region == 'WAS':
    lot = [(40,15),(40,50),(60,50),(60,15)]
  elif region == 'WNA':
    lot = [(-130,28.566),(-130,60),(-105,60),(-105,28.566)]
  
  xs = []
  ys = []
  for i in range(0,len(lot)-1):
    xs.append(lot[i][0])
    xs.append(lot[i+1][0])
    ys.append(lot[i][1])
    ys.append(lot[i+1][1])
  xs.append(lot[0][0])
  xs.append(lot[-1][0])
  ys.append(lot[0][1])
  ys.append(lot[-1][1])
  
  bmap.plot(xs,ys,latlon=True,color='red',zorder=5)

#Perform calculations
_,_,f_irr = calc_irr_diff(datasource,response,yr_start1,yr_end1,yr_start2,yr_end2)
df_irr = f_irr[1] - f_irr[0]

if method == 'threshold':
  if scaling == True:
    _,_, dT_irr = calc_irr_impact_thres(datasource,temp_product,response,t_res,thres_irr,u,yr_start1,yr_end1,yr_start2,yr_end2)
  else:
    dT_irr, _, _ = calc_irr_impact_thres(datasource,temp_product,response,t_res,thres_irr,u,yr_start1,yr_end1,yr_start2,yr_end2)
elif method == 'regression':
  dT_irr, _,_,_,_ = calc_irr_impact_regr(datasource,temp_product,response,t_res,thres_irr,def_regr,yr_start1,yr_end1,yr_start2,yr_end2)
  if len(dT_irr.shape) == 2:
    if response == "PC/PD":
      dT_irr[f_irr[1,:,:] - f_irr[0,:,:] < thres_irr] = np.nan
    elif response == "PD":
      dT_irr[f_irr[1,:,:] < thres_irr] = np.nan
  elif len(dT_irr.shape) == 3:
    if response == "PC/PD":
      dT_irr[np.tile(f_irr[1,:,:] - f_irr[0,:,:], (dT_irr.shape[0],1,1)) < thres_irr] = np.nan
    elif response == "PD":
      dT_irr[np.tile(f_irr[1,:,:], (dT_irr.shape[0],1,1)) < thres_irr] = np.nan
      
#Some stuff needed for discrete colorbar
bounds = np.linspace(-0.45,0.45,10)
norm = colors.BoundaryNorm(boundaries=bounds,ncolors=256)

#OUTPUTTING RESULTS ACCORDING TO USER-SPECIFIED OPTIONS
if(datasource == "CRU" and t_res == "monthly" and temp_product in ['tmp', 'dtr', 'tmn', 'tmx']):
  #Generate monthly boxplots to analyse average annual cyclicity
  #First the region is selected, after which the temperature output is stored
  #in the variable 'out' for all months 
  if(output_format == 'boxplot'):
    if SREX_region == 'GIL':
      out = np.zeros((12,dT_irr[0,:,:].ravel()[~np.isnan(dT_irr[0,:,:].ravel())].shape[0]))
      for i in range(0,12):
        out[i,:] = dT_irr[i,:,:].ravel()[~np.isnan(dT_irr[i,:,:].ravel())]
    elif SREX_region == 'CNA':
      out = np.zeros((12,dT_irr[0,237:280,150:190].ravel()[~np.isnan(dT_irr[0,237:280,150:190].ravel())].shape[0]))
      for i in range(0,12):
        out[i,:] = dT_irr[i,237:280,150:190].ravel()[~np.isnan(dT_irr[i,237:280,150:190].ravel())]
    elif SREX_region == 'EAS':
      out = np.zeros((12,dT_irr[0,220:280,560:650].ravel()[~np.isnan(dT_irr[0,220:280,560:650].ravel())].shape[0]))
      for i in range(0,12):
        out[i,:] = dT_irr[i,220:280,560:650].ravel()[~np.isnan(dT_irr[i,220:280,560:650].ravel())]
    elif SREX_region == 'MED':
      out = np.zeros((12,dT_irr[0,240:270,340:440].ravel()[~np.isnan(dT_irr[0,240:270,340:440].ravel())].shape[0]))
      for i in range(0,12):
        out[i,:] = dT_irr[i,240:270,340:440].ravel()[~np.isnan(dT_irr[i,240:270,340:440].ravel())]
    elif SREX_region == 'SAS':
      SAS1_shape = dT_irr[0,190:250,480:550].ravel()[~np.isnan(dT_irr[0,190:250,480:550].ravel())].shape[0]
      SAS2_shape = dT_irr[0,220:250,550:560].ravel()[~np.isnan(dT_irr[0,220:250,550:560].ravel())].shape[0]
      out = np.zeros((12,SAS1_shape+SAS2_shape))
      SAS1 = np.zeros((12,SAS1_shape))
      SAS2 = np.zeros((12,SAS2_shape))
      for i in range(0,12):
        SAS1[i,:] = dT_irr[i,190:250,480:550].ravel()[~np.isnan(dT_irr[i,190:250,480:550].ravel())]
        SAS2[i,:] = dT_irr[i,220:250,550:560].ravel()[~np.isnan(dT_irr[i,220:250,550:560].ravel())]
        out[i,:] = np.append(SAS1[i,:],SAS2[i,:])
    elif SREX_region == 'SEA':
      out = np.zeros((12,dT_irr[0,160:220,550:670].ravel()[~np.isnan(dT_irr[0,160:220,550:670].ravel())].shape[0]))
      for i in range(0,12):
        out[i,:] = dT_irr[i,160:220,550:670].ravel()[~np.isnan(dT_irr[i,160:220,550:670].ravel())]
    elif SREX_region == 'WAS':
      out = np.zeros((12,dT_irr[0,210:280,440:480].ravel()[~np.isnan(dT_irr[0,210:280,440:480].ravel())].shape[0]))
      for i in range(0,12):
        out[i,:] = dT_irr[i,210:280,440:480].ravel()[~np.isnan(dT_irr[i,210:280,440:480].ravel())]
    elif SREX_region == 'WNA':
      out = np.zeros((12,dT_irr[0,237:300,100:150].ravel()[~np.isnan(dT_irr[0,237:300,100:150].ravel())].shape[0]))
      for i in range(0,12):
        out[i,:] = dT_irr[i,237:300,100:150].ravel()[~np.isnan(dT_irr[i,237:300,100:150].ravel())]
    
    #Store the temperature output in a list
    data = [out[0],out[1],out[2],out[3],out[4],out[5],out[6],out[7],out[8],out[9],out[10],out[11]]
    
    #Define a figure and draw the boxplot
    plt.figure(figsize=(12,8))
    bp = plt.boxplot(data,patch_artist=True)
    
    #Perform statistical tests and color the boxplots differently based on the outcome
    for a in range(0,len(data)):
      p_ttest = stats.ttest_1samp(data[a],0)[1]
      p_shapiro = stats.shapiro(data[a])[1]
      if p_ttest <= p_value:
        bp['boxes'][a].set(edgecolor='green',linewidth=3)
      if p_shapiro <= p_value:
        bp['boxes'][a].set(hatch='//')
        bp['boxes'][a].set(facecolor='white')
    
    #Some commands to change the figure layout
    meanpatch = mpatches.Patch(edgecolor='green',facecolor='white',linewidth=3,label=r'reject $H_0: \mu = 0$')
    normalpatch = mpatches.Patch(facecolor='white',linewidth=0,hatch='//',label=r'reject $H_0: X \sim N(\mu,\sigma^2)$')
    plt.xticks([1,2,3,4,5,6,7,8,9,10,11,12],month_list)
    plt.xlabel('Month')
    plt.ylabel('T change [$^\circ$C]')
    plt.ylim([-1,1])
    plt.plot(np.arange(1,13),[x.mean() for x in data],'k--') #generates a line with the mean annual dT_irr signal
    plt.grid(True)
    plt.legend(handles=[meanpatch,normalpatch],loc='upper right')
    plt.title('method: %s, threshold: %1.2f, dataset: %s %s, region: %s, p: %s, \n ref: %i-%i'%(method,thres_irr,datasource,temp_product,SREX_region,p_value,yr_start1,yr_end1))
    if figshow == True:
      plt.show(block=False)
    
    #Save the figure
    if figsave == True:
      if(method == 'regression' and def_regr == True):
        figpath = '/net/exo/landclim/wthiery/observational_analysis/Figures/CRU/bxp_month/%s/dT_irr.bxp_m.%s.def.%1.2f.%s.%s.%s-%s.%s-%s.%s'%(SREX_region,method,thres_irr,temp_product,SREX_region,yr_start1,yr_end1,yr_start2,yr_end2,figformat)
      else:
        figpath = '/net/exo/landclim/wthiery/observational_analysis/Figures/CRU/bxp_month/%s/dT_irr.bxp_m.%s.%1.2f.%s.%s.%s-%s.%s-%s.%s'%(SREX_region,method,thres_irr,temp_product,SREX_region,yr_start1,yr_end1,yr_start2,yr_end2,figformat)
      if os.path.exists(figpath):
        print('Figure already exists at path: ' + figpath)
      else:
        print('Saving figure as:' + figpath)
        plt.savefig(figpath,bbox_inches='tight')
      
elif(datasource == "CRU" and temp_product in ['tmp', 'dtr', 'tmn', 'tmx'] and t_res == "seasonal"):    
  #Produce output based on user-specified options
  if(output_format in ['print', 'both']):
    if SREX_region == 'GIL':
      print np.nanmean(np.nanmean(dT_irr,axis=2),axis=1)
    elif SREX_region == 'CNA':
      print np.nanmean(np.nanmean(dT_irr[:,237:280,150:190],axis=2),axis=1)
    elif SREX_region == 'EAS':
      print np.nanmean(np.nanmean(dT_irr[:,220:280,560:650],axis=2),axis=1)
    elif SREX_region == 'MED':
      print np.nanmean(np.nanmean(dT_irr[:,240:270,340:440],axis=2),axis=1)
    elif SREX_region == 'SAS':
      c1 = np.count_nonzero(df_irr[190:250,480:550] >= thres_irr)
      c2 = np.count_nonzero(df_irr[220:250,550:560] >= thres_irr)
      if c2 == 0:
        print np.nanmean(np.nanmean(dT_irr[:,190:250,480:550],axis=2),axis=1)
      else:
        b1 = np.nanmean(np.nanmean(dT_irr[:,190:250,480:550],axis=2),axis=1)
        b2 = np.nanmean(np.nanmean(dT_irr[:,220:250,550:560],axis=2),axis=1)
        stack = np.concatenate((np.tile(b1,(c1,1)),np.tile(b2,(c2,1))),axis=0)
        print np.nanmean(stack,axis=0)
    elif SREX_region == 'SEA':
      print np.nanmean(np.nanmean(dT_irr[:,160:220,550:670],axis=2),axis=1)
    elif SREX_region == 'WAS':
      print np.nanmean(np.nanmean(dT_irr[:,210:280,440:480],axis=2),axis=1)
    elif SREX_region == 'WNA':
      print np.nanmean(np.nanmean(dT_irr[:,237:300,100:150],axis=2),axis=1)
      
  if(output_format in ['spatial plot', 'both']):
    #The procedure is the same for every region, only the lat/lon bounds and the region to be plotteddiffer
    if SREX_region == 'GIL':
      plt.figure(figsize=(12,8)) #Define a figure
      m = Basemap(projection='cyl',llcrnrlat=-60,urcrnrlat=90,llcrnrlon=-180,urcrnrlon=180,resolution='l') #Define the basemap
      m.drawcoastlines() #Draw the coastlines
      m.drawmapboundary(fill_color='lightgray') #Draw the map background
      m.fillcontinents(color='white',lake_color='lightgray') #Fill the continents
      m.imshow(dT_irr[plot_season,60:360,:],vmin=-0.5,vmax=0.5,norm=norm,cmap='coolwarm',interpolation='none',zorder=4) #Plot dT_irr
      plot_SREX_shape(m,'CNA') #Plot the SREX shapes with the function defined above
      plot_SREX_shape(m,'EAS')
      plot_SREX_shape(m,'MED')
      plot_SREX_shape(m,'SAS')
      plot_SREX_shape(m,'SEA')
      plot_SREX_shape(m,'WAS')
      plot_SREX_shape(m,'WNA')
      c = plt.colorbar(extend='both') #Plot the colorbar
      c.set_label('$\Delta T_{irr}$ [$^\circ$C]')
      plt.title('Region: %s, method: %s, threshold: %1.2f, dataset: %s %s, season: %s, \n ref: %i-%i'%(SREX_region,method,thres_irr,datasource,temp_product,season_list[plot_season],yr_start1,yr_end1))
      plt.show(block=False)
    elif SREX_region == 'CNA':
      plt.figure(figsize=(12,8))
      m = Basemap(projection='cyl',llcrnrlat=18.566,urcrnrlat=60,llcrnrlon=-115,urcrnrlon=-75,resolution='l')
      m.drawcoastlines()
      m.drawmapboundary(fill_color='lightgray')
      m.fillcontinents(color='white',lake_color='lightgray')
      m.imshow(dT_irr[plot_season,217:300,130:210],vmin=-0.5,vmax=0.5,norm=norm,cmap='coolwarm',interpolation='none',zorder=4)
      plot_SREX_shape(m,SREX_region)
      c = plt.colorbar(extend='both')
      c.set_label('$\Delta T_{irr}$ [$^\circ$C]')
      plt.title('Region: %s, method: %s, threshold: %1.2f, dataset: %s %s, season: %s, \n ref: %i-%i'%(SREX_region,method,thres_irr,datasource,temp_product,season_list[plot_season],yr_start1,yr_end1))
      plt.show(block=False)
    elif SREX_region == 'EAS':
      plt.figure(figsize=(12,8))
      m = Basemap(projection='cyl',llcrnrlat=10,urcrnrlat=60,llcrnrlon=90,urcrnrlon=155,resolution='l')
      m.drawcoastlines()
      m.drawmapboundary(fill_color='lightgray')
      m.fillcontinents(color='white',lake_color='lightgray')
      m.imshow(dT_irr[plot_season,200:300,540:670],vmin=-0.5,vmax=0.5,norm=norm,cmap='coolwarm',interpolation='none',zorder=4)
      plot_SREX_shape(m,SREX_region)
      c = plt.colorbar(extend='both')
      c.set_label('$\Delta T_{irr}$ [$^\circ$C]')
      plt.title('Region: %s, method: %s, threshold: %1.2f, dataset: %s %s, season: %s, \n ref: %i-%i'%(SREX_region,method,thres_irr,datasource,temp_product,season_list[plot_season],yr_start1,yr_end1))
      plt.show(block=False)
    elif SREX_region == 'MED':
      plt.figure(figsize=(12,8))
      m = Basemap(projection='cyl',llcrnrlat=20,urcrnrlat=55,llcrnrlon=-20,urcrnrlon=50,resolution='l')
      m.drawcoastlines()
      m.drawmapboundary(fill_color='lightgray')
      m.fillcontinents(color='white',lake_color='lightgray')
      m.imshow(dT_irr[plot_season,220:290,320:460],vmin=-0.5,vmax=0.5,norm=norm,cmap='coolwarm',interpolation='none',zorder=4)
      plot_SREX_shape(m,SREX_region)
      c = plt.colorbar(extend='both')
      c.set_label('$\Delta T_{irr}$ [$^\circ$C]')
      plt.title('Region: %s, method: %s, threshold: %1.2f, dataset: %s %s, season: %s, \n ref: %i-%i'%(SREX_region,method,thres_irr,datasource,temp_product,season_list[plot_season],yr_start1,yr_end1))
      plt.show(block=False)
    elif SREX_region == 'SAS':
      plt.figure(figsize=(12,8))
      m = Basemap(projection='cyl',llcrnrlat=-5,urcrnrlat=45,llcrnrlon=50,urcrnrlon=110,resolution='l')
      m.drawcoastlines()
      m.drawmapboundary(fill_color='lightgray')
      m.fillcontinents(color='white',lake_color='lightgray')
      plot_var = dT_irr[plot_season,170:270,460:580]
      plot_var[0:30,70:80] = np.nan
      m.imshow(plot_var,vmin=-0.5,vmax=0.5,norm=norm,cmap='coolwarm',interpolation='none',zorder=4)
      plot_SREX_shape(m,SREX_region)
      c = plt.colorbar(extend='both')
      c.set_label('$\Delta T_{irr}$ [$^\circ$C]')
      plt.title('Region: %s, method: %s, threshold: %1.2f, dataset: %s %s, season: %s, \n ref: %i-%i'%(SREX_region,method,thres_irr,datasource,temp_product,season_list[plot_season],yr_start1,yr_end1))
      plt.show(block=False)
    elif SREX_region == 'SEA':
      plt.figure(figsize=(12,8))
      m = Basemap(projection='cyl',llcrnrlat=-20,urcrnrlat=30,llcrnrlon=85,urcrnrlon=165,resolution='l')
      m.drawcoastlines()
      m.drawmapboundary(fill_color='lightgray')
      m.fillcontinents(color='white',lake_color='lightgray')
      m.imshow(dT_irr[plot_season,140:240,530:690],vmin=-0.5,vmax=0.5,norm=norm,cmap='coolwarm',interpolation='none',zorder=4)
      plot_SREX_shape(m,SREX_region)
      c = plt.colorbar(extend='both')
      c.set_label('$\Delta T_{irr}$ [$^\circ$C]')
      plt.title('Region: %s, method: %s, threshold: %1.2f, dataset: %s %s, season: %s, \n ref: %i-%i'%(SREX_region,method,thres_irr,datasource,temp_product,season_list[plot_season],yr_start1,yr_end1))
      plt.show(block=False)
    elif SREX_region == 'WAS':
      plt.figure(figsize=(12,8))
      m = Basemap(projection='cyl',llcrnrlat=5,urcrnrlat=60,llcrnrlon=30,urcrnrlon=70,resolution='l')
      m.drawcoastlines()
      m.drawmapboundary(fill_color='lightgray')
      m.fillcontinents(color='white',lake_color='lightgray')
      m.imshow(dT_irr[plot_season,190:300,420:500],vmin=-0.5,vmax=0.5,cmap='coolwarm',interpolation='none',zorder=4)
      plot_SREX_shape(m,SREX_region)
      c = plt.colorbar(extend='both')
      c.set_label('$\Delta T_{irr}$ [$^\circ$C]')
      plt.title('Region: %s, method: %s, threshold: %1.2f, dataset: %s %s, season: %s, \n ref: %i-%i'%(SREX_region,method,thres_irr,datasource,temp_product,season_list[plot_season],yr_start1,yr_end1))
      plt.show(block=False)
    elif SREX_region == 'WNA':
      plt.figure(figsize=(12,8))
      m = Basemap(projection='cyl',llcrnrlat=18.566,urcrnrlat=70,llcrnrlon=-140,urcrnrlon=-95,resolution='l')
      m.drawcoastlines()
      m.drawmapboundary(fill_color='lightgray')
      m.fillcontinents(color='white',lake_color='lightgray')
      m.imshow(dT_irr[plot_season,217:320,80:170],vmin=-0.5,vmax=0.5,norm=norm,cmap='coolwarm',interpolation='none',zorder=4)
      plot_SREX_shape(m,SREX_region)
      c = plt.colorbar(extend='both')
      c.set_label('$\Delta T_{irr}$ [$^\circ$C]')
      plt.title('Region: %s, method: %s, threshold: %1.2f, dataset: %s %s, season: %s, \n ref: %i-%i'%(SREX_region,method,thres_irr,datasource,temp_product,season_list[plot_season],yr_start1,yr_end1))
      plt.show(block=False)

  if output_format == 'boxplot':
    #A boxplot is generated for every region
    GIL = dT_irr[plot_season,:,:].ravel()[~np.isnan(dT_irr[plot_season,:,:].ravel())]
    CNA = dT_irr[plot_season,237:280,150:190].ravel()[~np.isnan(dT_irr[plot_season,237:280,150:190].ravel())]
    EAS = dT_irr[plot_season,220:280,560:650].ravel()[~np.isnan(dT_irr[plot_season,220:280,560:650].ravel())]
    MED = dT_irr[plot_season,240:270,340:440].ravel()[~np.isnan(dT_irr[plot_season,240:270,340:440].ravel())]
    SAS1 = dT_irr[plot_season,190:250,480:550].ravel()[~np.isnan(dT_irr[plot_season,190:250,480:550].ravel())]
    SAS2 = dT_irr[plot_season,220:250,550:560].ravel()[~np.isnan(dT_irr[plot_season,220:250,550:560].ravel())]
    SAS = np.append(SAS1,SAS2)
    SEA = dT_irr[plot_season,160:220,550:670].ravel()[~np.isnan(dT_irr[plot_season,160:220,550:670].ravel())]
    WAS = dT_irr[plot_season,210:280,440:480].ravel()[~np.isnan(dT_irr[plot_season,210:280,440:480].ravel())]
    WNA = dT_irr[plot_season,237:300,100:150].ravel()[~np.isnan(dT_irr[plot_season,237:300,100:150].ravel())]
    
    #Store all the regional data in a list
    data = [GIL,CNA,EAS,MED,SAS,SEA,WAS,WNA]
    plt.figure(figsize=(12,8))
    bp = plt.boxplot(data,patch_artist=True)
    
    #Perform statistical tests
    for a in range(0,len(data)):
      p_ttest = stats.ttest_1samp(data[a],0)[1]
      p_shapiro = stats.shapiro(data[a])[1]
      if p_ttest <= p_value:
        bp['boxes'][a].set(edgecolor='green',linewidth=3)
      if p_shapiro <= p_value:
        bp['boxes'][a].set(hatch='//')
        bp['boxes'][a].set(facecolor='white')
    
    #Some plotting commands
    meanpatch = mpatches.Patch(edgecolor='green',facecolor='white',linewidth=3,label=r'reject $H_0: \mu = 0$')
    normalpatch = mpatches.Patch(facecolor='white',linewidth=0,hatch='//',label=r'reject $H_0: X \sim N(\mu,\sigma^2)$')
    plt.xticks([1,2,3,4,5,6,7,8],['GIL','CNA','EAS','MED','SAS','SEA','WAS','WNA'])
    plt.xlabel('Region')
    plt.ylabel('T change [$^\circ$C]')
    plt.ylim([-2,2])
    
    lenlist = [x.shape[0] for x in data]
    for i in range(0,len(data)):
      plt.text(i+1,-1.9,str(lenlist[i]),ha='center',va='center')
    
    plt.grid(True)
    plt.legend(handles=[meanpatch,normalpatch],loc='upper right')
    plt.title('method: %s, threshold: %1.2f, dataset: %s %s, season: %s, p: %s, \n ref: %i-%i'%(method,thres_irr,datasource,temp_product,season_list[plot_season],p_value,yr_start1,yr_end1))
    plt.show(block=False)
  
  #Save Figure
  if figsave == True:
    if output_format in ['spatial plot', 'both']:
      if(method == 'regression' and def_regr == True):
        figpath = '/net/exo/landclim/wthiery/observational_analysis/Figures/CRU/%s/%s/dT_irr.%s.def.%1.2f.%s.%s.%s.%s-%s.%s-%s.%s'%(temp_product,SREX_region,method,thres_irr,temp_product,SREX_region,season_list[plot_season],yr_start1,yr_end1,yr_start2,yr_end2,figformat)
      else:
        figpath = '/net/exo/landclim/wthiery/observational_analysis/Figures/CRU/%s/%s/dT_irr.%s.alt.%1.2f.%s.%s.%s.%s-%s.%s-%s.%s'%(temp_product,SREX_region,method,thres_irr,temp_product,SREX_region,season_list[plot_season],yr_start1,yr_end1,yr_start2,yr_end2,figformat)
      if os.path.exists(figpath):
        print('Figure already exists at path: ' + figpath)
      else:  
        print('Saving figure as: ' + figpath)
        plt.savefig(figpath,bbox_inches='tight')
    elif output_format == 'boxplot':
      figpath = '/net/exo/landclim/wthiery/observational_analysis/Figures/CRU/boxplots/%s/dT_irr.boxplot.tests.%s.%1.2f.%s.%s.%s-%s.%s-%s.%s'%(temp_product,method,thres_irr,temp_product,season_list[plot_season],yr_start1,yr_end1,yr_start2,yr_end2,figformat)
      if os.path.exists(figpath):
        print('Figure already exists at path: ' + figpath)
      else:
        print('Saving figure as:' + figpath)
        plt.savefig(figpath,bbox_inches='tight')

elif(datasource == "CRU" and temp_product in ['tmp_max', 'dtr_max', 'tmx_max', 'tmn_min','tmn_max']):
  if(output_format in ['print', 'both']):
    if SREX_region == 'GIL':
      print np.nanmean(dT_irr)
    elif SREX_region == 'CNA':
      print np.nanmean(dT_irr[237:280,150:190])
    elif SREX_region == 'EAS':
      print np.nanmean(dT_irr[220:280,560:650])
    elif SREX_region == 'MED':
      print np.nanmean(dT_irr[240:270,340:440])
    elif SREX_region == 'SAS':
      c1 = np.count_nonzero(df_irr[190:250,480:550] >= thres_irr)
      c2 = np.count_nonzero(df_irr[220:250,550:560] >= thres_irr)
      if c2 == 0:
        print np.nanmean(dT_irr[190:250,480:550])
      else:
        b1 = np.nanmean(dT_irr[190:250,480:550])
        b2 = np.nanmean(dT_irr[220:250,550:560])
        stack = np.concatenate((np.tile(b1,(c1,1)),np.tile(b2,(c2,1))),axis=0)
        print np.nanmean(stack,axis=0)
    elif SREX_region == 'SEA':
      print np.nanmean(dT_irr[160:220,550:670])
    elif SREX_region == 'WAS':
      print np.nanmean(dT_irr[210:280,440:480])
    elif SREX_region == 'WNA':
      print np.nanmean(dT_irr[237:300,100:150])
    
  if(output_format in ['spatial plot', 'both']):
    #The procedure is the same for every region, only the lat/lon bounds and the region to be plotted differ
    if SREX_region == 'GIL':
      plt.figure(figsize=(12,8)) #Define figure
      m = Basemap(projection='cyl',llcrnrlat=-60,urcrnrlat=90,llcrnrlon=-180,urcrnrlon=180,resolution='l') #Define global map
      m.drawcoastlines() #Plot the coast lines
      m.drawmapboundary(fill_color='lightgray') #Plot the map background
      m.fillcontinents(color='white',lake_color='lightgray') #Fill the continents with a different color
      m.imshow(dT_irr[60:360,:],vmin=-0.5,vmax=0.5,norm=norm,cmap='coolwarm',interpolation='none',zorder=4) #Plot dT_irr on the map
      plot_SREX_shape(m,'CNA') #Plot the different SREX shapes
      plot_SREX_shape(m,'EAS')
      plot_SREX_shape(m,'MED')
      plot_SREX_shape(m,'SAS')
      plot_SREX_shape(m,'SEA')
      plot_SREX_shape(m,'WAS')
      plot_SREX_shape(m,'WNA')
      c = plt.colorbar(extend='both') #Plot the colorbar
      c.set_label('$\Delta T_{irr}$ [$^\circ$C]') 
      plt.title('Region: %s, method: %s, threshold: %1.2f, dataset: %s %s, \n ref: %i-%i'%(SREX_region,method,thres_irr,datasource,temp_product,yr_start1,yr_end1))
      plt.show(block=False)
    elif SREX_region == 'CNA':
      plt.figure(figsize=(12,8))
      m = Basemap(projection='cyl',llcrnrlat=18.566,urcrnrlat=60,llcrnrlon=-115,urcrnrlon=-75,resolution='l')
      m.drawcoastlines()
      m.drawmapboundary(fill_color='lightgray')
      m.fillcontinents(color='white',lake_color='lightgray')
      m.imshow(dT_irr[217:300,130:210],vmin=-0.5,vmax=0.5,norm=norm,cmap='coolwarm',interpolation='none',zorder=4)
      plot_SREX_shape(m,SREX_region)
      c = plt.colorbar(extend='both')
      c.set_label('$\Delta T_{irr}$ [$^\circ$C]')
      plt.title('Region: %s, method: %s, threshold: %1.2f, dataset: %s %s, \n ref: %i-%i'%(SREX_region,method,thres_irr,datasource,temp_product,yr_start1,yr_end1))
      plt.show(block=False)
    elif SREX_region == 'EAS':
      plt.figure(figsize=(12,8))
      m = Basemap(projection='cyl',llcrnrlat=10,urcrnrlat=60,llcrnrlon=90,urcrnrlon=155,resolution='l')
      m.drawcoastlines()
      m.drawmapboundary(fill_color='lightgray')
      m.fillcontinents(color='white',lake_color='lightgray')
      m.imshow(dT_irr[200:300,540:670],vmin=-0.5,vmax=0.5,norm=norm,cmap='coolwarm',interpolation='none',zorder=4)
      plot_SREX_shape(m,SREX_region)
      c = plt.colorbar(extend='both')
      c.set_label('$\Delta T_{irr}$ [$^\circ$C]')
      plt.title('Region: %s, method: %s, threshold: %1.2f, dataset: %s %s, \n ref: %i-%i'%(SREX_region,method,thres_irr,datasource,temp_product,yr_start1,yr_end1))
      plt.show(block=False)
    elif SREX_region == 'MED':
      plt.figure(figsize=(12,8))
      m = Basemap(projection='cyl',llcrnrlat=20,urcrnrlat=55,llcrnrlon=-20,urcrnrlon=50,resolution='l')
      m.drawcoastlines()
      m.drawmapboundary(fill_color='lightgray')
      m.fillcontinents(color='white',lake_color='lightgray')
      m.imshow(dT_irr[220:290,320:460],vmin=-0.5,vmax=0.5,norm=norm,cmap='coolwarm',interpolation='none',zorder=4)
      plot_SREX_shape(m,SREX_region)
      c = plt.colorbar(extend='both')
      c.set_label('$\Delta T_{irr}$ [$^\circ$C]')
      plt.title('Region: %s, method: %s, threshold: %1.2f, dataset: %s %s, \n ref: %i-%i'%(SREX_region,method,thres_irr,datasource,temp_product,yr_start1,yr_end1))
      plt.show(block=False)
    elif SREX_region == 'SAS':
      plt.figure(figsize=(12,8))
      m = Basemap(projection='cyl',llcrnrlat=-5,urcrnrlat=45,llcrnrlon=50,urcrnrlon=110,resolution='l')
      m.drawcoastlines()
      m.drawmapboundary(fill_color='lightgray')
      m.fillcontinents(color='white',lake_color='lightgray')
      plot_var = dT_irr[170:270,460:580]
      plot_var[0:30,70:80] = np.nan
      m.imshow(plot_var,vmin=-0.5,vmax=0.5,norm=norm,cmap='coolwarm',interpolation='none',zorder=4)
      plot_SREX_shape(m,SREX_region)
      c = plt.colorbar(extend='both')
      c.set_label('$\Delta T_{irr}$ [$^\circ$C]')
      plt.title('Region: %s, method: %s, threshold: %1.2f, dataset: %s %s, \n ref: %i-%i'%(SREX_region,method,thres_irr,datasource,temp_product,yr_start1,yr_end1))
      plt.show(block=False)
    elif SREX_region == 'SEA':
      plt.figure(figsize=(12,8))
      m = Basemap(projection='cyl',llcrnrlat=-20,urcrnrlat=30,llcrnrlon=85,urcrnrlon=165,resolution='l')
      m.drawcoastlines()
      m.drawmapboundary(fill_color='lightgray')
      m.fillcontinents(color='white',lake_color='lightgray')
      m.imshow(dT_irr[140:240,530:690],vmin=-0.5,vmax=0.5,norm=norm,cmap='coolwarm',interpolation='none',zorder=4)
      plot_SREX_shape(m,SREX_region)
      c = plt.colorbar(extend='both')
      c.set_label('$\Delta T_{irr}$ [$^\circ$C]')
      plt.title('Region: %s, method: %s, threshold: %1.2f, dataset: %s %s, \n ref: %i-%i'%(SREX_region,method,thres_irr,datasource,temp_product,yr_start1,yr_end1))
      plt.show(block=False)
    elif SREX_region == 'WAS':
      plt.figure(figsize=(12,8))
      m = Basemap(projection='cyl',llcrnrlat=5,urcrnrlat=60,llcrnrlon=30,urcrnrlon=70,resolution='l')
      m.drawcoastlines()
      m.drawmapboundary(fill_color='lightgray')
      m.fillcontinents(color='white',lake_color='lightgray')
      m.imshow(dT_irr[190:300,420:500],vmin=-0.5,vmax=0.5,norm=norm,cmap='coolwarm',interpolation='none',zorder=4)
      plot_SREX_shape(m,SREX_region)
      c = plt.colorbar(extend='both')
      c.set_label('$\Delta T_{irr}$ [$^\circ$C]')
      plt.title('Region: %s, method: %s, threshold: %1.2f, dataset: %s %s, \n ref: %i-%i'%(SREX_region,method,thres_irr,datasource,temp_product,yr_start1,yr_end1))
      plt.show(block=False)
    elif SREX_region == 'WNA':
      plt.figure(figsize=(12,8))
      m = Basemap(projection='cyl',llcrnrlat=18.566,urcrnrlat=70,llcrnrlon=-140,urcrnrlon=-95,resolution='l')
      m.drawcoastlines()
      m.drawmapboundary(fill_color='lightgray')
      m.fillcontinents(color='white',lake_color='lightgray')
      m.imshow(dT_irr[217:320,80:170],vmin=-0.5,vmax=0.5,norm=norm,cmap='coolwarm',interpolation='none',zorder=4)
      plot_SREX_shape(m,SREX_region)
      c = plt.colorbar(extend='both')
      c.set_label('$\Delta T_{irr}$ [$^\circ$C]')
      plt.title('Region: %s, method: %s, threshold: %1.2f, dataset: %s %s, \n ref: %i-%i'%(SREX_region,method,thres_irr,datasource,temp_product,yr_start1,yr_end1))
      plt.show(block=False)
      
  if output_format == 'boxplot':
    #This output option generates boxplots showing the result per region
    GIL = dT_irr.ravel()[~np.isnan(dT_irr.ravel())]
    CNA = dT_irr[237:280,150:190].ravel()[~np.isnan(dT_irr[237:280,150:190].ravel())]
    EAS = dT_irr[220:280,560:650].ravel()[~np.isnan(dT_irr[220:280,560:650].ravel())]
    MED = dT_irr[240:270,340:440].ravel()[~np.isnan(dT_irr[240:270,340:440].ravel())]
    SAS1 = dT_irr[190:250,480:550].ravel()[~np.isnan(dT_irr[190:250,480:550].ravel())]
    SAS2 = dT_irr[220:250,550:560].ravel()[~np.isnan(dT_irr[220:250,550:560].ravel())]
    SAS = np.append(SAS1,SAS2)
    SEA = dT_irr[160:220,550:670].ravel()[~np.isnan(dT_irr[160:220,550:670].ravel())]
    WAS = dT_irr[210:280,440:480].ravel()[~np.isnan(dT_irr[210:280,440:480].ravel())]
    WNA = dT_irr[237:300,100:150].ravel()[~np.isnan(dT_irr[237:300,100:150].ravel())]
    
    #Store the output in a list
    data = [GIL,CNA,EAS,MED,SAS,SEA,WAS,WNA]
    plt.figure(figsize=(12,8))
    bp = plt.boxplot(data,patch_artist=True)
    
    #Perform statistical tests: 1) Student's one-sample two-sided T-test to check
    #if the mean differs from 0, 2) Shapiro-Wilk test for normality to check if the
    #sample is normally distributed
    for a in range(0,len(data)):
      p_ttest = stats.ttest_1samp(data[a],0)[1]
      p_shapiro = stats.shapiro(data[a])[1]
      if p_ttest <= p_value:
        bp['boxes'][a].set(edgecolor='green',linewidth=3)
      if p_shapiro <= p_value:
        bp['boxes'][a].set(hatch='//')
      bp['boxes'][a].set(facecolor='white')
    
    #Some plotting commands
    meanpatch = mpatches.Patch(edgecolor='green',facecolor='white',linewidth=3,label=r'reject $H_0: \mu = 0$')
    normalpatch = mpatches.Patch(facecolor='white',linewidth=0,hatch='//',label=r'reject $H_0: X \sim N(\mu,\sigma^2)$')
    plt.xticks([1,2,3,4,5,6,7,8],['GIL','CNA','EAS','MED','SAS','SEA','WAS','WNA'])
    plt.xlabel('Region')
    plt.ylabel('T change [$^\circ$C]')
    plt.ylim([-2,2])
    
    #Plot the sample size on the bottom of the figure
    lenlist = [x.shape[0] for x in data]
    for i in range(0,len(data)):
      plt.text(i+1,-1.9,str(lenlist[i]),ha='center',va='center')
    
    plt.grid(True)
    plt.legend(handles=[meanpatch,normalpatch],loc='upper right')
    plt.title('method: %s, threshold: %1.2f, dataset: %s %s, p: %s, ref: %i-%i'%(method,thres_irr,datasource,temp_product,p_value,yr_start1,yr_end1))
    plt.show(block=False)      

  #Save Figure
  if figsave == True:
    if(method == 'regression' and def_regr == True):
      figpath = '/net/exo/landclim/wthiery/observational_analysis/Figures/CRU/%s/%s/dT_irr.%s.def.%1.2f.%s.%s.%s-%s.%s-%s.%s'%(temp_product,SREX_region,method,thres_irr,temp_product,SREX_region,yr_start1,yr_end1,yr_start2,yr_end2,figformat)
    else:
      figpath = '/net/exo/landclim/wthiery/observational_analysis/Figures/CRU/%s/%s/dT_irr.%s.alt.%1.2f.%s.%s.%s-%s.%s-%s.%s'%(temp_product,SREX_region,method,thres_irr,temp_product,SREX_region,yr_start1,yr_end1,yr_start2,yr_end2,figformat)
    if output_format in ['spatial plot', 'both']:
      if os.path.exists(figpath):
        print('Figure already exists at path: ' + figpath)
      else:  
        print('Saving figure as: ' + figpath)
        plt.savefig(figpath,bbox_inches='tight')
    elif output_format == 'boxplot':
      figpath = '/net/exo/landclim/wthiery/observational_analysis/Figures/CRU/boxplots/%s/dT_irr.boxplot.tests.%s.%1.2f.%s.%s-%s.%s-%s.%s'%(temp_product,method,thres_irr,temp_product,yr_start1,yr_end1,yr_start2,yr_end2,figformat)
      if os.path.exists(figpath):
        print('Figure already exists at path: ' + figpath)
      else:
        print('Saving figure as:' + figpath)
        plt.savefig(figpath,bbox_inches='tight')

elif(datasource in ["CESM","CRU_CESM"] or (datasource == "CRU" and temp_product in ['tmp_max', 'dtr_max', 'tmx_max', 'tmn_min','tmn_max'])):
  #Produce output based on user-specified options
  if(output_format in ['print', 'both']):
    if SREX_region == 'GIL':
      print np.nanmean(dT_irr)
    elif SREX_region == 'CNA':
      print np.nanmean(dT_irr[126:149,60:77])
    elif SREX_region == 'EAS':
      print np.nanmean(dT_irr[117:149,224:261])
    elif SREX_region == 'MED':
      print np.nanmean(dT_irr[128:144,136:177])
    elif SREX_region == 'SAS':
      c1 = np.count_nonzero(df_irr[101:133,192:220] >= thres_irr)
      c2 = np.count_nonzero(df_irr[117:133,220:225] >= thres_irr)
      if c2 == 0:
        print np.nanmean(dT_irr[101:133,192:220])
      else:
        b1 = np.nanmean(dT_irr[101:133,192:220])
        b2 = np.nanmean(dT_irr[117:133,220:225])
        stack = np.concatenate((np.tile(b1,(c1,1)),np.tile(b2,(c2,1))),axis=0)
        print np.nanmean(stack,axis=0)
    elif SREX_region == 'SEA':
      print np.nanmean(dT_irr[85:117,220:269])
    elif SREX_region == 'WAS':
      print np.nanmean(dT_irr[112:19,176:193])
    elif SREX_region == 'WNA':
      print np.nanmean(dT_irr[126:160,40:61])
 
  if(output_format in ['spatial plot', 'both']):
    #The plotting below follows the same procedure for every region (only differs for lat/lon bounds)
    #and is only explained with comments for GIL
    if SREX_region == 'GIL':
      plt.figure(figsize=(12,8)) #Define figure
      m = Basemap(projection='cyl',llcrnrlat=-60,urcrnrlat=90,llcrnrlon=-180,urcrnrlon=180,resolution='l') #Define global map
      m.drawcoastlines()
      m.drawmapboundary(fill_color='lightgray') #fill plot background
      m.fillcontinents(color='white',lake_color='lightgray') #color continents differently
      m.imshow(dT_irr[32::,:],vmin=-0.5,vmax=0.5,norm=norm,cmap='coolwarm',interpolation='none',zorder=4) #Plot dT_irr
      plot_SREX_shape(m,'CNA') #Plot SREX shapes
      plot_SREX_shape(m,'EAS')
      plot_SREX_shape(m,'MED')
      plot_SREX_shape(m,'SAS')
      plot_SREX_shape(m,'SEA')
      plot_SREX_shape(m,'WAS')
      plot_SREX_shape(m,'WNA')
      c = plt.colorbar(extend='both') #Plot colorbars
      c.set_label('$\Delta T_{irr}$ [$^\circ$C]')
      plt.title('Region: %s, method: %s, threshold: %1.2f, dataset: %s %s, \n ref: %i-%i'%(SREX_region,method,thres_irr,datasource,temp_product,yr_start1,yr_end1))
      plt.show(block=False)
    elif SREX_region == 'CNA':
      plt.figure(figsize=(12,8))
      m = Basemap(projection='cyl',llcrnrlat=18.566,urcrnrlat=60,llcrnrlon=-115,urcrnrlon=-75,resolution='l')
      m.drawcoastlines()
      m.drawmapboundary(fill_color='lightgray')
      m.fillcontinents(color='white',lake_color='lightgray')
      m.imshow(dT_irr[116:160,52:85],vmin=-0.5,vmax=0.5,norm=norm,cmap='coolwarm',interpolation='none',zorder=4)
      plot_SREX_shape(m,SREX_region)
      c = plt.colorbar(extend='both')
      c.set_label('$\Delta T_{irr}$ [$^\circ$C]')
      plt.title('Region: %s, method: %s, threshold: %1.2f, dataset: %s %s, \n ref: %i-%i'%(SREX_region,method,thres_irr,datasource,temp_product,yr_start1,yr_end1))
      plt.show(block=False)
    elif SREX_region == 'EAS':
      plt.figure(figsize=(12,8))
      m = Basemap(projection='cyl',llcrnrlat=10,urcrnrlat=60,llcrnrlon=90,urcrnrlon=155,resolution='l')
      m.drawcoastlines()
      m.drawmapboundary(fill_color='lightgray')
      m.fillcontinents(color='white',lake_color='lightgray')
      m.imshow(dT_irr[107:160,216:269],vmin=-0.5,vmax=0.5,norm=norm,cmap='coolwarm',interpolation='none',zorder=4)
      plot_SREX_shape(m,SREX_region)
      c = plt.colorbar(extend='both')
      c.set_label('$\Delta T_{irr}$ [$^\circ$C]')
      plt.title('Region: %s, method: %s, threshold: %1.2f, dataset: %s %s, \n ref: %i-%i'%(SREX_region,method,thres_irr,datasource,temp_product,yr_start1,yr_end1))
      plt.show(block=False)
    elif SREX_region == 'MED':
      plt.figure(figsize=(12,8))
      m = Basemap(projection='cyl',llcrnrlat=20,urcrnrlat=55,llcrnrlon=-20,urcrnrlon=50,resolution='l')
      m.drawcoastlines()
      m.drawmapboundary(fill_color='lightgray')
      m.fillcontinents(color='white',lake_color='lightgray')
      m.imshow(dT_irr[117:154,128:185],vmin=-0.5,vmax=0.5,norm=norm,cmap='coolwarm',interpolation='none',zorder=4)
      plot_SREX_shape(m,SREX_region)
      c = plt.colorbar(extend='both')
      c.set_label('$\Delta T_{irr}$ [$^\circ$C]')
      plt.title('Region: %s, method: %s, threshold: %1.2f, dataset: %s %s, \n ref: %i-%i'%(SREX_region,method,thres_irr,datasource,temp_product,yr_start1,yr_end1))
      plt.show(block=False)
    elif SREX_region == 'SAS':
      plt.figure(figsize=(12,8))
      m = Basemap(projection='cyl',llcrnrlat=-5,urcrnrlat=45,llcrnrlon=50,urcrnrlon=110,resolution='l')
      m.drawcoastlines()
      m.drawmapboundary(fill_color='lightgray')
      m.fillcontinents(color='white',lake_color='lightgray')
      plot_var = dT_irr[91:144,184:233]
      #plot_var[0:30,70:80] = np.nan
      m.imshow(plot_var,vmin=-0.5,vmax=0.5,norm=norm,cmap='coolwarm',interpolation='none',zorder=4)
      plot_SREX_shape(m,SREX_region)
      c = plt.colorbar(extend='both')
      c.set_label('$\Delta T_{irr}$ [$^\circ$C]')
      plt.title('Region: %s, method: %s, threshold: %1.2f, dataset: %s %s, \n ref: %i-%i'%(SREX_region,method,thres_irr,datasource,temp_product,yr_start1,yr_end1))
      plt.show(block=False)
    elif SREX_region == 'SEA':
      plt.figure(figsize=(12,8))
      m = Basemap(projection='cyl',llcrnrlat=-20,urcrnrlat=30,llcrnrlon=85,urcrnrlon=165,resolution='l')
      m.drawcoastlines()
      m.drawmapboundary(fill_color='lightgray')
      m.fillcontinents(color='white',lake_color='lightgray')
      m.imshow(dT_irr[75:128,212:277],vmin=-0.5,vmax=0.5,norm=norm,cmap='coolwarm',interpolation='none',zorder=4)
      plot_SREX_shape(m,SREX_region)
      c = plt.colorbar(extend='both')
      c.set_label('$\Delta T_{irr}$ [$^\circ$C]')
      plt.title('Region: %s, method: %s, threshold: %1.2f, dataset: %s %s, \n ref: %i-%i'%(SREX_region,method,thres_irr,datasource,temp_product,yr_start1,yr_end1))
      plt.show(block=False)
    elif SREX_region == 'WAS':
      plt.figure(figsize=(12,8))
      m = Basemap(projection='cyl',llcrnrlat=5,urcrnrlat=60,llcrnrlon=30,urcrnrlon=70,resolution='l')
      m.drawcoastlines()
      m.drawmapboundary(fill_color='lightgray')
      m.fillcontinents(color='white',lake_color='lightgray')
      m.imshow(dT_irr[101:160,168:201],vmin=-0.5,vmax=0.5,norm=norm,cmap='coolwarm',interpolation='none',zorder=4)
      plot_SREX_shape(m,SREX_region)
      c = plt.colorbar(extend='both')
      c.set_label('$\Delta T_{irr}$ [$^\circ$C]')
      plt.title('Region: %s, method: %s, threshold: %1.2f, dataset: %s %s, \n ref: %i-%i'%(SREX_region,method,thres_irr,datasource,temp_product,yr_start1,yr_end1))
      plt.show(block=False)
    elif SREX_region == 'WNA':
      plt.figure(figsize=(12,8))
      m = Basemap(projection='cyl',llcrnrlat=18.566,urcrnrlat=70,llcrnrlon=-140,urcrnrlon=-95,resolution='l')
      m.drawcoastlines()
      m.drawmapboundary(fill_color='lightgray')
      m.fillcontinents(color='white',lake_color='lightgray')
      m.imshow(dT_irr[116:170,32:69],vmin=-0.5,vmax=0.5,norm=norm,cmap='coolwarm',interpolation='none',zorder=4)
      plot_SREX_shape(m,SREX_region)
      c = plt.colorbar(extend='both')
      c.set_label('$\Delta T_{irr}$ [$^\circ$C]')
      plt.title('Region: %s, method: %s, threshold: %1.2f, dataset: %s %s, \n ref: %i-%i'%(SREX_region,method,thres_irr,datasource,temp_product,yr_start1,yr_end1))
      plt.show(block=False)
      
  if output_format == 'boxplot':
    #Store regional dT_irr output in 1D arrays
    GIL = dT_irr.ravel()[~np.isnan(dT_irr.ravel())]
    CNA = dT_irr[126:149,60:77].ravel()[~np.isnan(dT_irr[126:149,60:77].ravel())]
    EAS = dT_irr[117:149,224:261].ravel()[~np.isnan(dT_irr[117:149,224:261].ravel())]
    MED = dT_irr[128:144,136:177].ravel()[~np.isnan(dT_irr[128:144,136:177].ravel())]
    SAS1 = dT_irr[101:133,192:220].ravel()[~np.isnan(dT_irr[101:133,192:220].ravel())]
    SAS2 = dT_irr[117:133,220:225].ravel()[~np.isnan(dT_irr[117:133,220:225].ravel())]
    SAS = np.append(SAS1,SAS2)
    SEA = dT_irr[85:117,220:269].ravel()[~np.isnan(dT_irr[85:117,220:269].ravel())]
    WAS = dT_irr[112:149,176:193].ravel()[~np.isnan(dT_irr[112:149,176:193].ravel())]
    WNA = dT_irr[126:160,40:61].ravel()[~np.isnan(dT_irr[126:160,40:61].ravel())]
    
    #Arrange the regional arrays in a list
    data = [GIL,CNA,EAS,MED,SAS,SEA,WAS,WNA]
    plt.figure(figsize=(12,8))
    bp = plt.boxplot(data,patch_artist=True)
    
    #Perform statistical tests
    for a in range(0,len(data)):
      p_ttest = stats.ttest_1samp(data[a],0)[1]
      p_shapiro = stats.shapiro(data[a])[1]
      if p_ttest <= p_value:
        bp['boxes'][a].set(edgecolor='green',linewidth=3)
      if p_shapiro <= p_value:
        bp['boxes'][a].set(hatch='//')
      bp['boxes'][a].set(facecolor='white')
    
    #Some plotting commands
    meanpatch = mpatches.Patch(edgecolor='green',facecolor='white',linewidth=3,label=r'reject $H_0: \mu = 0$')
    normalpatch = mpatches.Patch(facecolor='white',linewidth=0,hatch='//',label=r'reject $H_0: X \sim N(\mu,\sigma^2)$')
    plt.xticks([1,2,3,4,5,6,7,8],['GIL','CNA','EAS','MED','SAS','SEA','WAS','WNA'])
    plt.xlabel('Region')
    plt.ylabel('T change [$^\circ$C]')
    plt.ylim([-2,2])
    
    lenlist = [x.shape[0] for x in data]
    for i in range(0,len(data)):
      plt.text(i+1,-1.9,str(lenlist[i]),ha='center',va='center')
    
    plt.grid(True)
    plt.legend(handles=[meanpatch,normalpatch],loc='upper right')
    plt.title('method: %s, threshold: %1.2f, dataset: %s %s, p: %s, ref: %i-%i'%(method,thres_irr,datasource,temp_product,p_value,yr_start1,yr_end1))
    plt.show(block=False)      

  #Save Figure
  if figsave == True:
    if(method == 'regression' and def_regr == True):
      figpath = '/net/exo/landclim/wthiery/observational_analysis/Figures/CRU/%s/%s/dT_irr.%s.def.%1.2f.%s.%s.%s-%s.%s-%s.%s'%(temp_product,SREX_region,method,thres_irr,temp_product,SREX_region,yr_start1,yr_end1,yr_start2,yr_end2,figformat)
    else:
      figpath = '/net/exo/landclim/wthiery/observational_analysis/Figures/CRU/%s/%s/dT_irr.%s.alt.%1.2f.%s.%s.%s-%s.%s-%s.%s'%(temp_product,SREX_region,method,thres_irr,temp_product,SREX_region,yr_start1,yr_end1,yr_start2,yr_end2,figformat)
    if output_format in ['spatial plot', 'both']:
      if os.path.exists(figpath):
        print('Figure already exists at path: ' + figpath)
      else:  
        print('Saving figure as: ' + figpath)
        plt.savefig(figpath,bbox_inches='tight')
    elif output_format == 'boxplot':
      figpath = '/net/exo/landclim/wthiery/observational_analysis/Figures/CRU/boxplots/%s/dT_irr.boxplot.tests.%s.%1.2f.%s.%s-%s.%s-%s.%s'%(temp_product,method,thres_irr,temp_product,yr_start1,yr_end1,yr_start2,yr_end2,figformat)
      if os.path.exists(figpath):
        print('Figure already exists at path: ' + figpath)
      else:
        print('Saving figure as:' + figpath)
        plt.savefig(figpath,bbox_inches='tight')



elif datasource == "E-OBS":
  #Print mean dT_irr
  if output_format in ['print', 'both']:
    if EOBS_region == 'EUR':
      print np.nanmean(dT_irr)
    elif EOBS_region == 'MED':
      print np.nanmean(dT_irr[19:79,122:322])
  
  #Plot results
  if output_format in ['spatial plot', 'both']:
    if EOBS_region == 'EUR':
      plt.figure(figsize=(12,8))
      m = Basemap(projection='cyl',llcrnrlat=25.375,urcrnrlat=75.375,llcrnrlon=-40.375,urcrnrlon=75.375,resolution='l')
      m.drawcoastlines()
      m.drawmapboundary(fill_color='lightgray')
      m.fillcontinents(color='white',lake_color='lightgray')
      m.imshow(dT_irr,vmin=-0.5,vmax=0.5,norm=norm,cmap='coolwarm',interpolation='none',zorder=4)
      plot_SREX_shape(m,'MED')
      c = plt.colorbar(extend='both')
      c.set_label('$\Delta T_{irr}$ [$^\circ$C]')
      plt.title('Region: %s, method: %s, threshold: %1.2f, dataset: %s %s'%(EOBS_region,method,thres_irr,datasource,temp_product))
      plt.show(block=False)
    elif EOBS_region == 'MED':
      plt.figure(figsize=(12,8))
      m = Basemap(projection='cyl',llcrnrlat=28.125,urcrnrlat=46.875,llcrnrlon=-11.875,urcrnrlon=41.875,resolution='l')
      m.drawcoastlines()
      m.drawmapboundary(fill_color='lightgray')
      m.fillcontinents(color='white',lake_color='lightgray')
      m.imshow(dT_irr[11:87,114:330],vmin=-0.5,vmax=0.5,norm=norm,cmap='coolwarm',interpolation='none',zorder=4)
      plot_SREX_shape(m,EOBS_region)
      c = plt.colorbar(extend='both')
      c.set_label('$\Delta T_{irr}$ [$^\circ$C]')
      plt.title('Region: %s, method: %s, threshold: %1.2f, dataset: %s %s'%(EOBS_region,method,thres_irr,datasource,temp_product))
      plt.show(block=False)
      
  if output_format == 'boxplot':
    EUR = dT_irr.ravel()[~np.isnan(dT_irr.ravel())]
    MED = dT_irr[19:79,122:322].ravel()[~np.isnan(dT_irr[19:79,122:322].ravel())]
    
    data = [EUR, MED]
    plt.figure(figsize=(12,8))
    bp = plt.boxplot(data,patch_artist=True)
    for a in range(0,len(data)):
      p_ttest = stats.ttest_1samp(data[a],0)[1]
      p_shapiro = stats.shapiro(data[a])[1]
      if p_ttest <= p_value:
        bp['boxes'][a].set(edgecolor='green',linewidth=3)
      if p_shapiro <= p_value:
        bp['boxes'][a].set(hatch='//')
      bp['boxes'][a].set(facecolor='white')
      
    meanpatch = mpatches.Patch(edgecolor='green',facecolor='white',linewidth=3,label=r'reject $H_0: \mu = 0$')
    normalpatch = mpatches.Patch(facecolor='white',linewidth=0,hatch='//',label=r'reject $H_0: X \sim N(\mu,\sigma^2)$')
    plt.xticks([1,2],['EUR','MED'])
    plt.xlabel('Region')
    plt.ylabel('T change [$^\circ$C]')
    plt.ylim([-1,1])
    
    lenlist = [x.shape[0] for x in data]
    for i in range(0,len(data)):
      plt.text(i+1,-1.9,str(lenlist[i]),ha='center',va='center')
    
    plt.grid(True)
    plt.legend(handles=[meanpatch,normalpatch],loc='upper right')
    plt.title('method: %s, threshold: %1.2f, dataset: %s %s, p: %s,ref: %i-%i'%(method,thres_irr,datasource,temp_product,p_value,yr_start1,yr_end1))
    plt.show(block=False)      

  #Save Figure
  if figsave == True:
    figpath = '/net/exo/landclim/wthiery/observational_analysis/Figures/E-OBS/%s/%s/dT_irr.%s.%1.2f.%s.%s.%s-%s.%s-%s.%s'%(temp_product,EOBS_region,method,thres_irr,temp_product,EOBS_region,yr_start1,yr_end1,yr_start2,yr_end2,figformat)
    if output_format in ['spatial plot','both']:
      if os.path.exists(figpath):
        print('Figure already exists at path: ' + figpath)
      else:
        print('Saving figure as:' + figpath)
        plt.savefig(figpath,bbox_inches='tight')
    elif output_format == 'boxplot':
      figpath = '/net/exo/landclim/wthiery/observational_analysis/Figures/E-OBS/boxplots/%s/dT_irr.boxplot.tests.%s.%1.2f.%s.%s-%s.%s-%s.%s'%(temp_product,method,thres_irr,temp_product,yr_start1,yr_end1,yr_start2,yr_end2,figformat)
      if os.path.exists(figpath):
        print('Figure already exists at path: ' + figpath)
      else:
        print('Saving figure as:' + figpath)
        plt.savefig(figpath,bbox_inches='tight')
    
elif datasource == "HadEX2":
  if output_format in ['print', 'both']:
    print np.nanmean(dT_irr)
  
  if output_format in ['spatial plot', 'both']:
    plt.figure(figsize=(12,8))
    m = Basemap(projection='cyl',llcrnrlat=-60,urcrnrlat=90,llcrnrlon=-180,urcrnrlon=176.25,resolution='l')
    m.drawcoastlines()
    m.drawmapboundary(fill_color='lightgray')
    m.fillcontinents(color='white',lake_color='lightgray')
    m.imshow(dT_irr[12::,:],vmin=-0.5,vmax=0.5,norm=norm,cmap='coolwarm',interpolation='none',zorder=4)
    plot_SREX_shape(m,'CNA')
    plot_SREX_shape(m,'EAS')
    plot_SREX_shape(m,'MED')
    plot_SREX_shape(m,'SAS')
    plot_SREX_shape(m,'SEA')
    plot_SREX_shape(m,'WAS')
    plot_SREX_shape(m,'WNA')
    c = plt.colorbar(extend='both')
    c.set_label('$\Delta T_{irr}$ [$^\circ$C]')
    plt.title('Region: %s, method: %s %s, threshold: %1.2f, dataset: %s %s, \n ref: %i-%i'%(SREX_region,method,response,thres_irr,datasource,temp_product,yr_start1,yr_end1))
    plt.show(block=False)
    
  if figsave == True:
    figpath = '/net/exo/landclim/wthiery/observational_analysis/Figures/HadEX2/%s/dT_irr.%s.alt.%s.%1.2f.%s.%s.%s.%s-%s.%s-%s.%s'%(SREX_region,method,response,thres_irr,datasource,temp_product,SREX_region,yr_start1,yr_end1,yr_start2,yr_start2,figformat)
    if os.path.exists(figpath):
      print('Figure already exists at path: ' + figpath)
    else:
      print('Saving figure as:' + figpath)
      plt.savefig(figpath,bbox_inches='tight')
  
    
