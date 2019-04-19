"""calc_irr_impact_regr.py

author: Auke Visser
date: 27.10.2016

This script calculates the impact of irrigation on temperature
following the algorithm by Kumar et al. (2013) and Lejeune et 
al. (2016, in rev.)


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

execfile('extract_T_irr.py')
execfile('calc_irr_diff.py')
execfile('calc_tf_diff.py')

def calc_irr_impact_regr(datasource,temp_product,response,t_res,thres_irr,def_regr,yr_start1,yr_end1,yr_start2,yr_end2):
  T_out = extract_T_irr(datasource,temp_product,response,t_res,yr_start1,yr_end1,yr_start2,yr_end2)
  lat,lon,f_irr = calc_irr_diff(datasource,response,yr_start1,yr_end1,yr_start2,yr_end2)

  #Define the big box size for the different datasets
  if datasource == "HadEX2":
    bl = 2
    br = bl+1
  elif datasource == "CRU":
    bl = 7
    br = bl+1
  elif datasource in ["CESM","CRU_CESM","CRUv4_CESM","Berkeley_CESM"]:
    bl = 7
    br = bl+1
  
  #Define the pixel selection criterion based on present-day f_irr (PD)
  #or changes in irrigated fraction
  if response == 'PC/PD':
    f_irr_crit = f_irr[1,:,:] - f_irr[0,:,:]
  elif response == 'PD':
    f_irr_crit = f_irr[1,:,:]
    
  df_irr = f_irr[1,:,:] - f_irr[0,:,:]
    
  print('Done importing changes in temperature and irrigated fraction')
  
  if datasource == 'CRU':
    if len(T_out.shape) == 4:
      
      #Define the difference in temperature as the change in 
      #f_irr between the PC and PD period
      T_diff = T_out[1,:,:,:] - T_out[0,:,:,:]
      del T_out
      
      #Define some variables
      irr_impact  = np.zeros((T_diff.shape[0],df_irr.shape[0],df_irr.shape[1])) #where dT_irr is stored
      cod         = np.zeros((T_diff.shape[0],df_irr.shape[0],df_irr.shape[1])) #r^2 for the multi-linear regression model
      local_eff   = np.zeros((T_diff.shape[0],df_irr.shape[0],df_irr.shape[1])) #dT_tot at the pixel of interest
      region_eff  = np.zeros((T_diff.shape[0],df_irr.shape[0],df_irr.shape[1])) #average dT_tot in the big box
      lats        = np.tile(lat,(lon.shape[0],1)).T
      lons        = np.tile(lon,(lat.shape[0],1))
      defrate     = calc_tf_diff(datasource,response,yr_start1,yr_end1,yr_start2,yr_end2) #
      or_file     = nc.Dataset('/net/exo/landclim/wthiery/observational_analysis/Data/CLMdata_topography.0.5deg.nc','r') #Find elevation file
      orog        = or_file.variables['HSURF'][:] #extract elevation and store in array
      
      #Set all values to NaN, so that all values that are not selected for
      #analysis are discarded automatically when calculating results
      irr_impact[:,:,:] = np.nan
      cod[:,:,:]        = np.nan
      local_eff[:,:,:]  = np.nan
      region_eff[:,:,:] = np.nan
      ind_values = []
      
      for k in range(0,T_diff.shape[0]): #loop over seasons/months
        for i in range(bl,irr_impact.shape[1]-bl):
          for j in range(bl,irr_impact.shape[2]-bl):
            
            #Pixel selection: at least 60% valid data points and 8% points where f_irr (or df_irr,
            #depending on 'response' defined when calling the function) exceeds 0.
            if(	np.count_nonzero(~np.isnan(f_irr_crit[i-bl:i+br,j-bl:j+br])) >= np.ceil(0.6*(bl+br)**2) and
            np.count_nonzero(~np.isnan(T_diff[k,i-bl:i+br,j-bl:j+br])) >= np.ceil(0.6*(bl+br)**2) and
            np.count_nonzero(f_irr_crit[i-bl:i+br,j-bl:j+br] != 0) >= np.ceil(0.08*(bl+br)**2) and
            f_irr_crit[i,j] >= thres_irr
            ):
	      
              if def_regr == True:
                #Define [n_nonan x 5] array with effects of irr_rate, orog, lat and lon
                df_irr_box = df_irr[i-bl:i+br,j-bl:j+br]
                dT_box = T_diff[k,i-bl:i+br,j-bl:j+br]
                df_def_box = defrate[i-bl:i+br,j-bl:j+br]
                irr_box = df_irr_box.ravel()[~np.isnan(df_irr_box.ravel()) & ~np.isnan(dT_box.ravel()) & ~np.isnan(df_def_box.ravel())]
                dT_nonan = dT_box.ravel()[~np.isnan(df_irr_box.ravel()) & ~np.isnan(dT_box.ravel()) & ~np.isnan(df_def_box.ravel())]
		
                arr_eff = np.zeros((irr_box.shape[0],5))
                arr_eff[:,0] = irr_box
                arr_eff[:,1] = orog[i-bl:i+br,j-bl:j+br].ravel()[~np.isnan(df_irr_box.ravel()) & ~np.isnan(dT_box.ravel()) & ~np.isnan(df_def_box.ravel())]
                arr_eff[:,2] = lats[i-bl:i+br,j-bl:j+br].ravel()[~np.isnan(df_irr_box.ravel()) & ~np.isnan(dT_box.ravel()) & ~np.isnan(df_def_box.ravel())]
                arr_eff[:,3] = lons[i-bl:i+br,j-bl:j+br].ravel()[~np.isnan(df_irr_box.ravel()) & ~np.isnan(dT_box.ravel()) & ~np.isnan(df_def_box.ravel())]
                arr_eff[:,4] = defrate[i-bl:i+br,j-bl:j+br].ravel()[~np.isnan(df_irr_box.ravel()) & ~np.isnan(dT_box.ravel()) & ~np.isnan(df_def_box.ravel())]
                
                ind_values.append(np.vstack((arr_eff[:,0],dT_nonan)))
                
                #Define the multi-linear regression model
                reg = linear_model.LinearRegression()
                reg.fit(arr_eff,dT_nonan)
                
                #Extract the regression coefficient for irrigation and multiply it by the 
                #defrate in the cell of interest
                irr_impact[k,i,j] = reg.coef_[0] * df_irr[i,j]
                cod[k,i,j] = reg.score(arr_eff,dT_nonan)
                local_eff[k,i,j] = T_diff[k,i,j]
                region_eff[k,i,j] = np.nanmean(dT_nonan)
                
              elif def_regr == False:
                #Define [n_nonan x 4] array with effects of irr_rate, orog, lat and lon
                df_irr_box = df_irr[i-bl:i+br,j-bl:j+br]
                dT_box = T_diff[k,i-bl:i+br,j-bl:j+br]
                irr_box = df_irr_box.ravel()[~np.isnan(df_irr_box.ravel()) & ~np.isnan(dT_box.ravel())]
                dT_nonan = dT_box.ravel()[~np.isnan(df_irr_box.ravel()) & ~np.isnan(dT_box.ravel())]
                
                arr_eff = np.zeros((irr_box.shape[0],4))
                arr_eff[:,0] = irr_box
                arr_eff[:,1] = orog[i-bl:i+br,j-bl:j+br].ravel()[~np.isnan(df_irr_box.ravel()) & ~np.isnan(dT_box.ravel())]
                arr_eff[:,2] = lats[i-bl:i+br,j-bl:j+br].ravel()[~np.isnan(df_irr_box.ravel()) & ~np.isnan(dT_box.ravel())]
                arr_eff[:,3] = lons[i-bl:i+br,j-bl:j+br].ravel()[~np.isnan(df_irr_box.ravel()) & ~np.isnan(dT_box.ravel())]
                
                ind_values.append(np.vstack((arr_eff[:,0],dT_nonan)))
                
                #Define the multi-linear regression model
                reg = linear_model.LinearRegression()
                reg.fit(arr_eff,dT_nonan)
                
                #Extract the regression coefficient for irrigation and multiply it by the 
                #defrate in the cell of interest
                irr_impact[k,i,j] = reg.coef_[0] * df_irr[i,j]
                cod[k,i,j] = reg.score(arr_eff,dT_nonan)
                local_eff[k,i,j] = T_diff[k,i,j]
                region_eff[k,i,j] = np.nanmean(dT_nonan)
                
    elif len(T_out.shape) == 3:
      
      T_diff = T_out[1,:,:] - T_out[0,:,:]
      
      del T_out
      irr_impact = np.zeros((df_irr.shape[0],df_irr.shape[1]))
      cod = np.zeros((df_irr.shape[0],df_irr.shape[1]))
      local_eff = np.zeros((df_irr.shape[0],df_irr.shape[1]))
      region_eff = np.zeros((df_irr.shape[0],df_irr.shape[1]))
      lats = np.tile(lat,(lon.shape[0],1)).T
      lons = np.tile(lon,(lat.shape[0],1))
      defrate = calc_tf_diff(datasource,response,yr_start1,yr_end1,yr_start2,yr_end2)
      or_file = nc.Dataset('/net/exo/landclim/wthiery/observational_analysis/Data/CLMdata_topography.0.5deg.nc','r')
      orog = or_file.variables['HSURF'][:]
      
      irr_impact[:,:] = np.nan
      cod[:,:] = np.nan
      local_eff[:,:] = np.nan
      region_eff[:,:] = np.nan
      ind_values = []
      
      for i in range(bl,irr_impact.shape[0]-bl):
        for j in range(bl,irr_impact.shape[1]-bl):
        
          if(	np.count_nonzero(~np.isnan(df_irr[i-bl:i+br,j-bl:j+br])) >= np.ceil(0.6*(bl+br)**2) and
            np.count_nonzero(~np.isnan(T_diff[i-bl:i+br,j-bl:j+br])) >= np.ceil(0.6*(bl+br)**2) and
            np.count_nonzero(df_irr[i-bl:i+br,j-bl:j+br] != 0) >= np.ceil(0.08*(bl+br)**2) and
            f_irr_crit[i,j] >= thres_irr
            ):
            
            if def_regr == True:
              #Define [n_nonan x 5] array with effects of irr_rate, orog, lat and lons
              df_irr_box = df_irr[i-bl:i+br,j-bl:j+br]
              dT_box = T_diff[i-bl:i+br,j-bl:j+br]
              df_def_box = defrate[i-bl:i+br,j-bl:j+br]
              irr_box = df_irr_box.ravel()[~np.isnan(df_irr_box.ravel()) & ~np.isnan(dT_box.ravel()) & ~np.isnan(df_def_box.ravel())]
              dT_nonan = dT_box.ravel()[~np.isnan(df_irr_box.ravel()) & ~np.isnan(dT_box.ravel()) & ~np.isnan(df_def_box.ravel())]
              
              arr_eff = np.zeros((irr_box.shape[0],5))
              arr_eff[:,0] = irr_box
              arr_eff[:,1] = orog[i-bl:i+br,j-bl:j+br].ravel()[~np.isnan(df_irr_box.ravel()) & ~np.isnan(dT_box.ravel()) & ~np.isnan(df_def_box.ravel())]
              arr_eff[:,2] = lats[i-bl:i+br,j-bl:j+br].ravel()[~np.isnan(df_irr_box.ravel()) & ~np.isnan(dT_box.ravel()) & ~np.isnan(df_def_box.ravel())]
              arr_eff[:,3] = lons[i-bl:i+br,j-bl:j+br].ravel()[~np.isnan(df_irr_box.ravel()) & ~np.isnan(dT_box.ravel()) & ~np.isnan(df_def_box.ravel())]
              arr_eff[:,4] = defrate[i-bl:i+br,j-bl:j+br].ravel()[~np.isnan(df_irr_box.ravel()) & ~np.isnan(dT_box.ravel()) & ~np.isnan(df_def_box.ravel())]
              
              #ind_values.append(np.vstack((arr_eff[:,0],dT_nonan)))
              
              #Define the multi-linear regression model
              if arr_eff.shape[0] >= np.ceil(0.6*(bl+br)**2):
                reg = linear_model.LinearRegression()
                reg.fit(arr_eff,dT_nonan)
                
                #Extract the regression coefficient for irrigation and multiply it by the 
                #irr_rate in the cell of interest
                irr_impact[i,j] = reg.coef_[0] * df_irr[i,j]
                cod[i,j] = reg.score(arr_eff,dT_nonan)
                local_eff[i,j] = T_diff[i,j]
                region_eff[i,j] = np.nanmean(dT_nonan)
        
            elif def_regr == False:
              #Define [n_nonan x 4] array with effects of irr_rate, orog, lat and lons
              df_irr_box = df_irr[i-bl:i+br,j-bl:j+br]
              dT_box = T_diff[i-bl:i+br,j-bl:j+br]
              irr_box = df_irr_box.ravel()[~np.isnan(df_irr_box.ravel()) & ~np.isnan(dT_box.ravel())]
              dT_nonan = dT_box.ravel()[~np.isnan(df_irr_box.ravel()) & ~np.isnan(dT_box.ravel())]
              
              arr_eff = np.zeros((irr_box.shape[0],4))
              arr_eff[:,0] = irr_box
              arr_eff[:,1] = orog[i-bl:i+br,j-bl:j+br].ravel()[~np.isnan(df_irr_box.ravel()) & ~np.isnan(dT_box.ravel())]
              arr_eff[:,2] = lats[i-bl:i+br,j-bl:j+br].ravel()[~np.isnan(df_irr_box.ravel()) & ~np.isnan(dT_box.ravel())]
              arr_eff[:,3] = lons[i-bl:i+br,j-bl:j+br].ravel()[~np.isnan(df_irr_box.ravel()) & ~np.isnan(dT_box.ravel())]
            
              ind_values.append(np.vstack((arr_eff[:,0],dT_nonan)))
	      
              #Define the multi-linear regression model
              if arr_eff.shape[0] >= np.ceil(0.6*(bl+br)**2):

                reg = linear_model.LinearRegression()
                reg.fit(arr_eff,dT_nonan)
                
                #Extract the regression coefficient for irrigation and multiply it by the 
                #irr_rate in the cell of interest
                irr_impact[i,j] = reg.coef_[0] * df_irr[i,j]
                cod[i,j] = reg.score(arr_eff,dT_nonan)
                local_eff[i,j] = T_diff[i,j]
                region_eff[i,j] = np.nanmean(dT_nonan)
                
  elif datasource == 'E-OBS':
    
    T_diff = T_out[1,:,:] - T_out[0,:,:]
      
    del T_out
    irr_impact = np.zeros((df_irr.shape[0],df_irr.shape[1]))
    cod = np.zeros((df_irr.shape[0],df_irr.shape[1]))
    local_eff = np.zeros((df_irr.shape[0],df_irr.shape[1]))
    region_eff = np.zeros((df_irr.shape[0],df_irr.shape[1]))                                                   
    lats = np.tile(lat,(lon.shape[0],1)).T
    lons = np.tile(lon,(lat.shape[0],1))
    #!!! File not yet available, implement tomorrow!!!
    or_file = nc.Dataset('/net/exo/landclim/wthiery/observational_analysis/Data/CLMdata_topography_regridded.0.25deg.nc','r')
    orog = or_file.variables['HSURF'][:]
    
    irr_impact[:,:] = np.nan
    cod[:,:] = np.nan
    local_eff[:,:] = np.nan
    region_eff[:,:] = np.nan
    ind_values = []
    
    for i in range(bl,irr_impact.shape[0]-bl):
      for j in range(bl,irr_impact.shape[1]-bl):
        
        if(	np.count_nonzero(~np.isnan(f_irr_crit[i-bl:i+br,j-bl:j+br])) >= np.ceil(0.6*(bl+br)**2) and
          np.count_nonzero(~np.isnan(T_diff[i-bl:i+br,j-bl:j+br])) >= np.ceil(0.6*(bl+br)**2) and
          np.count_nonzero(f_irr_crit[i-bl:i+br,j-bl:j+br] != 0) >= np.ceil(0.08*(bl+br)**2) and
          f_irr_crit[i,j] >= thres_irr
          ):
            
          #Define [n_nonan x 4] array with effects of irr_rate, orog, lat and lons
          df_irr_box = df_irr[i-bl:i+br,j-bl:j+br]
          dT_box = T_diff[i-bl:i+br,j-bl:j+br]
          irr_box = df_irr_box.ravel()[~np.isnan(df_irr_box.ravel()) & ~np.isnan(dT_box.ravel())]
          dT_nonan = dT_box.ravel()[~np.isnan(df_irr_box.ravel()) & ~np.isnan(dT_box.ravel())]
          
          arr_eff = np.zeros((irr_box.shape[0],4))
          arr_eff[:,0] = irr_box
          arr_eff[:,1] = orog[i-bl:i+br,j-bl:j+br].ravel()[~np.isnan(df_irr_box.ravel()) & ~np.isnan(dT_box.ravel())]
          arr_eff[:,2] = lats[i-bl:i+br,j-bl:j+br].ravel()[~np.isnan(df_irr_box.ravel()) & ~np.isnan(dT_box.ravel())]
          arr_eff[:,3] = lons[i-bl:i+br,j-bl:j+br].ravel()[~np.isnan(df_irr_box.ravel()) & ~np.isnan(dT_box.ravel())]
          
          ind_values.append(np.vstack((arr_eff[:,0],dT_nonan)))
          
          #Define the multi-linear regression model
          reg = linear_model.LinearRegression()
          reg.fit(arr_eff,dT_nonan)
          
          #Extract the regression coefficient for irrigation and multiply it by the 
          #irr_rate in the cell of interest
          irr_impact[i,j] = reg.coef_[0] * df_irr[i,j]
          cod[i,j] = reg.score(arr_eff,dT_nonan)
          local_eff[i,j] = T_diff[i,j]
          region_eff[i,j] = np.nanmean(dT_nonan)
          
  elif datasource == 'HadEX2':
    
    T_diff = T_out[1,:,:] - T_out[0,:,:]
    
    del T_out
    irr_impact = np.zeros((df_irr.shape[0],df_irr.shape[1]))
    cod = np.zeros((df_irr.shape[0],df_irr.shape[1]))
    local_eff = np.zeros((df_irr.shape[0],df_irr.shape[1]))
    region_eff = np.zeros((df_irr.shape[0],df_irr.shape[1]))
    lats = np.tile(lat,(lon.shape[0],1)).T
    lons = np.tile(lon,(lat.shape[0],1))
    #!!! File not yet available, implement tomorrow!!!
    or_file = nc.Dataset('/net/exo/landclim/wthiery/observational_analysis/Data/CLMdata_topography.2.5x3.75deg.nc','r')
    orog = or_file.variables['HSURF'][:]
    
    irr_impact[:,:] = np.nan
    cod[:,:] = np.nan
    local_eff[:,:] = np.nan
    region_eff[:,:] = np.nan
    ind_values = []
    
    for i in range(bl,irr_impact.shape[0]-bl):
      for j in range(bl,irr_impact.shape[1]-bl):
        
        if(	np.count_nonzero(~np.isnan(f_irr_crit[i-bl:i+br,j-bl:j+br])) >= np.ceil(0.6*(bl+br)**2) and
            np.count_nonzero(~np.isnan(T_diff[i-bl:i+br,j-bl:j+br])) >= np.ceil(0.6*(bl+br)**2) and
            np.count_nonzero(f_irr_crit[i-bl:i+br,j-bl:j+br] != 0) >= np.ceil(0.08*(bl+br)**2) and
            f_irr_crit[i,j] >= thres_irr
          ):
          
          #Define [n_nonan x 4] array with effects of irr_rate, orog, lat and lons
          df_irr_box = df_irr[i-bl:i+br,j-bl:j+br]
          dT_box = T_diff[i-bl:i+br,j-bl:j+br]
          irr_box = df_irr_box.ravel()[~np.isnan(df_irr_box.ravel()) & ~np.isnan(dT_box.ravel())]
          dT_nonan = dT_box.ravel()[~np.isnan(df_irr_box.ravel()) & ~np.isnan(dT_box.ravel())]
          
          arr_eff = np.zeros((irr_box.shape[0],4))
          arr_eff[:,0] = irr_box
          arr_eff[:,1] = orog[i-bl:i+br,j-bl:j+br].ravel()[~np.isnan(df_irr_box.ravel()) & ~np.isnan(dT_box.ravel())]
          arr_eff[:,2] = lats[i-bl:i+br,j-bl:j+br].ravel()[~np.isnan(df_irr_box.ravel()) & ~np.isnan(dT_box.ravel())]
          arr_eff[:,3] = lons[i-bl:i+br,j-bl:j+br].ravel()[~np.isnan(df_irr_box.ravel()) & ~np.isnan(dT_box.ravel())]
          
          ind_values.append(np.vstack((arr_eff[:,0],dT_nonan)))
          
          #Define the multi-linear regression model
          reg = linear_model.LinearRegression()
          reg.fit(arr_eff,dT_nonan)
          
          ##Extract the regression coefficient for irrigation and multiply it by the 
          #irr_rate in the cell of interest
          irr_impact[i,j] = reg.coef_[0] * df_irr[i,j]
          cod[i,j] = reg.score(arr_eff,dT_nonan)
          local_eff[i,j] = T_diff[i,j]
          region_eff[i,j] = np.nanmean(dT_nonan)
          
  elif datasource in ["CESM","CRU_CESM","CRUv4_CESM","Berkeley_CESM"]:
    
    T_diff = T_out[1,:,:] - T_out[0,:,:]
    
    del T_out
    irr_impact = np.zeros((df_irr.shape[0],df_irr.shape[1]))
    cod = np.zeros((df_irr.shape[0],df_irr.shape[1]))
    local_eff = np.zeros((df_irr.shape[0],df_irr.shape[1]))
    region_eff = np.zeros((df_irr.shape[0],df_irr.shape[1]))
    lats = np.tile(lat,(lon.shape[0],1)).T
    lons = np.tile(lon,(lat.shape[0],1))
    or_file = nc.Dataset('/net/cfc/landclim1/wthiery/cesm_output/f.e122.F1850PDC5.f09_g16.control-io192.001/lnd/hist/year1976/f.e122.F1850PDC5.f09_g16.control-io192.001.clm2.h1.1976-01-01-00000.nc','r')
    orog = or_file.variables['topo'][:]
    
    irr_impact[:,:] = np.nan
    cod[:,:] = np.nan
    local_eff[:,:] = np.nan
    region_eff[:,:] = np.nan
    ind_values = []
    
    for i in range(bl,irr_impact.shape[0]-bl):
      for j in range(br,irr_impact.shape[1]-bl):
        
        if( np.count_nonzero(~np.isnan(f_irr_crit[i-bl:i+br,j-bl:j+br])) >= np.ceil(0.6*(bl+br)**2) and
            np.count_nonzero(~np.isnan(T_diff[i-bl:i+br,j-bl:j+br])) >= np.ceil(0.6*(bl+br)**2) and
            np.count_nonzero(f_irr_crit[i-bl:i+br,j-bl:j+br] != 0) >= np.ceil(0.08*(bl+br)**2) and
            f_irr_crit[i,j] >= thres_irr
          ):
          
          #Define [n_nonan x 4] array with effects of irr_rate, orog, lat and lon
          df_irr_box = df_irr[i-bl:i+br,j-bl:j+br]
          dT_box = T_diff[i-bl:i+br,j-bl:j+br]
          irr_box = df_irr_box.ravel()[~np.isnan(df_irr_box.ravel()) & ~np.isnan(dT_box.ravel())]
          dT_nonan = dT_box.ravel()[~np.isnan(df_irr_box.ravel()) & ~np.isnan(dT_box.ravel())]
          
          arr_eff = np.zeros((irr_box.shape[0],4))
          arr_eff[:,0] = irr_box
          arr_eff[:,1] = orog[i-bl:i+br,j-bl:j+br].ravel()[~np.isnan(df_irr_box.ravel()) & ~np.isnan(dT_box.ravel())]
          arr_eff[:,2] = lats[i-bl:i+br,j-bl:j+br].ravel()[~np.isnan(df_irr_box.ravel()) & ~np.isnan(dT_box.ravel())]
          arr_eff[:,3] = lons[i-bl:i+br,j-bl:j+br].ravel()[~np.isnan(df_irr_box.ravel()) & ~np.isnan(dT_box.ravel())]
          
          ind_values.append(np.vstack((arr_eff[:,0],dT_nonan)))
          
          #Define the multi-linear regression model
          if arr_eff.shape[0] >= np.ceil(0.6*(bl+br)**2):
            
            reg = linear_model.LinearRegression()
            reg.fit(arr_eff,dT_nonan)
            
            #Extract the regression coefficient for irrigation and multiply it by the irr_rate in the cell of interest
            irr_impact[i,j] = reg.coef_[0] * df_irr[i,j]
            cod[i,j] = reg.score(arr_eff,dT_nonan)
            local_eff[i,j] = T_diff[i,j]
            region_eff[i,j] = np.nanmean(dT_nonan)
  
  return(irr_impact,ind_values,cod,local_eff,region_eff)

#The lines below are used to test the function, and should be commented out during regular
#execution of the algorithm.

#dT_irr, indvals = calc_irr_impact_regr('CRU','tmp','PC/PD','seasonal',0.1,True,1901,1930,1981,2010)
#dT_irr_regr = calc_irr_impact_regr('E-OBS','tx','PC/PD','seasonal',0.1,True,1951,1980,1981,2010)
#dT_irr_regr = calc_irr_impact_regr('CRU','dtr_max','PC/PD','seasonal',0.1,True,1901,1930,1981,2010)
#dT_irr, _,_,_,_ = calc_irr_impact_regr('CRU','tmx','PD','monthly',0.4,False,1901,1930,1981,2010)
#dT_irr, _,_,_,_ = calc_irr_impact_regr('CRU','tmx_max','PC/PD','monthly',0.5,True,1901,1930,1981,2010)
#dT_irr,_,_,_,_ = calc_irr_impact_regr('CRU','tmx_max','PD','monthly',0.3,False,1901,1930,1981,2010)
#dT_irr,_,cod = calc_irr_impact_regr('HadEX2','Ann','PD','seasonal',0.05,True,1956,1980,1985,2009)
#dT_irr, _,_,_,_ = calc_irr_impact_regr('CESM','TREFHTMX','PD','seasonal',0.2,False,1901,1930,1981,2010)
#dT_irr, _,_,_,_ = calc_irr_impact_regr('CRU_CESM','tmx_max','PD','seasonal',0.2,False,1901,1930,1981,2010)
