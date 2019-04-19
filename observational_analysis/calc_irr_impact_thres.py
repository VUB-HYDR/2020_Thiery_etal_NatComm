"""calc_irr_impact_thres.py

author: Auke Visser
date: 13.10.2016

This script calculates the impact on temperature
following the algorithm by Kumar et al. (2013) and 
Lejeune et al. (2016, in rev.)

This code is inspired by Quentin Lejeune's NCL version of the
Kumar algorithm.

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

execfile('extract_T_irr.py')
execfile('calc_irr_diff.py')

#Select big box size
#Choose an uneven number C as the number of grid cells in one direction, then:
#bl1 = (C-1)/2
#
#bl2 and bl3 should be chosen such that they are two larger than bl1 and bl2, respectively.

#Factor by which thres_irr should be multiplied to serve as criterion for df_irr splitting 
#in the big box
#u = 0.5

def calc_irr_impact_thres(datasource,temp_product,response,t_res,thres_irr,u,yr_start1,yr_end1,yr_start2,yr_end2):
  T_out = extract_T_irr(datasource,temp_product,response,t_res,yr_start1,yr_end1,yr_start2,yr_end2)
  lat,lon,f_irr = calc_irr_diff(datasource,response,yr_start1,yr_end1,yr_start2,yr_end2)

  #This script uses a variable big box size depending on the data source,
  #which is specified below
  if datasource == "HadEX2":
    bl1 = 2
    br1 = bl1+1
    bl2 = 3
    br2 = bl2+1
    bl3 = 4 #Only for lon biggest box
    br3 = bl3+1 #Only for lat
  elif datasource == "CRU":
    bl1 = 7
    br1 = bl1+1
    bl2 = 9
    br2 = bl2+1
    bl3 = 11 #Only for lon biggest box
    br3 = bl3+1 #Only for lat
  elif datasource in ["CESM", "CRU_CESM"]:
    bl1 = 5
    br1 = bl1+1
    bl2 = 6
    br2 = bl2+1
    bl3 = 7
    br3 = bl3+1
  
  #Define pixel selection criteria for the number of 
  #valid data entries inside the big box
  thres_nonan = np.ceil(0.3*(bl1+bl2)**2)
  thres_valid = np.ceil(0.1*(bl1+bl2)**2)
  
  #Pixel selection is performed differently for
  #PC/PD (based on change in irrigated fraction)
  #or PD (based on present-day irrigated fraction)
  if response == 'PC/PD':
    f_irr_thres = f_irr[1,:,:] - f_irr[0,:,:]
  elif response == 'PD':
    f_irr_thres = f_irr[1,:,:]
  
  #df_irr is the change in irrigated fraction, which
  #is used in executing the method
  df_irr = f_irr[1,:,:] - f_irr[0,:,:]
    
  #irrdata = nc.Dataset('/net/exo/landclim/wthiery/observational_analysis/Data/HID_regridded/hid_v1.0.f_irr_hyde_final_ir.0.25deg.10y.1900-2005.nc','r')
  #irr = irrdata.variables['aei_hyde_final_ir'][:] #IMPORTANT: remove when applying method to irrigation changes!!!
  #df_irr = irr[0,:,:] #IMPORTANT: remove when applying method to irrigation changes!!!
  print('Done importing changes in temperature and irrigated fraction')
  
  if(datasource == "CRU" and temp_product in ['tmp','dtr','tmn','tmx']):
    T_diff = T_out[1,:,:,:] - T_out[0,:,:,:]
    del T_out,lat,lon
    irr_impact = np.zeros((T_diff.shape[0],df_irr.shape[0],df_irr.shape[1]))
    dT_GC = np.zeros((T_diff.shape[0],df_irr.shape[0],df_irr.shape[1]))
    CtSt = np.zeros((T_diff.shape[0],df_irr.shape[0],df_irr.shape[1]))
    BB = 0
    count = 0
    
    irr_impact[:,:,:] = np.nan
    dT_GC[:,:,:] = np.nan
    CtSt[:,:,:] = np.nan
    if response == "PD":
      thres_irr *= u
    
    for k in range(0,T_diff.shape[0]): #loop over seasons/months
      for i in range(bl2,df_irr.shape[0]-bl2):
        for j in range(bl3,df_irr.shape[1]-bl3):
        
          #Start the analysis if a cell has experienced a change in irrigated fraction exceeding the threshold
          if(~np.isnan(T_diff[k,i,j]) and ~np.isnan(df_irr[i,j]) and f_irr_thres[i,j] >= thres_irr):
            irr_window_1 = df_irr[i-bl1:i+br1,j-bl1:j+br1]
            CtSt1 = float(np.count_nonzero(irr_window_1 <= thres_irr)) / float(np.count_nonzero(~np.isnan(irr_window_1)))
	    
            #The smallest big box size is selected if the number of cells exceeding or below the irrigation
            #threshold is not too skewed, and if there are sufficient valid data entries, as well as sufficient
            #cells above and below the f_irr threshold
            if(	0.35 <= CtSt1 <= 0.65 and np.count_nonzero(~np.isnan(irr_window_1)) >= thres_nonan and
                np.count_nonzero(irr_window_1[~np.isnan(irr_window_1)] <= thres_irr) >= thres_valid and
                np.count_nonzero(irr_window_1[~np.isnan(irr_window_1)] >= thres_irr) >= thres_valid
              ):
              
              BB = 1
	      
            else:
              irr_window_2 = df_irr[i-bl2:i+br2,j-bl2:j+br2]
              CtSt2 = float(np.count_nonzero(irr_window_2 <= thres_irr)) / float(np.count_nonzero(~np.isnan(irr_window_2)))
              
              #If the selection criteria are not met for the smallest big box, increase the search window size
              #and try again
              if(	0.35 <= CtSt2 <= 0.65 and np.count_nonzero(~np.isnan(irr_window_2)) > thres_nonan and
                  np.count_nonzero(irr_window_2[~np.isnan(irr_window_2)] <= thres_irr) >= thres_valid and
                  np.count_nonzero(irr_window_2[~np.isnan(irr_window_2)] >= thres_irr) >= thres_valid
                ):
                
                BB = 2
                
              else:
              
                irr_window_3 = df_irr[i-bl2:i+br2,j-bl3:j+br3]
                CtSt3 = float(np.count_nonzero(irr_window_3 <= thres_irr)) / float(np.count_nonzero(~np.isnan(irr_window_3)))
                
                #If the pixel selection criteria are not met again, increase the search window size once again
                #and select the size for which the ratio of cells above and below the threshold is closest to
                #1, still ensuring all the other criteria are met.
                if(CtSt3 < 0.35 or CtSt3 > 0.65):
                  
                  best_CtSt = np.nanmin([abs(CtSt1-0.5),abs(CtSt2-0.5),abs(CtSt3-0.5)])
                  
                  if(	abs(CtSt1-0.5) == best_CtSt and np.count_nonzero(~np.isnan(irr_window_1)) >= thres_nonan and
                      np.count_nonzero(irr_window_1[~np.isnan(irr_window_1)] <= thres_irr) >= thres_valid and
                      np.count_nonzero(irr_window_1[~np.isnan(irr_window_1)] >= thres_irr) >= thres_valid
                    ):
                    
                    BB = 1
                    
                  elif(	abs(CtSt2-0.5) == best_CtSt and np.count_nonzero(~np.isnan(irr_window_2)) >= thres_nonan and
                        np.count_nonzero(irr_window_2[~np.isnan(irr_window_2)] <= thres_irr) >= thres_valid and
                        np.count_nonzero(irr_window_2[~np.isnan(irr_window_2)] >= thres_irr) >= thres_valid 
                      ):
                    
                    BB = 2
                    
                  elif(	abs(CtSt3-0.5) == best_CtSt and np.count_nonzero(~np.isnan(irr_window_3)) >= thres_nonan and
                        np.count_nonzero(irr_window_3[~np.isnan(irr_window_3)] <= thres_irr) >= thres_valid and
                        np.count_nonzero(irr_window_3[~np.isnan(irr_window_3)] >= thres_irr) >= thres_valid
                      ):
                    
                    BB = 3
                    
                  elif(	0.35 <= CtSt3 <= 0.65 and np.count_nonzero(~np.isnan(irr_window_3)) >= thres_nonan and
                        np.count_nonzero(irr_window_3[~np.isnan(irr_window_3)] <= thres_irr) >= thres_valid and
                        np.count_nonzero(irr_window_3[~np.isnan(irr_window_3)] >= thres_irr) >= thres_valid
                      ):
                    
                    BB = 3
		  
	  #Calculate irrigation impact as:
	  #The mean of the cells that have experienced a decrease in irrigated fraction minus
	  #The mean of the cells that have not experienced a decrease in irrigated fraction
	  
          if BB == 3:
            T_window = T_diff[k,i-bl2:i+br2,j-bl3:j+br3]
            df = np.nanmean(irr_window_3[irr_window_3 >= thres_irr]) - np.nanmean(irr_window_3[irr_window_3 < thres_irr])
            CtSt[k,i,j] = CtSt3
            irr_impact[k,i,j] = np.nanmean(T_window[irr_window_3 >= thres_irr]) - np.nanmean(T_window[irr_window_3 < thres_irr])
            dT_GC[k,i,j] = (irr_impact[k,i,j] / df) * df_irr[i,j]
            count += 1
            
          elif BB == 2:
            T_window = T_diff[k,i-bl2:i+br2,j-bl2:j+br2]
            df = np.nanmean(irr_window_2[irr_window_2 >= thres_irr]) - np.nanmean(irr_window_2[irr_window_2 < thres_irr])
            CtSt[k,i,j] = CtSt2
            irr_impact[k,i,j] = np.nanmean(T_window[irr_window_2 >= thres_irr]) - np.nanmean(T_window[irr_window_2 < thres_irr])
            dT_GC[k,i,j] = (irr_impact[k,i,j] / df) * df_irr[i,j]
            count += 1
            
          elif BB == 1:
            T_window = T_diff[k,i-bl1:i+br1,j-bl1:j+br1]
            df = np.nanmean(irr_window_1[irr_window_1 >= thres_irr]) - np.nanmean(irr_window_1[irr_window_1 < thres_irr])
            CtSt[k,i,j] = CtSt1
            irr_impact[k,i,j] = np.nanmean(T_window[irr_window_1 >= thres_irr]) - np.nanmean(T_window[irr_window_1 < thres_irr])
            dT_GC[k,i,j] = (irr_impact[k,i,j] / df) * df_irr[i,j]
            count += 1
            
          elif BB == 0:
            irr_impact[k,i,j] = np.nan
            CtSt[k,i,j] = np.nan
            
          #Reinitialization
          BB = 0
      
      #print count
      count = 0
  
  #The procedure is mostly repeated for other variables, with the only difference being that 
  #there is just one temperature product to be compared (as opposed to 4 or 12 for seasonal
  #or monthly analysis with CRU). 
  elif((datasource == "CRU" and temp_product in ["tmp_max", "dtr_max", "tmx_max", "tmn_min","tmn_max"]) or datasource in ["E-OBS","HadEX2","CESM","CRU_CESM"]):
    T_diff = T_out[1,:,:] - T_out[0,:,:]
    del T_out,lat,lon
    irr_impact = np.zeros((df_irr.shape[0],df_irr.shape[1]))
    dT_GC = np.zeros((df_irr.shape[0],df_irr.shape[1]))
    CtSt = np.zeros((df_irr.shape[0],df_irr.shape[1]))
    BB = 0
    count = 0
    
    irr_impact[:,:] = np.nan
    dT_GC[:,:] = np.nan
    CtSt[:,:] = np.nan
    if response == "PD":
      thres_irr *= u

    
    for i in range(bl2,df_irr.shape[0]-bl2):
      for j in range(bl3,df_irr.shape[1]-bl3):
        
        #Start the analysis if a cell has experienced a change in irrigated fraction exceeding the threshold
        if(~np.isnan(T_diff[i,j]) and ~np.isnan(df_irr[i,j]) and f_irr_thres[i,j] >= thres_irr/u):
          irr_window_1 = df_irr[i-bl1:i+br1,j-bl1:j+br1]
          CtSt1 = float(np.count_nonzero(irr_window_1 <= thres_irr)) / float(np.count_nonzero(~np.isnan(irr_window_1)))

          if(	0.35 <= CtSt1 <= 0.65 and np.count_nonzero(~np.isnan(irr_window_1)) >= thres_nonan and
              np.count_nonzero(irr_window_1[~np.isnan(irr_window_1)] <= thres_irr) >= thres_valid and
              np.count_nonzero(irr_window_1[~np.isnan(irr_window_1)] >= thres_irr) >= thres_valid
            ):
            
            BB = 1
            
          else:
            irr_window_2 = df_irr[i-bl2:i+br2,j-bl2:j+br2]
            CtSt2 = float(np.count_nonzero(irr_window_2 <= thres_irr)) / float(np.count_nonzero(~np.isnan(irr_window_2)))
            
            if(	0.35 <= CtSt2 <= 0.65 and np.count_nonzero(~np.isnan(irr_window_2)) > thres_nonan and
                np.count_nonzero(irr_window_2[~np.isnan(irr_window_2)] <= thres_irr) >= thres_valid and
                np.count_nonzero(irr_window_2[~np.isnan(irr_window_2)] >= thres_irr) >= thres_valid
              ):
              
              BB = 2
              
            else:
              
              irr_window_3 = df_irr[i-bl2:i+br2,j-bl3:j+br3]
              CtSt3 = float(np.count_nonzero(irr_window_3 <= thres_irr)) / float(np.count_nonzero(~np.isnan(irr_window_3)))
              
              if(CtSt3 < 0.35 or CtSt3 > 0.65):
                
                best_CtSt = np.nanmin([abs(CtSt1-0.5),abs(CtSt2-0.5),abs(CtSt3-0.5)])
                
                if(	abs(CtSt1-0.5) == best_CtSt and np.count_nonzero(~np.isnan(irr_window_1)) >= thres_nonan and
                    np.count_nonzero(irr_window_1[~np.isnan(irr_window_1)] <= thres_irr) >= thres_valid and
                    np.count_nonzero(irr_window_1[~np.isnan(irr_window_1)] >= thres_irr) >= thres_valid
                  ):
                  
                  BB = 1
                  
                elif(	abs(CtSt2-0.5) == best_CtSt and np.count_nonzero(~np.isnan(irr_window_2)) >= thres_nonan and
                      np.count_nonzero(irr_window_2[~np.isnan(irr_window_2)] <= thres_irr) >= thres_valid and
                      np.count_nonzero(irr_window_2[~np.isnan(irr_window_2)] >= thres_irr) >= thres_valid 
                    ):
                  
                  BB = 2
                  
                elif(	abs(CtSt3-0.5) == best_CtSt and np.count_nonzero(~np.isnan(irr_window_3)) >= thres_nonan and
                      np.count_nonzero(irr_window_3[~np.isnan(irr_window_3)] <= thres_irr) >= thres_valid and
                      np.count_nonzero(irr_window_3[~np.isnan(irr_window_3)] >= thres_irr) >= thres_valid
                    ):
                  
                  BB = 3
                  
                elif(	0.35 <= CtSt3 <= 0.65 and np.count_nonzero(~np.isnan(irr_window_3)) >= thres_nonan and
                      np.count_nonzero(irr_window_3[~np.isnan(irr_window_3)] <= thres_irr) >= thres_valid and
                      np.count_nonzero(irr_window_3[~np.isnan(irr_window_3)] >= thres_irr) >= thres_valid
                    ):
                  
                  BB = 3
                  
        #Calculate irrigation impact as:
        #The mean of the cells that have experienced a decrease in irrigated fraction minus
        #The mean of the cells that have not experienced a decrease in irrigated fraction
        
        if BB == 3:
          T_window = T_diff[i-bl2:i+br2,j-bl3:j+br3]
          df = np.nanmean(irr_window_3[irr_window_3 >= thres_irr]) - np.nanmean(irr_window_3[irr_window_3 < thres_irr])
          #print np.nanmean(irr_window_3[irr_window_3 >= thres_irr]) == np.nanmean(irr_window_3[(irr_window_3 >= thres_irr) & ~np.isnan(irr_window_3)])
          #print np.nanmean(irr_window_3[irr_window_3 < thres_irr]) == np.nanmean(irr_window_3[(irr_window_3 < thres_irr) & ~np.isnan(irr_window_3)])
          #print np.count_nonzero(irr_window_3 >= thres_irr)
          CtSt[i,j] = CtSt3
          irr_impact[i,j] = np.nanmean(T_window[irr_window_3 >= thres_irr]) - np.nanmean(T_window[irr_window_3 < thres_irr])
          dT_GC[i,j] = (irr_impact[i,j] / df) * df_irr[i,j]
          count += 1
          
        elif BB == 2:
          T_window = T_diff[i-bl2:i+br2,j-bl2:j+br2]
          df = np.nanmean(irr_window_2[irr_window_2 >= thres_irr]) - np.nanmean(irr_window_2[irr_window_2 < thres_irr])
          #print np.nanmean(irr_window_2[irr_window_2 >= thres_irr]) == np.nanmean(irr_window_2[(irr_window_2 >= thres_irr) & ~np.isnan(irr_window_2)])
          #print np.nanmean(irr_window_2[irr_window_2 < thres_irr]) == np.nanmean(irr_window_2[(irr_window_2 < thres_irr) & ~np.isnan(irr_window_2)])
          #print np.count_nonzero(irr_window_2 >= thres_irr)
          CtSt[i,j] = CtSt2
          irr_impact[i,j] = np.nanmean(T_window[irr_window_2 >= thres_irr]) - np.nanmean(T_window[irr_window_2 < thres_irr])
          dT_GC[i,j] = (irr_impact[i,j] / df) * df_irr[i,j]
          count += 1
      
        elif BB == 1:
          T_window = T_diff[i-bl1:i+br1,j-bl1:j+br1]
          df = np.nanmean(irr_window_1[irr_window_1 >= thres_irr]) - np.nanmean(irr_window_1[irr_window_1 < thres_irr])
          #print np.nanmean(irr_window_1[irr_window_1 >= thres_irr]) == np.nanmean(irr_window_1[(irr_window_1 >= thres_irr) & ~np.isnan(irr_window_1)])
          #print np.nanmean(irr_window_1[irr_window_1 < thres_irr]) == np.nanmean(irr_window_1[(irr_window_1 < thres_irr) & ~np.isnan(irr_window_1)])
          #print np.count_nonzero(irr_window_1 >= thres_irr)
          CtSt[i,j] = CtSt1
          irr_impact[i,j] = np.nanmean(T_window[irr_window_1 >= thres_irr]) - np.nanmean(T_window[irr_window_1 < thres_irr])
          dT_GC[i,j] = (irr_impact[i,j] / df) * df_irr[i,j]
          count += 1
          
        elif BB == 0:
          irr_impact[i,j] = np.nan
          CtSt[i,j] = np.nan
          
        #Reinitialization
        BB = 0
    
    #print count
    count = 0
    
  return(irr_impact,CtSt,dT_GC)

#The lines below are used to test the function, and should be commented out during regular
#execution of the algorithm.
      
#dT_irr, CtSt = calc_irr_impact_thres('CRU','tmp','PC/PD','seasonal',0.15,1901,1930,1981,2010)
#dT_irr, CtSt = calc_irr_impact_thres('E-OBS','tx','PC/PD','seasonal',0.1,1951,1980,1981,2010)
#dT_irr, CtSt = calc_irr_impact_thres('CRU','tmp_max','PC/PD','seasonal',0.15,1901,1930,1981,2010)
#dT_irr, _, dT_GC = calc_irr_impact_thres('CRU','tmp','PC/PD','monthly',0.15,1901,1930,1981,2010)
#dT_irr, _, dT_GC = calc_irr_impact_thres('CRU','tmx_max','PD','seasonal',0.15,1901,1930,1981,2010)
#dT_irr, _, dT_GC = calc_irr_impact_thres('CESM','TREFHTMX','PD','monthly',0.25,0.5,1901,1930,1981,2010)
#dT_irr, _, dT_GC = calc_irr_impact_thres('CRU_CESM','tmx_max','PD','monthly',0.25,0.5,1901,1930,1981,2010)
