"""plot_T_out.py

author: Auke Visser
date: 24.10.2016

This script plots the output of extract_T_irr.py

"""

import netCDF4 as nc
import numpy as np
import scipy
import os
import matplotlib
matplotlib.rcParams['backend'] = "Qt4Agg"
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import maskoceans
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt

execfile('extract_T_irr.py')

T_out = extract_T_irr('CRU','Tmean','PC/PD',1901,1930,1981,2010)
plot_var = T_out[1,2,:,:] - T_out[0,2,:,:]
Tfile = nc.Dataset('/net/exo/landclim/data/dataset/CRUTS/v3.22/0.5deg_lat-lon_1m/original/cru_ts3.22.1901.2013.tmp.dat.nc','r')
lat = Tfile.variables['lat'][:]
lon = Tfile.variables['lon'][:]
lats = np.tile(lat,(lon.shape[0],1)).T
lons = np.tile(lon,(lat.shape[0],1))
m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,llcrnrlon=-180,urcrnrlon=180,resolution='l')

pv_oceanmask = maskoceans(lons,lats,plot_var,inlands=True,resolution='l',grid=10)

plt.figure()
m.drawcoastlines()
m.drawmapboundary(fill_color='lightgray')
m.imshow(pv_oceanmask,cmap="RdBu_r",vmin=-3,vmax=3)
c = plt.colorbar()
c.set_label('T change [$^\circ$ C]')
plt.show(block=False)