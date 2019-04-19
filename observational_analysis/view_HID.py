"""view_HID.py

author: Auke Visser
date: 6.10.2016

This script is used for a first exploratory look on 
the Historical Irrigation Dataset"""

import netCDF4 as nc
import numpy as np
import scipy
import os
import matplotlib
matplotlib.rcParams['backend'] = "Qt4Agg"
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import matplotlib.colors as colors

# Declare path and other variables
path = '/net/exo/landclim/wthiery/observational_analysis/Data/HID_regridded/'
#product = 'aei_earthstat_ir' 
product = 'aei_hyde_final_ir'
filename = 'hid_v1.0.f_irr_hyde_final_ir.0.5deg.10y.1900-2005.nc'

yr_start=1951
yr_end=1980
figsave=True
figformat='pdf'

# Access variables
HID = nc.Dataset(path+filename,'r')
irr_data = HID.variables[product][:]

#Some stuff needed for discrete colorbar
bounds = np.linspace(0,1,11)
norm = colors.BoundaryNorm(boundaries=bounds,ncolors=256)

f_irr = np.concatenate((np.tile(irr_data[0,:,:],(5,1,1)),	\
			np.tile(irr_data[1,:,:],(10,1,1)),	\
			np.tile(irr_data[2,:,:],(10,1,1)),	\
			np.tile(irr_data[3,:,:],(10,1,1)),	\
			np.tile(irr_data[4,:,:],(10,1,1)),	\
			np.tile(irr_data[5,:,:],(10,1,1)),	\
			np.tile(irr_data[6,:,:],(10,1,1)),	\
			np.tile(irr_data[7,:,:],(10,1,1)),	\
			np.tile(irr_data[8,:,:],(8,1,1)),	\
			np.tile(irr_data[9,:,:],(10,1,1)),	\
			np.tile(irr_data[10,:,:],(10,1,1)),\
			np.tile(irr_data[11,:,:],(10,1,1)),\
			np.tile(irr_data[12,:,:],(10,1,1)),\
			np.tile(irr_data[13,:,:],(8,1,1))),axis=0)

irr_plot = np.nanmean(f_irr[yr_start-1900:yr_end-1899,:,:],axis=0)

# Define basemap and some plotting features
m = Basemap(projection='cyl',llcrnrlat=-60,urcrnrlat=90,\
            llcrnrlon=-180,urcrnrlon=180,resolution='c')
# Plot data
plt.figure(figsize=(18,6))
m.drawcoastlines()
m.drawmapboundary(fill_color='lightgray')
m.imshow(irr_plot[60::,:],norm=norm,cmap='Blues',interpolation='none')
c = plt.colorbar()
c.set_label('Irrigated fraction [-]')
plt.title('%i-%i, resolution: 0.5$^\circ$'%(yr_start,yr_end))
plt.show(block=False)

if figsave == True:
  figpath = '/net/exo/landclim/wthiery/observational_analysis/Figures/HID/f_irr.0.5deg.%i-%i.%s'%(yr_start,yr_end,figformat)
  if os.path.exists(figpath):
    print('Figure already exists at path: ' + figpath)
  else:
    print('Saving figure as: ' + figpath)
    plt.savefig(figpath,bbox_inches='tight')
