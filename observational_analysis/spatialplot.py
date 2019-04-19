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

def spatialplot(m,plot_var,vmin,vmax,figtitle):
    m = Basemap(projection='cyl',llcrnrlat=-60,urcrnrlat=90,llcrnrlon=-180,urcrnrlon=180,resolution='l')
    
    fig = plt.figure(figsize=(12,4))
    sp1 = fig.add_subplot(111)
    m.drawcoastlines(zorder=3)
    m.drawmapboundary(fill_color='lightgray')
    m.imshow(plot_var,cmap='bwr',vmin=vmin,vmax=vmax,interpolation='none',zorder=2)
    c = plt.colorbar()
    c.set_label('$\Delta$T [K]')
    plt.title(figtitle)
    
    plt.show()