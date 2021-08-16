import os
import numpy as np
import matplotlib
#matplotlib.use('Qt5Agg')
from matplotlib import pyplot as plt
from matplotlib import gridspec, ticker, patches
import cartopy.crs as ccrs
import cartopy
import cartopy.io.shapereader as shpreader
import cartopy.feature as cf
from matplotlib import cm
from matplotlib.colors import ListedColormap
from glob import glob

"""
Set colorbar for slip
"""
#cD = {'red': [[0.0, 1.0, 1.0],
#  [0.13278008298755187, 1.0, 1.0],
#  [0.34024896265560167, 0.0, 0.0],
#  [0.4979253112033195, 0.0, 0.0],
#  [0.6597510373443983, 1.0, 1.0],
#  [0.8381742738589212, 1.0, 1.0],
#  [1.0, 1.0, 1.0]],
# 'green': [[0.0, 1.0, 1.0],
#  [0.13278008298755187, 1.0, 1.0],
#  [0.34024896265560167, 1.0, 1.0],
#  [0.4979253112033195, 1.0, 1.0],
#  [0.6597510373443983, 1.0, 1.0],
#  [0.8381742738589212, 0.6666666666666666, 0.6666666666666666],
#  [1.0, 0.0, 0.0]],
# 'blue': [[0.0, 1.0, 1.0],
#  [0.13278008298755187, 1.0, 1.0],
#  [0.34024896265560167, 1.0, 1.0],
#  [0.4979253112033195, 0.0, 0.0],
#  [0.6597510373443983, 0.0, 0.0],
#  [0.8381742738589212, 0.0, 0.0],
#  [1.0, 0.0, 0.0]]}
#slipcpt = matplotlib.colors.LinearSegmentedColormap('slipcpt',cD)

rm = 100 # amoount of lines to remove on black end of magma_r
magma_cpt = cm.get_cmap('magma_r', 512) #start with magma_r
white_bit = np.array([255/256, 250/256, 250/256, 1]) #create array of white
slip_cpt = magma_cpt(np.linspace(0,1,512)) #initialize slip_cpt
slip_cpt[rm:,:] = slip_cpt[0:-rm,:] #move beginning up to remove black end
r_s = np.linspace(white_bit[0],slip_cpt[rm][0],rm) # gradient from white to beginning of new magma
g_s = np.linspace(white_bit[1],slip_cpt[rm][1],rm)
b_s = np.linspace(white_bit[2],slip_cpt[rm][2],rm)
slip_cpt[:rm,:][:,0] = r_s
slip_cpt[:rm,:][:,1] = g_s
slip_cpt[:rm,:][:,2] = b_s
slipcpt = ListedColormap(slip_cpt)

def set_map_cartopy(ax, margins, tectonic=None, countries=None, bathymetry=None):
    """
    """
    ax.coastlines(resolution='10m', zorder=3)
    ax.spines["bottom"].set_linewidth(10)
    ax.spines["top"].set_linewidth(10)
    ax.spines["left"].set_linewidth(10)
    ax.spines["right"].set_linewidth(10)
    #ax.stock_img()
    gl = ax.gridlines(linewidth=1, color='black', alpha=0.3, draw_labels=True)
    gl.top_labels = False
    gl.right_labels = False
    #ax.add_feature(cartopy.feature.OCEAN)
    ocean = cartopy.feature.NaturalEarthFeature('physical', 'ocean', \
        scale='110m', edgecolor='none', facecolor=cartopy.feature.COLORS['water'])
#    ax.add_feature(ocean)
    if tectonic:
        ax.add_feature(tectonic)
    if countries:
        ax.add_feature(countries)
    if bathymetry:
        ax.add_feature(bathymetry)
    faults = glob('*.fault')
    if len(faults) > 0:
        for kfault in range(len(faults)):
            fault_trace = np.genfromtxt(faults[kfault])
            ax.plot(fault_trace[:,0],fault_trace[:,1],'k',zorder=100)
    aftershocks = glob('*aftershocks*')
    if len(aftershocks) > 0:
        for kafter in range(len(aftershocks)):
            print('...Adding aftershocks from: '+str(aftershocks[kafter]))
            aftershock = np.genfromtxt(aftershocks[kafter], delimiter="\t")
            aftershock_lat = aftershock[:,3]
            aftershock_lon = aftershock[:,4]
            aftershock_mag = aftershock[:,6]
            ax.scatter(aftershock_lon, aftershock_lat, s=aftershock_mag*10, c='0.5', zorder = 200)

    min_lon, max_lon, min_lat, max_lat = margins
    ax.set_xlim(min_lon, max_lon)
    ax.set_ylim(min_lat, max_lat)
    return ax


def plot_map(ax, latitudes, longitudes, values, min_val=None, max_val=None,
             transform=None):
    """
    """
    min_val = min([np.amin(value) for value in values]) if not min_val else min_val
    max_val = max([np.amax(value) for value in values]) if not max_val else max_val
    zipped = zip(longitudes, latitudes, values)
    for longitude, latitude, value in zipped:
        if np.prod(longitude.shape) > 1:
            cs = ax.pcolormesh(
                longitude, latitude, value, zorder=3, vmin=min_val, cmap=slipcpt,
                vmax=max_val, edgecolor='0.5', lw=0.5, transform=transform)
    return ax, cs


if __name__ == '__main__':
    print(1)
