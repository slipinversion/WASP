import os
import numpy as np
#import matplotlib
#matplotlib.use('Qt5Agg')
from matplotlib import pyplot as plt
from matplotlib import gridspec, ticker, patches
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import cartopy.feature as cf


def set_map_cartopy(ax, margins, tectonic=None, countries=None):
    """
    """
    ax.coastlines(resolution='10m', zorder=3)
    gl = ax.gridlines(linewidth=1, color='black', alpha=0.3, draw_labels=True)
    gl.xlabels_bottom = False
    gl.ylabels_right = False
    if tectonic:
        ax.add_feature(tectonic)
    if countries:
        ax.add_feature(countries)
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
                longitude, latitude, value, zorder=3, vmin=min_val, cmap='jet',
                vmax=max_val, edgecolor='none', transform=transform)
    return ax, cs


if __name__ == '__main__':
    print(1)
