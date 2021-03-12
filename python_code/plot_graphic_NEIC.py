# -*- coding: utf-8 -*-
"""The routines here allow to plot the solution model FFM modelling, as well
as moment rate function, and waveform fits. 
"""


import argparse
from matplotlib import pyplot as plt
from matplotlib import gridspec, ticker, patches
import matplotlib
from scipy.interpolate import griddata
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import cartopy.feature as cf
from obspy.imaging.beachball import beach
import numpy as np
import os
import get_outputs
import load_ffm_model
import json
import errno
from datetime import datetime
#from clawpack.geoclaw import dtopotools
#
# local modules
#
import fault_plane as pf
import velocity_models as mv
import plane_management as pl_mng
import seismic_tensor as tensor
import shakemap_tools as shakemap
from waveform_plots_NEIC import plot_waveform_fits
from plot_maps_NEIC import plot_map, set_map_cartopy
from matplotlib import cm
from matplotlib.colors import ListedColormap
"""
Set colorbar for slip
"""
rm = 100 # amount of lines to remove on black end of magma_r
ad = 50 # how much at the zero end should be *just* white before transitioning to meet colors
magma_cpt = cm.get_cmap('magma_r', 512) #start with magma_r
white_bit = np.array([255/256, 250/256, 250/256, 1]) #create array of white
slip_cpt = magma_cpt(np.linspace(0,1,512)) #initialize slip_cpt
slip_cpt[rm:,:] = slip_cpt[0:-rm,:] #move beginning up to remove black end
r_s = np.linspace(white_bit[0],slip_cpt[rm][0],rm-ad) # gradient from white to beginning of new magma
g_s = np.linspace(white_bit[1],slip_cpt[rm][1],rm-ad)
b_s = np.linspace(white_bit[2],slip_cpt[rm][2],rm-ad)
slip_cpt[ad:rm,:][:,0] = r_s
slip_cpt[ad:rm,:][:,1] = g_s
slip_cpt[ad:rm,:][:,2] = b_s
slip_cpt[:ad,:] = white_bit
slipcpt = ListedColormap(slip_cpt)

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

def plot_ffm_sol(tensor_info, segments_data, point_sources, shear, solution,
                 vel_model, default_dirs, autosize=False, mr_time=False, evID=None):
    """Main routine. Allows to coordinate execution of different plotting
    routines.
    
    :param tensor_info: dictionary with moment tensor information
    :param segments: list of dictionaries with properties of fault segments
    :param point_sources: properties of point sources of the fault plane
    :param shear: values of shear modulous for each subfault
    :param solution: dictionary with output kinematic model properties
    :param vel_model: dictionary with velocity model properties
    :param default_dirs: dictionary with default directories to be used
    :type default_dirs: dict
    :type tensor_info: dict
    :type segments: list
    :type point_sources: array
    :type shear: array
    :type solution: dict
    :type vel_model: dict
    
    .. rubric:: Example:
   
    First, we load necessary modules.
    
    >>> import json
    >>> import get_outputs # Allows us to get properties of inverted model 
    >>> import management as mng # Allows us to load location of plotting files
    >>> import fault_plane as pf
    >>> import plane_management as pl_mng

    Next, we load necessary data for plots.

    >>> vel_model = json.load(open('velmodel_data.json')) # Assume vel_model stored in file 'velmodel_data.json'
    >>> segments, rise_time, point_sources = pl_mng.__read_planes_info() # Loads point sources and segments information
    >>> solution = get_outputs.read_solution_static_format(segments, point_sources)
    >>> shear = pf.shear_modulous(point_sources, velmodel=vel_model)           
    >>> tensor_info = {
            'moment_mag': 7 * 10 ** 27,
            'date_origin': UTCDateTime(2014, 04, 01, 23, 46, 47)
            'lat': -19.5,
            'lon': -70.5,
            'depth': 25,
            'time_shift': 44,
            'half_duration': 40,
            'centroid_lat': -21,
            'centroid_lon': -70,
            'centroid_depth': 35
        }

    Next, we plot solution

    >>> default_dirs = mng.default_dirs()
    >>> plot_ffm_sol(tensor_info, segments, point_sources, shear, solution,    
    >>>              vel_model, default_dirs)
    
    .. note::
        
        To plot the results of the FFM modelling, we need to run this code 
        in a folder whih contains files Solucion.txt, Fault.time, Fault.pos,
        Event_mult.in, and some among the files synm.tele, synm.str_low,
        synm.str and synm.cgps.
        
    .. note::
    
        When running this code manually, it is good idea to check if
        the information among files Solucion.txt, Fault.pos, Fault.time,
        and Event_mult.in is consistent.
    
    """
    #_plot_vel_model(vel_model, point_sources)
    _plot_moment_rate_function(segments_data, shear, point_sources, mr_time=mr_time)
    #_PlotRiseTime(segments, point_sources, solution)
    #_PlotRuptTime(segments, point_sources, solution)
    _PlotSlipDistribution(segments, point_sources, solution, autosize=autosize)
    _PlotMap(tensor_info, segments, point_sources, solution, default_dirs)
    _PlotSlipTimes(segments, point_sources, solution)

def plot_misfit(used_data_type, forward=False):
    """Plot misfit of observed and synthetic data
    
    :param used_data_type: list with data types used in modelling
    :param forward: whether model is result of kinematic modelling or not
    :type used_data_type: list
    :type forward: bool, optional
    """
    if 'tele_body' in used_data_type:
        if not os.path.isfile('tele_waves.json'):
            raise FileNotFoundError(
                    errno.ENOENT, os.strerror(errno.ENOENT), 'tele_waves.json')
        traces_info = json.load(open('tele_waves.json'))
#        _plot_beachballs(tensor_info, segments, files=traces_info, phase='P')
#        _plot_beachballs(tensor_info, segments, files=traces_info, phase='SH')
        traces_info = get_outputs.get_data_dict(
                traces_info, syn_file='synm.tele')
        values = [['BHZ'], ['SH']]
        for components in values:
            plot_waveform_fits(traces_info, components, 'tele_body',
                               forward=forward)
    if 'surf_tele' in used_data_type:
        if not os.path.isfile('surf_waves.json'):
            raise FileNotFoundError(
                    errno.ENOENT, os.strerror(errno.ENOENT), 'surf_waves.json')
        traces_info = json.load(open('surf_waves.json'))
        traces_info = get_outputs.get_data_dict(
                traces_info, syn_file='synm.str_low')
        values = [['BHZ'], ['SH']]
        for components in values:
            plot_waveform_fits(traces_info, components, 'surf_tele',
                               start_margin='custom', forward=forward)
    if 'strong_motion' in used_data_type:
        if not os.path.isfile('strong_motion_waves.json'):
            raise FileNotFoundError(
                    errno.ENOENT, os.strerror(errno.ENOENT),
                    'strong_motion_waves.json')
        traces_info = json.load(open('strong_motion_waves.json'))
        traces_info = get_outputs.get_data_dict(
                traces_info, syn_file='synm.str')#, observed=False)
        values = [['HLZ', 'HNZ'], ['HLE', 'HNE'], ['HLN', 'HNN']]
        for components in values:
            plot_waveform_fits(traces_info, components, 'strong_motion',
                               forward=forward, start_margin=10)
    if 'cgps' in used_data_type:
        if not os.path.isfile('cgps_waves.json'):
            raise FileNotFoundError(
                    errno.ENOENT, os.strerror(errno.ENOENT), 'cgps_waves.json')
        traces_info = json.load(open('cgps_waves.json'))
        traces_info = get_outputs.get_data_dict(
                traces_info, syn_file='synm.cgps')
        values = [['LXZ', 'LHZ', 'LYZ'], ['LXE', 'LHE', 'LYE'], ['LXN', 'LHN', 'LYN']]
        for components in values:
            plot_waveform_fits(traces_info, components, 'cgps',
                                forward=forward, start_margin=0)
    return


def _plot_vel_model(velmodel, point_sources):
    """We plot the seismic velocity model as a function of depth
    """
    print('Creating Velocity Model Plot...')
    max_depth = [max(ps_segment[:, :, :, :, 2].flatten())\
        for ps_segment in point_sources]
    max_depth = max(max_depth)
    p_vel = np.array(velmodel['p_vel']).astype(np.float)
    sh_vel = np.array(velmodel['s_vel']).astype(np.float)
    thick = np.array(velmodel['thick']).astype(np.float)
    
    depths = np.zeros(len(thick) + 1)
    
    depths[1:] = np.cumsum(thick)
    depths = np.array([depth for depth in depths if depth < 70])
    depths = np.append([depths], [70])#[max_depth])
    plt.plot((p_vel[0], p_vel[0]), (depths[0], depths[1]), 'b-', label='P')
    plt.plot((sh_vel[0], sh_vel[0]), (depths[0], depths[1]), 'r-', label='SH')
    j = len(depths) - 2
    for i in range(j):
        plt.plot((p_vel[i], p_vel[i]), (depths[i], depths[i + 1]), 'b-')
        plt.plot(
            (p_vel[i], p_vel[i + 1]), (depths[i + 1], depths[i + 1]), 'b-')
        plt.plot((sh_vel[i], sh_vel[i]), (depths[i], depths[i + 1]), 'r-')
        plt.plot(
            (sh_vel[i], sh_vel[i + 1]), (depths[i + 1], depths[i + 1]), 'r-')
    
    plt.plot((p_vel[j], p_vel[j]), (depths[j], depths[j + 1]), 'b-')
    plt.plot((sh_vel[j], sh_vel[j]), (depths[j], depths[j + 1]), 'r-')

    plt.title('Crust model for north of Chile')#'Body wave velocity model')
    plt.xlabel('Body wave velocity $(km/s)$')
    plt.ylabel('Depth $(km)$')
    plt.legend(loc='upper right')

    ax = plt.gca()
    ax.invert_yaxis()
    plt.savefig('crust_body_wave_vel_model.png', bbox_inches='tight')
    plt.close()


def _PlotRuptTime(segments, point_sources, solution):
    """We plot time distribution based on the FFM solution model
    """
    print('Creating Rupture Time Plot...')
    rupt_time = solution['rupture_time']
    max_rupt_time = [np.max(rupt_time_seg.flatten()) for rupt_time_seg in rupt_time]
    max_rupt_time = np.max(max_rupt_time)
    slip = solution['slip']
    max_slip = [np.max(slip_seg.flatten()) for slip_seg in slip]
    max_slip = np.max(max_slip)
    x_label = 'Distance along strike $(km)$'
    y_label = 'Distance along dip $(km)$'
    for i_segment, (segment, rupt_time_seg, slip_seg, ps_seg)\
    in enumerate(zip(segments, rupt_time, slip, point_sources)):
#
# Plot the slip distribution
#
        #indexes = np.where(slip_seg < 0.1 * max_slip)
        #rupt_time_seg[indexes] = 0
        fig = plt.figure(figsize=(15, 8))
        ax = fig.add_subplot(111)
        ax.set_ylabel(y_label)
        ax.set_xlabel(x_label)
        fig.subplots_adjust(right=0.75)
        if i_segment == 0:
            ax.plot(0, 0, 'w*', ms=20)
        ax, im = __several_axes(
            rupt_time_seg, segment, ps_seg, ax, max_val=max_rupt_time)
        cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
        cb = fig.colorbar(im, cax=cbar_ax)
        cb.set_label('Rupt_time (s)')
        plt.savefig('RuptTime_plane{}.png'.format(i_segment),
                    bbox_inches='tight')
        plt.close()
    return


def _PlotRiseTime(segments, point_sources, solution):
    """We plot rise time distribution based on the FFM solution model
    """
    print('Creating Rise Time Plot...')
    rise_time = solution['trise']
    max_trise = [np.max(trise_seg.flatten()) for trise_seg in rise_time]
    max_trise = np.max(max_trise)
    fall_time = solution['tfall']
    max_tfall = [np.max(tfall_seg.flatten()) for tfall_seg in fall_time]
    max_tfall = np.max(max_tfall)
    slip = solution['slip']
    max_slip = [np.max(slip_seg.flatten()) for slip_seg in slip]
    max_slip = np.max(max_slip)
    x_label = 'Distance along strike $(km)$'
    y_label = 'Distance along dip $(km)$'
    for i_segment, (segment, trise_seg, tfall_seg, slip_seg, ps_seg)\
    in enumerate(zip(segments, rise_time, fall_time, slip, point_sources)):
#
# Plot the slip distribution
#
        indexes = np.where(slip_seg < 0.1 * max_slip)
        trise_seg[indexes] = 0
        tfall_seg[indexes] = 0
        fig, axes = plt.subplots(1, 2, figsize=(20, 10), sharex=True, sharey=True)
        fig.subplots_adjust(bottom=0.15)
        axes[0].set_ylabel(y_label)
        axes[0].set_xlabel(x_label)
        if i_segment == 0:
            axes[0].plot(0, 0, 'w*', ms=20)
        axes[0], im = __several_axes(
            trise_seg, segment, ps_seg, axes[0], max_val=max_trise,
            autosize=False)
        axes[1].set_ylabel(y_label)
        axes[1].set_xlabel(x_label)
        if i_segment == 0:
            axes[1].plot(0, 0, 'w*', ms=20)
        axes[1], im = __several_axes(
            tfall_seg, segment, ps_seg, axes[1], max_val=max_tfall,
            autosize=False)
        cbar_ax = fig.add_axes([0.1, 0.05, 0.8, 0.05])
        cb = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
        cb.set_label('Rise_time (s)')
        plt.savefig('RiseTime_plane{}.png'.format(i_segment),
                    bbox_inches='tight')
        plt.close()
    return


def _PlotSlipDistribution(segments, point_sources, solution, autosize=False):
    """We plot slip distribution based on the FFM solution model
    """
    print('Creating Slip Distribution Plot...')
    slip = solution['slip']
    rake = solution['rake']
    rupt_time = solution['rupture_time']
    max_slip = [np.max(slip_seg.flatten()) for slip_seg in slip]
    max_slip = np.max(max_slip)
    x_label = 'Distance Along Strike (km)'
    y_label = 'Dist Along Dip (km)'
    for i_segment, (segment, slip_seg, rake_seg, rupttime_seg, ps_seg)\
    in enumerate(zip(segments, slip, rake, rupt_time, point_sources)):
        max_slip_seg = np.max(slip_seg.flatten())
        u = slip_seg * np.cos(rake_seg * np.pi / 180.0) / max_slip_seg
        v = slip_seg * np.sin(rake_seg * np.pi / 180.0) / max_slip_seg
#
# Plot the slip distribution
#
        plt.rc('axes', titlesize=20)
        plt.rc('axes', labelsize=20)
        plt.rc('xtick', labelsize=16)
        plt.rc('ytick',labelsize=16)
        plt.rc('font', size=20)
        fig = plt.figure(figsize=(9, 7))
        ax = fig.add_subplot(111)
        ax.set_ylabel(y_label,fontsize=16)
        ax.set_xlabel(x_label,fontsize=16)
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position('top')
        ax.spines["bottom"].set_linewidth(3)
        ax.spines["top"].set_linewidth(3)
        ax.spines["left"].set_linewidth(3)
        ax.spines["right"].set_linewidth(3)
        fig.subplots_adjust(right=0.85, top=0.85, bottom=0.3)
        stk_subfaults, dip_subfaults, delta_strike, delta_dip, hyp_stk, hyp_dip\
                = pl_mng.__unpack_plane_data(segment)
        x = np.arange(stk_subfaults) * delta_strike - hyp_stk * delta_strike
        y = np.arange(dip_subfaults) * delta_dip - hyp_dip * delta_dip
        X = np.linspace(min(x),max(x),20*len(x))
        Y = np.linspace(min(y),max(y),20*len(y))
        z = (u**2+v**2)**0.5
        xcols, yrows = np.meshgrid(x,y)
        XCOLS, YROWS = np.meshgrid(X,Y)
        orig_grid = np.transpose(np.array([xcols.flatten(),yrows.flatten()]))
        new_grid = np.transpose(np.array([XCOLS.flatten(),YROWS.flatten()]))
        grid_rupttime = griddata(orig_grid, rupttime_seg.flatten(), new_grid, method='linear')
        grid_rupttime_reshape = grid_rupttime.reshape((np.shape(XCOLS)))
        contplot = ax.contour(XCOLS, YROWS, grid_rupttime_reshape, colors='0.75', linestyles='dashed', levels = range(10,500,10), linewidths=1.0)
        grid_z = griddata(orig_grid,z.flatten(),new_grid, method='linear')
        grid_z_reshape = grid_z.reshape((np.shape(XCOLS)))
        ax.quiver(x, y, u, v, scale=15.0, width=0.002, color='0.5', clip_on=False)
        if i_segment == 0:
            ax.plot(0, 0, 'w*', ms=15, markeredgewidth=1.5, markeredgecolor='k')
        ax, im = __several_axes(
                grid_z_reshape, segment, ps_seg, ax, max_val=1., autosize=autosize)
        plt.clabel(contplot, fmt = '%.0f', inline = True, fontsize=14, colors='k')
        cbar_ax = fig.add_axes([0.125, 0.15, 0.5, 0.07])
        sm = plt.cm.ScalarMappable(cmap=slipcpt,norm=plt.Normalize(vmin=0.,vmax=max_slip/100.))
        cb = fig.colorbar(sm, cax=cbar_ax, orientation='horizontal')
        cb.outline.set_linewidth(3)
        cb.set_label('Slip (m)',fontsize=18)
        strk = segment['strike']
        dip = segment['dip']
        ax.text(-0.1, 1.22, 'Strike = '+ str(int(strk)), fontsize=15,
                     fontweight='bold', transform=ax.transAxes, va='top', ha='left')
        ax.text(-0.1, 1.13, 'Dip = ' + str(int(dip)), fontsize=15,
                     fontweight='bold', transform=ax.transAxes, va='top', ha='left')
        ax.text(0, -.04, 'Rupture Front Contours Plotted Every 10 s', fontsize=15,
                     fontweight='bold', transform=ax.transAxes, va='top', ha='left')
        plt.savefig('SlipDist_plane{}.png'.format(i_segment), dpi=300)
        plt.close()
    return

def _PlotSlipTimes(segments, point_sources, solution):
    """We plot slip distribution based on the FFM solution model
    """
    print('Creating Slip Times Plot...')
    slip = solution['slip']
    rake = solution['rake']
    rupt_time = solution['rupture_time']
    rise_time = solution['trise']
    max_slip = [np.max(slip_seg.flatten()) for slip_seg in slip]
    max_slip = np.max(max_slip)
    x_label = 'Distance Along Strike (km)'
    y_label = 'Dist Along Dip (km)'
    for i_segment, (segment, slip_seg, rake_seg, rupttime_seg, risetime_seg, ps_seg)\
    in enumerate(zip(segments, slip, rake, rupt_time, rise_time, point_sources)):
        max_slip_seg = np.max(slip_seg.flatten())
        u = slip_seg * np.cos(rake_seg * np.pi / 180.0) / max_slip_seg
        v = slip_seg * np.sin(rake_seg * np.pi / 180.0) / max_slip_seg
#
# Plot the slip distribution
#
        plt.rc('axes', titlesize=16)
        plt.rc('axes', labelsize=16)
        plt.rc('xtick', labelsize=14)
        plt.rc('ytick',labelsize=14)
        plt.rc('font', size=16)
        fig = plt.figure(figsize=(8, 13))

        ax1 = fig.add_subplot(311)
        ax1.set_ylabel(y_label,fontsize=16)
        ax1.set_xlabel(x_label,fontsize=16)
        ax1.xaxis.tick_top()
        ax1.xaxis.set_label_position('top')
        ax1.spines["bottom"].set_linewidth(3)
        ax1.spines["top"].set_linewidth(3)
        ax1.spines["left"].set_linewidth(3)
        ax1.spines["right"].set_linewidth(3)
        fig.subplots_adjust(right=0.75, hspace=0.38)
        stk_subfaults, dip_subfaults, delta_strike, delta_dip, hyp_stk, hyp_dip\
                = pl_mng.__unpack_plane_data(segment)
        x = np.arange(stk_subfaults) * delta_strike - hyp_stk * delta_strike
        y = np.arange(dip_subfaults) * delta_dip - hyp_dip * delta_dip
        X = np.linspace(min(x),max(x),20*len(x))
        Y = np.linspace(min(y),max(y),20*len(y))
        z = (u**2+v**2)**0.5
        xcols, yrows = np.meshgrid(x,y)
        XCOLS, YROWS = np.meshgrid(X,Y)
        orig_grid = np.transpose(np.array([xcols.flatten(),yrows.flatten()]))
        new_grid = np.transpose(np.array([XCOLS.flatten(),YROWS.flatten()]))
        grid_rupttime = griddata(orig_grid, rupttime_seg.flatten(), new_grid, method='linear')
        grid_rupttime_reshape = grid_rupttime.reshape((np.shape(XCOLS)))
        contplot = ax1.contour(XCOLS, YROWS, grid_rupttime_reshape, colors='0.75', linestyles='dashed', levels = range(10,500,10), linewidths=1.0)
        ax1.quiver(x, y, u, v, scale=15.0, width=0.002, color='0.5', clip_on=False)
        if i_segment == 0:
            ax1.plot(0, 0, 'w*', ms=15, markeredgewidth=1.5, markeredgecolor='k')
       # ax1, im = __several_axes(
       #         grid_rupttime_reshape, segment, ps_seg, ax1, autosize=False, cmap='plasma')
        ax1, im = __several_axes(
                rupttime_seg, segment, ps_seg, ax1, autosize=False, cmap='plasma')
        plt.clabel(contplot, fmt = '%.0f', inline = True, fontsize=11, colors='k')
        cbar_ax1 = fig.add_axes([0.85, 0.678, 0.03, 0.2])
        sm = plt.cm.ScalarMappable(cmap='plasma',norm=plt.Normalize(vmin=0.,vmax=max(grid_rupttime)))
        cb = fig.colorbar(sm, cax=cbar_ax1, orientation='vertical')
        cb.outline.set_linewidth(3)
        cb.set_label('Rupture Time (s)',fontsize=16)
        strk = segment['strike']
        dip = segment['dip']
        ax1.text(-0.1, 1.25, 'Strike = '+ str(int(strk)), fontsize=13,
                     fontweight='bold', transform=ax1.transAxes, va='top', ha='left')
        ax1.text(-0.1, 1.19, 'Dip = ' + str(int(dip)), fontsize=13,
                     fontweight='bold', transform=ax1.transAxes, va='top', ha='left')
        ax1.text(0, -.04, 'Rupture Front Contours Plotted Every 10 s', fontsize=13,
                     fontweight='bold', transform=ax1.transAxes, va='top', ha='left')

        ax2 = fig.add_subplot(312)
        ax2.set_ylabel(y_label,fontsize=16)
        ax2.set_xlabel(x_label,fontsize=16)
        ax2.xaxis.tick_top()
        ax2.xaxis.set_label_position('top')
        ax2.spines["bottom"].set_linewidth(3)
        ax2.spines["top"].set_linewidth(3)
        ax2.spines["left"].set_linewidth(3)
        ax2.spines["right"].set_linewidth(3)
        vrup_ref = segment['rupture_vel']
        
        rupt_vel = rupttime_seg
        idx0 = np.where(rupt_vel<0.1)
        rupt_vel[idx0] = 1. #placeholder
        diag_dist = np.sqrt(xcols**2+yrows**2) 
        rupt_vel = diag_dist/rupt_vel
        rupt_vel[idx0] = vrup_ref #referencevel
        idx10 = np.where(slip_seg>0.05*max_slip)
        mean_vrup = np.mean(rupt_vel[idx10].flatten())
        ## failed attempt to do a subfault-by-subfault vrup calculation. maybe later...
        #rupt_grad_x = np.gradient(rupttime_seg, axis=0)
        #rupt_grad_y = np.gradient(rupttime_seg, axis=1)
        #rupt_grad_full = np.sqrt(rupt_grad_x**2+rupt_grad_y**2)
        #rupt_vel = (np.sqrt((x[1]-x[0])**2+(y[1]-y[0])**2))/rupt_grad_full
        #rupt_vel[idx0] = vrup_ref

        #print('Min Vr: ' +str(min(rupt_vel.flatten()))+' km/s')
        #print('Max Vr: ' +str(max(rupt_vel.flatten()))+' km/s')
        #print('Avg. Vr: '+str(np.mean(rupt_vel.flatten()))+' km/s')
        contplot = ax2.contour(XCOLS, YROWS, grid_rupttime_reshape, colors='0.9', linestyles='dashed', levels = range(10,500,10), linewidths=1.0)
        ax2.quiver(x, y, u, v, scale=15.0, width=0.002, color='0.5', clip_on=False)
        if i_segment == 0:
            ax2.plot(0, 0, 'w*', ms=15, markeredgewidth=1.5, markeredgecolor='k')
        ax2, im = __several_axes(
                rupt_vel, segment, ps_seg, ax2, autosize=False, cmap='viridis_r', min_val=min(rupt_vel.flatten()), max_val=max(rupt_vel.flatten()))
        plt.clabel(contplot, fmt = '%.0f', inline = True, fontsize=11, colors='k')
        sm = plt.cm.ScalarMappable(cmap='viridis_r',norm=plt.Normalize(vmin=0.,vmax=max(rupt_vel.flatten())))
        cbar_ax2 = fig.add_axes([0.85, 0.396, 0.03, 0.2])
        cb = fig.colorbar(sm, cax=cbar_ax2, orientation='vertical')
        cb.outline.set_linewidth(3)
        cb.set_label('Rupture Velocity (km/s)',fontsize=16)
        ax2.text(0, -.04, 'Average* Rupture Velocity to Subfault: '+"{:.2f}".format(mean_vrup)+' km/s',
                 fontsize=13,fontweight='bold', transform=ax2.transAxes, va='top', ha='left')

        ax3 = fig.add_subplot(313)
        ax3.set_ylabel(y_label,fontsize=16)
        ax3.set_xlabel(x_label,fontsize=16)
        ax3.xaxis.tick_top()
        ax3.xaxis.set_label_position('top')
        ax3.spines["bottom"].set_linewidth(3)
        ax3.spines["top"].set_linewidth(3)
        ax3.spines["left"].set_linewidth(3)
        ax3.spines["right"].set_linewidth(3)
        #grid_risetime = griddata(orig_grid, risetime_seg.flatten(), new_grid, method='linear')
        #grid_risetime_reshape = grid_risetime.reshape((np.shape(XCOLS)))
        contplot = ax3.contour(XCOLS, YROWS, grid_rupttime_reshape, colors='0.5', linestyles='dashed', levels = range(10,500,10), linewidths=1.0)
        ax3.quiver(x, y, u, v, scale=15.0, width=0.002, color='0.5', clip_on=False)
        if i_segment == 0:
            ax3.plot(0, 0, 'w*', ms=15, markeredgewidth=1.5, markeredgecolor='k')
       # ax3, im = __several_axes(
       #         grid_risetime_reshape, segment, ps_seg, ax3, autosize=False, cmap='magma_r')
        ax3, im = __several_axes(
                risetime_seg, segment, ps_seg, ax3, autosize=False, cmap='cividis_r')
        plt.clabel(contplot, fmt = '%.0f', inline = True, fontsize=11, colors='k')
        cbar_ax3 = fig.add_axes([0.85, 0.115, 0.03, 0.2]) 
        sm = plt.cm.ScalarMappable(cmap='cividis_r',norm=plt.Normalize(vmin=0.,vmax=max(risetime_seg.flatten())))
        cb = fig.colorbar(sm, cax=cbar_ax3, orientation='vertical')
        cb.outline.set_linewidth(3)
        cb.set_label('Slip Duration (s)',fontsize=16)
        mean_rise = np.mean(risetime_seg[idx10])
        ax3.text(0, -.04, 'Average* Slip Duration (Rise Time): '+"{:.2f}".format(mean_rise)+' s',
                 fontsize=13,fontweight='bold', transform=ax3.transAxes, va='top', ha='left')
        ax3.text(0,-0.12, '*Includes subfaults with >5% maximum slip',
                 fontsize=13,fontweight='bold', transform=ax3.transAxes, va='top', ha='left')
        #plt.gcf().subplots_adjust(bottom=0.3, top=0.85, right=0.75)
        plt.savefig('SlipTime_plane{}.png'.format(i_segment), dpi=300)#,
        #            pad_inches=5,bbox_inches='tight')
        plt.close()
    return

def _PlotCumulativeSlip(segments, point_sources, solution, tensor_info, evID=None):
    print('Creating Along Strike/Dip Plot...')
    slip = solution['slip']
    max_slip = [np.max(slip_seg.flatten()) for slip_seg in slip]
    max_slip = np.max(max_slip)
    x_label = 'Distance Along Strike (km)'
    y_label = 'Dist Along Dip (km)'
    for i_segment, (segment, slip_seg, ps_seg)\
    in enumerate(zip(segments, slip, point_sources)):
        max_slip_seg = np.max(slip_seg.flatten())
#
# Plot the slip distribution
#
        plt.rc('axes', titlesize=11)
        plt.rc('axes', labelsize=11)
        plt.rc('xtick', labelsize=11)
        plt.rc('ytick',labelsize=11)
        plt.rc('font', size=11)
        fig = plt.figure(figsize=(8, 8))


        stk_subfaults, dip_subfaults, delta_strike, delta_dip, hyp_stk, hyp_dip\
                = pl_mng.__unpack_plane_data(segment)
        x = np.arange(stk_subfaults) * delta_strike - hyp_stk * delta_strike
        y = np.arange(dip_subfaults) * delta_dip - hyp_dip * delta_dip

        ### PLOT IT ###
        grid = plt.GridSpec(4, 4, hspace=0.2, wspace=0.2)
        main_ax = fig.add_subplot(grid[1:3, 1:4])
        main_ax.spines["bottom"].set_linewidth(3)
        main_ax.spines["top"].set_linewidth(3)
        main_ax.spines["left"].set_linewidth(3)
        main_ax.spines["right"].set_linewidth(3)
        if i_segment == 0:
            main_ax.plot(0,0,'w*',ms=15, markeredgewidth=1.5, markeredgecolor='k')
        main_ax, im = __several_axes(
                slip_seg, segment, ps_seg, main_ax, max_val=max_slip, autosize=False, depth_ax=False)

        ### Slip_threshold to ignore low-slip region ###
        slip_threshold = 0.1*(max_slip/100.)
        if slip_threshold < 1.0:
            slip_threshold = 1.0
        if max_slip/100. < 3.0:
            slip_threshold = 0.5
        if max_slip/100. < 1.0:
            slip_threshold = 0.2
        slip_seg0 = slip_seg
        cumslip_seg = slip_seg
        idxslip = np.where(cumslip_seg < 100.*slip_threshold)
        cumslip_seg[idxslip] = 0.
        along_strike = np.sum(cumslip_seg, axis = 0)
        along_dip = np.sum(cumslip_seg, axis = 1)
        ### Calculate equivalent slip length in along-strike and along-dip directions ###
        eq_len_AS = shakemap.equivalent_slip_length(along_strike/100, delta_strike)
        left_edge_AS = shakemap.locate_equivalent_slip(along_strike/100, delta_strike, eq_len_AS)
        left_edge_AS = left_edge_AS + x[0]
        eq_len_AD = shakemap.equivalent_slip_length(along_dip/100, delta_dip)
        left_edge_AD = shakemap.locate_equivalent_slip(along_dip/100, delta_dip, eq_len_AD)
        left_edge_AD = left_edge_AD + y[0]

        main_ax.hlines(left_edge_AD, left_edge_AS, left_edge_AS + eq_len_AS, 'k', ls='--')
        main_ax.hlines(left_edge_AD + eq_len_AD, left_edge_AS, left_edge_AS + eq_len_AS, 'k', ls='--')
        main_ax.vlines(left_edge_AS, left_edge_AD, left_edge_AD + eq_len_AD, 'k', ls='--')
        main_ax.vlines(left_edge_AS + eq_len_AS, left_edge_AD, left_edge_AD + eq_len_AD, 'k', ls='--')
        main_ax.set_xlabel(x_label, fontsize=11)
        main_ax.set_ylabel(y_label, fontsize=11)
        main_ax.xaxis.tick_bottom()
        main_ax.xaxis.set_label_position('bottom')
        main_ax.yaxis.tick_right()
        main_ax.yaxis.set_label_position('right')
        main_ax.text(1, 1.05, '*Fault Dimensions Not To Scale', fontsize=9, fontweight='bold',
                     transform=main_ax.transAxes, va='top', ha='right')

        ax1 = fig.add_subplot(grid[0,1:4], sharex=main_ax)
        ax1.xaxis.tick_top()
        ax1.xaxis.set_label_position('top')
        ax1.spines["bottom"].set_linewidth(3)
        ax1.spines["top"].set_linewidth(3)
        ax1.spines["left"].set_linewidth(3)
        ax1.spines["right"].set_linewidth(3)
        ax1.plot(x, along_strike/100,'r', lw=2)
        ax1.scatter(x, along_strike/100, c='0.5', s=25, edgecolor='k', lw=1, zorder=5)
        ax1.set_xlabel('Distance Along Strike (km)', fontsize=11)
        ax1.set_ylabel('Cumulative Slip (m)', fontsize=11) 
        ax1.invert_yaxis()
        ax1.yaxis.tick_right()
        ax1.yaxis.set_label_position('right')
        ax1.set_ylim(0,max(along_strike/100)+.1*max(along_strike/100))
        ax1.vlines(left_edge_AS, 0, max(along_strike/100)+.1*max(along_strike/100), 'k', ls='--')
        ax1.vlines(left_edge_AS + eq_len_AS, 0, max(along_strike/100)+.1*max(along_strike/100), 'k', ls='--')

        ax2 = fig.add_subplot(grid[1:3,0],sharey=main_ax)
        ax2.yaxis.tick_left()
        ax2.xaxis.set_label_position('bottom')
        ax2.yaxis.set_label_position('left')
        ax2.spines["bottom"].set_linewidth(3)
        ax2.spines["top"].set_linewidth(3)
        ax2.spines["left"].set_linewidth(3)
        ax2.spines["right"].set_linewidth(3)
        ax2.plot(along_dip/100, y,'r', lw=2)
        ax2.scatter(along_dip/100, y, c='0.5', s=25, edgecolor='k', lw=1, zorder=5)
        ax2.set_ylabel('Distance Along Dip (km)', fontsize=11)#, rotation=270)
        ax2.set_xlabel('Cumulative Slip (m)', fontsize=11)
        ax2.set_xlim(0,max(along_dip/100)+.1*max(along_dip/100))
        ax2.hlines(left_edge_AD, 0, max(along_dip/100)+.1*max(along_dip/100), 'k', ls='--')
        ax2.hlines(left_edge_AD + eq_len_AD, 0, max(along_dip/100)+.1*max(along_dip/100), 'k', ls='--')
        ax2.invert_xaxis()
        ###LEGEND###
        cbar_ax = fig.add_axes([0.124, 0.712, 0.169, 0.025])
        sm = plt.cm.ScalarMappable(cmap=slipcpt,norm=plt.Normalize(vmin=0.,vmax=max_slip/100.))
        cb = fig.colorbar(sm, cax=cbar_ax, orientation='horizontal')
        cbar_ax.xaxis.set_ticks_position('top')
        cbar_ax.xaxis.set_label_position('top')
        cb.outline.set_linewidth(3)
        cb.set_label('Slip (m)',fontsize=10, fontweight='bold')
        strk = segment['strike']
        dip = segment['dip']
        ax1.text(-0.2, 0.85, 'Strike = '+ str(int(strk)), fontsize=11,
                     transform=ax1.transAxes, va='top', ha='center')
        ax1.text(-0.2, 0.7, 'Dip = ' + str(int(dip)), fontsize=11,
                     transform=ax1.transAxes, va='top', ha='center')
        if evID!=None:
            ax1.text(-0.2, 1.1, evID, fontsize=14, fontweight='bold',
                     transform=ax1.transAxes, va='top', ha='center')

        plt.savefig('CumulativeSlip_plane{}.png'.format(i_segment), dpi=300)
        plt.close()

        ### recover edge lat/lon/dep from _AS and _AD info: ###
        hyp_lon = tensor_info['lon']
        hyp_lat = tensor_info['lat']
        hyp_dep = tensor_info['depth']
        corner_1, corner_2, corner_3, corner_4 = shakemap.translate_xy_to_latlondep(segment, hyp_lon, hyp_lat, hyp_dep, 
                                              eq_len_AS, eq_len_AD, left_edge_AS, left_edge_AD)
    return corner_1, corner_2, corner_3, corner_4

def _PlotSlipDist_Compare(segments, point_sources, input_model,
                          solution):
    """We plot slip distribution based on the FFM solution model
    """
    slip = solution['slip']
    rake = solution['rake']
    slip2 = input_model['slip']
    rake2 = input_model['rake']
    max_slip = [np.max(slip_seg.flatten()) for slip_seg in slip]
    max_slip = np.max(max_slip)
    max_slip2 = [np.max(slip_seg2.flatten()) for slip_seg2 in slip2]
    max_slip2 = np.max(max_slip2)
    max_slip = np.maximum(max_slip, max_slip2)
    x_label = 'Distance along strike $(km)$'
    y_label = 'Distance along dip $(km)$'
    zipped = zip(segments, slip, rake, slip2, rake2, point_sources)
    for i_segment, (segment, slip_seg, rake_seg, slip_seg2, rake_seg2, ps_seg)\
    in enumerate(zipped):
        max_slip_seg = np.max(slip_seg.flatten())
        max_slip_seg2 = np.max(slip_seg2.flatten())
        max_slip_seg = np.maximum(max_slip_seg, max_slip_seg2)
        u = slip_seg * np.cos(rake_seg * np.pi / 180.0) / max_slip_seg
        v = slip_seg * np.sin(rake_seg * np.pi / 180.0) / max_slip_seg
        u2 = slip_seg2 * np.cos(rake_seg2 * np.pi / 180.0) / max_slip_seg
        v2 = slip_seg2 * np.sin(rake_seg2 * np.pi / 180.0) / max_slip_seg
#
# Plot the slip distribution
#
        fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(30, 8))
        ax0.set_ylabel(y_label)
        ax0.set_xlabel(x_label)
        ax1.set_ylabel(y_label)
        ax1.set_xlabel(x_label)
        fig.subplots_adjust(right=0.75)
        stk_subfaults, dip_subfaults, delta_strike, delta_dip, hyp_stk, hyp_dip\
                = pl_mng.__unpack_plane_data(segment)
        x = np.arange(stk_subfaults) * delta_strike - hyp_stk * delta_strike
        y = np.arange(dip_subfaults) * delta_dip - hyp_dip * delta_dip
        ax0.quiver(x, y, u, v, scale=15.0, width=0.003)
        ax1.quiver(x, y, u2, v2, scale=15.0, width=0.003)
        if i_segment == 0:
            ax0.plot(0, 0, 'w*', ms=20)
            ax1.plot(0, 0, 'w*', ms=20)
        ax0, im = __several_axes(slip_seg, segment, ps_seg, ax0,
                                 max_val=max_slip, autosize=False)
        ax1, im = __several_axes(slip_seg2, segment, ps_seg, ax1,
                                 max_val=max_slip, autosize=False)
        ax0.set_title('Inverted model', fontsize=20)
        ax1.set_title('Original model', fontsize=20)
        cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
        cb = fig.colorbar(im, cax=cbar_ax)
        cb.set_label('Slip (cm)')
        plt.savefig('SlipDist_Compare_plane{}.png'.format(i_segment),
                    bbox_inches='tight')
        plt.close()
    return


def _PlotMap(tensor_info, segments, point_sources, solution, default_dirs,
             convex_hulls=[], files_str=None, stations_gps=None,
             option='Solucion.txt', max_slip=None):
    """We plot slip map.
    """
    import collections
    print('Creating Slip Map...')
    plane_info = segments[0]
    stk_subfaults, dip_subfaults, delta_strike, delta_dip, hyp_stk, hyp_dip\
            = pl_mng.__unpack_plane_data(plane_info)
    x = np.arange(stk_subfaults) * delta_strike - hyp_stk * delta_strike
    y = np.arange(dip_subfaults) * delta_dip - hyp_dip * delta_dip
    slip = solution['slip']
#
# accurate plot coordinates
#
    segments_lats, segments_lons, segments_deps = __redefine_lat_lon(segments, point_sources)
    min_lats = [min(segment_lat.flatten()) for segment_lat in segments_lats]
    max_lats = [max(segment_lat.flatten()) for segment_lat in segments_lats]
    min_lons = [min(segment_lon.flatten()) for segment_lon in segments_lons]
    max_lons = [max(segment_lon.flatten()) for segment_lon in segments_lons]
    min_lat = np.min(min_lats)# - 0.5
    max_lat = np.max(max_lats)# + 0.5
    min_lon = np.min(min_lons)# - 0.5
    max_lon = np.max(max_lons)# + 0.5

    # outline the segments in black #
    for segment in range(len(segments_lats)):
        depths = segments_deps[segment].flatten()
        min_lats_idx = np.where(segments_lats[segment].flatten()==min_lats[segment])[0][0]
        cornerA = [segments_lons[segment].flatten()[min_lats_idx], min_lats[segment], depths[min_lats_idx], 0]
        min_lons_idx = np.where(segments_lons[segment].flatten()==min_lons[segment])[0][0]
        cornerB = [min_lons[segment], segments_lats[segment].flatten()[min_lons_idx], depths[min_lons_idx], 1]
        max_lats_idx = np.where(segments_lats[segment].flatten()==max_lats[segment])[0][0]
        cornerC = [segments_lons[segment].flatten()[max_lats_idx], max_lats[segment], depths[max_lats_idx], 2]
        max_lons_idx = np.where(segments_lons[segment].flatten()==max_lons[segment])[0][0]
        cornerD = [max_lons[segment], segments_lats[segment].flatten()[max_lons_idx], depths[max_lons_idx], 3]
        corners = np.c_[cornerA, cornerB, cornerC, cornerD]
        max_dep = max(corners[2,:])
        idx_max = np.where(np.array(corners[2,:]) == max_dep)[0][0]
        min1 = min2 = corners[:,idx_max]
        for i in range(len(corners)):
            if corners[2,i] < min1[2]:
                min2 = min1
                min1 = corners[:,i]
            elif corners[2,i] < min2[2] and collections.Counter(corners[:,i]) != collections.Counter(min1):
                min2 = corners[:,i]

        updip = np.c_[min1,min2]
        corners = np.c_[corners,cornerA]

    margin = 1.3 * (stk_subfaults * delta_strike) / 111.19#min(3 * (n_sub_x * delta_x) / 111.19, 10)
    lat0 = tensor_info['lat']
    lon0 = tensor_info['lon']
 #   margins = [min_lon, max_lon, min_lat, max_lat]
    tectonic = '{}.shp'.format(default_dirs['trench_graphics'])
    dictn = {
            'projection': ccrs.PlateCarree(),
            #'facecolor': '#eafff5'
            'facecolor': 'None'
    }

    fig, ax = plt.subplots(1, 1, figsize=(15, 15), subplot_kw=dictn)
    #ax.outline_patch.set_linewidth(2)
    ax.spines['geo'].set_linewidth(2)
    fig.subplots_adjust(hspace=0, wspace=0, top=0.9, bottom=0.1, right=0.8)
    from cartopy.io.img_tiles import Stamen
    tiler = Stamen('terrain-background')
    ax.add_image(tiler, 10)

#    bathymetry = cf.NaturalEarthFeature('physical','bathymetry','10m')
#    tectonic = cf.ShapelyFeature(
#            shpreader.Reader(tectonic).geometries(), ccrs.PlateCarree(),
#            edgecolor='white', facecolor=(198/255., 236/255., 253/255.), lw=4)
    tectonic = cf.ShapelyFeature(
            shpreader.Reader(tectonic).geometries(), ccrs.PlateCarree(),
            edgecolor='white', facecolor="None", lw=4)
    shpfilename = shpreader.natural_earth(
                resolution='10m', category='cultural', name='admin_0_countries')
#    countries = cf.ShapelyFeature(
#            shpreader.Reader(shpfilename).geometries(), ccrs.PlateCarree(), 
#            edgecolor='black', facecolor='lightgray')
    countries = cf.ShapelyFeature(
            shpreader.Reader(shpfilename).geometries(), ccrs.PlateCarree(), 
            edgecolor='black', facecolor='None')
#    ax = set_map_cartopy(ax, margins, tectonic=tectonic, countries=countries)
#    ax.plot(lon0, lat0, 'w*', markersize=15, transform=ccrs.PlateCarree(),
#            zorder=4)
    if files_str is not None:
        for file in files_str:
            name = file['name']
            latp, lonp = file['location']
            min_lat = min(min_lat, latp)
            max_lat = max(max_lat, latp)
            min_lon = min(min_lon, lonp)
            max_lon = max(max_lon, lonp)
            distance = max(np.abs(latp - lat0), np.abs(lonp - lon0))
            margin = max(margin, 1.2 * distance)
            ax.plot(lonp, latp, 'wo', markersize=10,
                    transform=ccrs.PlateCarree(), zorder=4)
            ax.text(lonp + 0.1, latp + 0.1, '{}'.format(name),
                    transform=ccrs.PlateCarree(), zorder=4)
    if stations_gps is not None:
        max_obs = np.zeros(3)
        stations_gps2 = []
        for name, sta_lat, sta_lon, obs, syn, error in stations_gps:
            min_lat = min(min_lat, sta_lat)
            max_lat = max(max_lat, sta_lat)
            min_lon = min(min_lon, sta_lon)
            max_lon = max(max_lon, sta_lon)
            stations_gps2 = stations_gps2\
                + [[name, sta_lat, sta_lon, obs, syn, error]]
            max_obs = np.maximum([abs(float(v)) for v in obs], max_obs)
            distance = max(np.abs(sta_lat - lat0), np.abs(sta_lon - lon0))
            margin = max(margin, 1.2 * distance)
        max_obs = np.max(max_obs)
        plt.text(lon0 + margin - 2, lat0 + margin - 0.25,
                 '{:.2f} cm'.format(max_obs), transform=ccrs.PlateCarree())
        plt.text(lon0 + margin - 2, lat0 + margin - 0.45,
                 '{:.2f} cm'.format(max_obs), transform=ccrs.PlateCarree())
        for name, sta_lat, sta_lon, obs, syn, error in stations_gps2:
            #plt.plot(sta_lon, sta_lat, 'ks', transform=ccrs.PlateCarree(),
            #         markersize=14)
            gps_z, gps_n, gps_e = syn
            east_west = float(gps_e) / max_obs #/ 100
            north_south = float(gps_n) / max_obs #/ 100
            plt.arrow(sta_lon, sta_lat, east_west, north_south, color='r',
                      zorder=3, linewidth=2, head_width=0.05, head_length=0.05,
                      transform=ccrs.PlateCarree())
            up_down = float(gps_z) / max_obs#/ 100
            #plt.arrow(sta_lon, sta_lat, 0.0, up_down, color='r', zorder=3,
            #          linewidth=2, head_width=0.05, head_length=0.05,
            #          transform=ccrs.PlateCarree())
            gps_z, gps_n, gps_e = obs
            east_west = float(gps_e) / max_obs #/ 100
            north_south = float(gps_n) / max_obs #/ 100
            plt.arrow(sta_lon, sta_lat, east_west, north_south, zorder=3,
                      linewidth=2, head_width=0.05, head_length=0.05,
                      transform=ccrs.PlateCarree())
            up_down = float(gps_z) / max_obs#/ 100
            #plt.arrow(sta_lon, sta_lat, 0.0, up_down, zorder=3,
            #          linewidth=2, head_width=0.05, head_length=0.05,
            #          transform=ccrs.PlateCarree())
            plt.text(sta_lon + 0.02, sta_lat + 0.02, '{}'.format(name),
                     transform=ccrs.PlateCarree())
            err_z, err_n, err_e = error
            width = float(err_e) / max_obs#/ 100
            height = float(err_n) / max_obs#/ 100
        #    ellipse = patches.Ellipse(
        #            (sta_lon + east_west, sta_lat + north_south), width,
        #            height, zorder=4, color='k', linewidth=10,
        #            transform=ccrs.PlateCarree())
        #    plt.gca().add_patch(ellipse)
        #plt.arrow(lon0 + margin - 0.2, lat0 + margin - 0.2, -1, 0, color='r',
        #          zorder=3, linewidth=2, head_width=0.05, head_length=0.05,
        #          transform=ccrs.PlateCarree())
        #plt.arrow(lon0 + margin - 0.2, lat0 + margin - 0.4, -1, 0, color='k',
        #          zorder=3, linewidth=2, head_width=0.05, head_length=0.05,
        #          transform=ccrs.PlateCarree())
    max_slip = max([np.amax(slip_fault) for slip_fault in slip])\
        if not max_slip else max_slip
    margins = [min_lon - 2, max_lon + 2, min_lat - 2, max_lat + 2]
    #ax = set_map_cartopy(ax, margins, tectonic=tectonic, countries=countries)
    ax = set_map_cartopy(ax, margins, tectonic=tectonic, countries=None, bathymetry=None)
    ax.plot(lon0, lat0, '*', markersize=15, transform=ccrs.PlateCarree(),
            zorder=4, markerfacecolor="None", markeredgecolor='k', markeredgewidth=1.5)
    ax.plot(corners[0,:],corners[1,:], 'k',lw=5)
    ax.plot(updip[0,:],updip[1,:], 'r', lw=5)
#
# plot slip map
#
    ax, cs = plot_map(ax, segments_lats, segments_lons, slip, max_val=max_slip, 
                      transform=dictn['projection'])
#    cbar_ax = fig.add_axes([0.17, 0.12, 0.3, 0.02])
    sm = plt.cm.ScalarMappable(cmap=slipcpt,norm=plt.Normalize(vmin=0.,vmax=max_slip/100.))
    #from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    #cb_ax = inset_axes(ax, width = "30%", height = "5%", loc = 'lower left')
    #from matplotlib.axes.Axes import inset_axes
    cb_ax = ax.inset_axes([0.04, 0.04, 0.3, 0.02])
    #cb_ax = fig.add_axes([0.04, 0.04, 0.3, 0.02])
    cbar = plt.colorbar(sm, cax=cb_ax, orientation='horizontal')
#    cbar = plt.colorbar(sm, cax=cbar_ax, orientation='horizontal')
    cbar.outline.set_linewidth(3)
    cbar.set_label('Slip (m)')
    cbar.ax.xaxis.set_ticks_position('top')
    cbar.ax.xaxis.set_label_position('top')
    ax.set_aspect('auto', adjustable=None)
    plt.savefig('Map.png', dpi=300, bbox_inches='tight')
    plt.close()
    return


def _PlotComparisonMap(tensor_info, segments, point_sources, input_model,
                       solution):
    """We plot slip map.
    """
    input_slip = input_model['slip']
    plane_info = segments[0]
    stk_subfaults, dip_subfaults, delta_strike, delta_dip, hyp_stk, hyp_dip\
        = pl_mng.__unpack_plane_data(plane_info)
    if stk_subfaults * dip_subfaults == 1:
        return
    slip = solution['slip']
#
# accurate plot coordinates
#
    segments_lats, segments_lons, segments_deps = __redefine_lat_lon(segments, point_sources)

    margin = min(1.5 * (stk_subfaults * delta_strike) / 111.19, 10)#3
    lat0 = tensor_info['lat']
    lon0 = tensor_info['lon']
    margins = [lon0 - margin, lon0 + margin, lat0 - margin, lat0 + margin]
    dictn = {'projection': ccrs.PlateCarree(), 'facecolor': '#eafff5'}
    shpfilename = shpreader.natural_earth(
        resolution='10m', category='cultural', name='admin_0_countries')
    countries = cf.ShapelyFeature(
        shpreader.Reader(shpfilename).geometries(), ccrs.PlateCarree(), 
        edgecolor='black', facecolor='lightgray')

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(30, 15), subplot_kw=dictn)
    ax1.set_title('Inverted model', fontsize=22)
    ax2.set_title('Original model', fontsize=22)
    fig.subplots_adjust(hspace=0, wspace=0.1, top=0.9, bottom=0.3)
    for ax in [ax1, ax2]:
        ax = set_map_cartopy(ax, margins, countries=countries)
        ax.plot(lon0, lat0, 'w*', markersize=15, transform=ccrs.PlateCarree(),
            zorder=4)
    max_slip = max([np.amax(slip_fault) for slip_fault in slip])
    max_slip2 = max([np.amax(input_slip2) for input_slip2 in input_slip])
    max_slip = max(max_slip, max_slip2)
#
# plot slip map
#
    ax1, cs1 = plot_map(ax1, segments_lats, segments_lons, slip, max_val=max_slip,
                       transform=dictn['projection'])
    ax2, cs2 = plot_map(ax2, segments_lats, segments_lons, input_slip,
                        max_val=max_slip, transform=dictn['projection'])
    fig.subplots_adjust(bottom=0.15)
    cbar_ax = fig.add_axes([0.1, 0.05, 0.8, 0.05])
    cb = fig.colorbar(cs2, cax=cbar_ax, orientation='horizontal')
    cb.set_label('Slip (cm)')
    plt.savefig('Comparison.png', bbox_inches='tight')
    plt.close()
    return


def __redefine_lat_lon(segments, point_sources):
    """
    """
    segments_lats = [[]] * len(segments)
    segments_lons = [[]] * len(segments)
    segments_deps = [[]] * len(segments)
    for i, point_sources_seg in enumerate(point_sources):
        lat = point_sources_seg[:, :, :, :, 0]
        lon = point_sources_seg[:, :, :, :, 1]
        dep = point_sources_seg[:, :, :, :, 2]
        ny, nx, a, b = lat.shape
        new_lat = np.zeros((ny + 1, nx + 1))
        new_lon = np.zeros((ny + 1, nx + 1))
        new_dep = np.zeros((ny + 1, nx + 1))
        for j in range(ny):
            for k in range(nx):
                new_lat[j, k] = lat[j, k, 0, 0]
                new_lon[j, k] = lon[j, k, 0, 0]
                new_dep[j, k] = dep[j, k, 0, 0]
        for k in range(nx):
            new_lat[-1, k] = lat[-1, k, -1, 0]
            new_lon[-1, k] = lon[-1, k, -1, 0]
            new_dep[-1, k] = dep[-1, k, -1, 0]
        for j in range(ny):
            new_lat[j, -1] = lat[j, -1, 0, -1]
            new_lon[j, -1] = lon[j, -1, 0, -1]
            new_dep[j, -1] = dep[j, -1, 0, -1]
        new_lat[-1, -1] = lat[-1, -1, -1, -1]
        new_lon[-1, -1] = lon[-1, -1, -1, -1]
        new_dep[-1, -1] = dep[-1, -1, -1, -1]
        segments_lats[i] = new_lat
        segments_lons[i] = new_lon
        segments_deps[i] = new_dep
    return segments_lats, segments_lons, segments_deps


def _plot_moment_rate_function(segments_data, shear, point_sources, mr_time=None):
    """We plot moment rate function
    """
    print('Creating Moment Rate Plot...')
    segments = segments_data['segments']
    plane_info = segments[0]
    dt = 0.01
    model = load_ffm_model.load_ffm_model(segments_data, point_sources)
    slip = model['slip']
    trup = model['trup']
    tl = model['trise']
    tr = model['tfall']
    properties = pl_mng.__unpack_plane_data(plane_info)
    delta_strike, delta_dip = [properties[2], properties[3]]
    stk_subfaults, dip_subfaults = [properties[0], properties[1]]
    tmax = 1.5 * np.max([np.amax(trup_seg + tl_seg + tr_seg)\
                         for trup_seg, tl_seg, tr_seg in zip(trup, tl, tr)])
    tmax = tmax if stk_subfaults * dip_subfaults > 1\
        else (tl[0][0, 0] + tr[0][0, 0]) * 8.0
    nmax = int((tmax/dt + 1))
    mr = np.zeros(nmax)
    seismic_moment = 0
    
    for segment, slip_seg, trup_seg, trise_seg, tfall_seg, shear_seg,\
    point_sources_seg in zip(segments, slip, trup, tl, tr, shear, point_sources):
        ny_total, nx_total = np.shape(slip_seg)
        moment_rate = np.zeros(nmax)
        for iy in range(ny_total):
            for ix in range(nx_total):
                rupt_vel = segment['rupture_vel']
                rise_time = np.zeros(nmax)
                tfall = tfall_seg[iy, ix]
                trise = trise_seg[iy, ix]
                array1 = np.arange(0, trise, dt)
                tmid = len(array1)
                rise_time[:tmid]\
                    = (1 - np.cos(np.pi * array1 / trise)) / (trise + tfall)
                array2 = np.arange(0, tfall, dt)
                tend = tmid + len(array2)
                rise_time[tmid:tend]\
                    = (1 + np.cos(np.pi * array2 / tfall)) / (tfall + trise)
                duration = int(max(delta_strike, delta_dip) / dt / rupt_vel)
                source_dur = np.ones(duration) / duration
                start_index = max(0, int(trup_seg[iy, ix] / dt))
#                source_dur[start_index:start_index + duration] = np.ones(duration)
                product = slip_seg[iy, ix] * shear_seg[iy, ix] / 100 / 10
                sub_rise_time = rise_time[:tend]
                convolve = np.convolve(source_dur, sub_rise_time)
                moment_rate[start_index:start_index + len(convolve)]\
                = moment_rate[start_index:start_index + len(convolve)]\
                + convolve * product

        seismic_moment = seismic_moment\
            + np.sum((slip_seg / 100) * (shear_seg / 10)\
                     * (delta_strike * 1000) * (delta_dip * 1000))
        
#    
# find moment rate function
#
        for i in range(nmax):
            time = i * dt
            mr[i] = mr[i]\
                + moment_rate[i] * (delta_strike * 1000) * (delta_dip * 1000)
    
    time = np.arange(nmax) * dt 
    with open('STF.txt', 'w') as outf:
        outf.write('dt: {}\n'.format(dt))
        outf.write('Time[s]     Moment_Rate [Nm]\n')
        for t, val in zip(time, mr):
            outf.write('{:8.2f}:   {:8.4e}\n'.format(t, val))
    
    seismic_moment = np.trapz(mr, dx=0.01)
    magnitude = 2.0 * (np.log10(seismic_moment * 10 ** 7) - 16.1) / 3.0
    fig = plt.figure(figsize=(10,8))
    plt.rc('axes', titlesize=16)
    plt.rc('axes', labelsize=16)
    plt.rc('xtick', labelsize=14)
    plt.rc('ytick',labelsize=14)
    plt.rc('font', size=16)
    ax = fig.add_subplot(111)
    ax.spines["bottom"].set_linewidth(3)
    ax.spines["top"].set_linewidth(3)
    ax.spines["left"].set_linewidth(3)
    ax.spines["right"].set_linewidth(3)

    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Relative Moment Rate (Nm/s)')
    rel_mr = mr/(max(mr))
    plt.text(
        0.99 * max(time), 0.95 * max(rel_mr),
        'Max Mr: {:.2E} Nm/sec'.format(max(mr)), ha='right',fontweight='bold')    
    plt.text(
        0.99 * max(time), 0.90 * max(rel_mr),
        'M$_0$: {:.2E} Nm'.format(seismic_moment), ha='right',fontweight='bold')
    plt.text(
        0.99 * max(time), 0.85 * max(rel_mr), 'M$_w$: {:.2f}'.format(magnitude), ha='right',fontweight='bold')
    plt.grid('on', ls='dotted')
    #plt.minorticks_on()
    plt.grid(which='minor',linestyle='dotted', color='0.5')
    plt.grid(which='major',linestyle='dotted', color='0.5')
    ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.2))
    ax.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(0.1))
    ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(40))
    ax.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(20))
    ax.fill_between(time, rel_mr, color='0.9')
    ax.plot(time, rel_mr, 'k', lw=2)
    if mr_time == None:
        tenth = np.ones(len(time))*0.1
        idx = np.argwhere(np.diff(np.sign(rel_mr-tenth))).flatten()
        ax.vlines(time[idx][-1],0,1,'r',linestyle='dashed',lw=2, dashes=(0,(5,3)))
    else:
        ax.vlines(mr_time,0,1,'r',linestyle='dashed',lw=2,dashes=(9,(5,3)))
    ax.set_ylim([0,1])
    ax.set_xlim([0,max(time)])
    plt.savefig('MomentRate.png')#, bbox_inches='tight')
    plt.close()
    return


def _PlotSnapshotSlip(tensor_info, segments, point_sources, solution):
    """we plot snapshots of the rupture process.
    """
    plane_info = segments[0]
    dt = 0.01
    slip = solution['slip']
    trup = solution['rupture_time']
    tl = solution['t_rise']
    tr = solution['t_fall']
    stk_subfaults, dip_subfaults, delta_strike, delta_dip, hyp_stk, hyp_dip\
        = pl_mng.__unpack_plane_data(plane_info)
    if [stk_subfaults, dip_subfaults] == [1, 1]:
        return
    tmid = trup + tl
    tstop = trup + tl + tr
    srmax = slip / (tr + tl) * 2.
#
# Define the vector field of the slip and plotting parameters.
#   
    x = np.arange(stk_subfaults) * delta_strike - hyp_stk * delta_strike
    y = np.arange(dip_subfaults) * delta_dip - hyp_dip * delta_dip
    ratio = max(1,
        (9 / 16) * ((stk_subfaults * delta_strike) / (dip_subfaults * delta_dip))) 
    vmax = np.amax(slip)
    tmax = np.amax(trup)
    step = int((tmax/dt + 1) / 9.)
#
# snapshots of rupture process
#
    fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(9, 9))
    i = 1

    for ax in axes.flat:
        time = i * step * dt
        srate, cslip, broken = __rupture_process(
            time, slip, srmax, trup, tl, tr, tmid, tstop)
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        if np.max(broken) > np.min(broken): 
            ax.contour(-broken, 1, colors='k',
                       extent=__extent_plot(plane_info))
        ax.contourf(cslip, cmap='jet', vmin=0, vmax=vmax,
                    extent=__extent_plot(plane_info))
        ax.plot(0, 0, 'r*', ms=20)
        ax.invert_yaxis()
        ax.set_aspect(ratio)
        ax.set_title('Time: {0:.2f} s'.format(time))
        if i == 9: im = ax.contourf(x, y, cslip, cmap='jet')
        i = i + 1
            
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    cb = fig.colorbar(im, cax=cbar_ax)
    cb.set_label('Slip (cm)')
    fig.text(0.1, 0.1, 'CSN Automatic \nSolution', fontsize=50, color='gray',
             ha='left', va='bottom', alpha=0.5, wrap=True)
    plt.savefig('SlipSnapshot.png', bbox_inches='tight')
    plt.close()
    return

def shakemap_polygon(segments, point_sources, solution, tensor_info, evID):
    corner_1, corner_2, corner_3, corner_4 = _PlotCumulativeSlip(segments, point_sources, solution, tensor_info, evID=evID)
    # tranlate above output into lon/lat/dep for txt output #
    print('Writing shakemap_polygon.txt to file...')
    ### Create shakemap_polygon.txt file ###
    now = datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S")
    corners = []
    txtout = open("shakemap_polygon.txt",'w')
    txtout.write('#Source: USGS NEIC Rapid Finite Fault \n')
    txtout.write('#Event ID: '+evID +'\n')
    txtout.write('#Model created: '+now+'\n')
    txtout.write(corner_1+'\n')
    txtout.write(corner_2+'\n')
    txtout.write(corner_3+'\n')
    txtout.write(corner_4+'\n')
    txtout.write(corner_1+'\n')
    txtout.close()

def plot_beachball(tensor_info, segments, files=None, phase=None):
    """Here we plot the beachball for the event. Optionally, we add the
    location of teleseismic data used in the FFM modelling.
    """
#
# Get the focal mechanism
#
    plane_info = segments[0]
    strike = plane_info['strike']
    dip = plane_info['dip']
    rake = plane_info['rake']
    fm = [strike, dip, rake]
#
# Plot the beach ball
#
    fig = plt.figure(figsize=(7, 7))
    bb = beach(fm, width=3.0, zorder=1)
    ax = plt.gca()
    ax.add_collection(bb)
    plt.plot(0, 0, 'b*', markersize=20)
    ax.set_aspect('equal')
#
# Load the station informations
#
    if files:
        plt.plot(0, 0, 'r*', markersize=20)
        for file in files:
            comp = file['component']
            if comp == 'BHZ': comp = 'P'
            dist = file['distance']
            az = file['azimuth']
            nam = file['name']
            if comp == phase:
                r = dist * (3.0 / 2.) / 90.
                x1 = np.sin(az * np.pi / 180.) * r
                y1 = np.cos(az * np.pi / 180.) * r
                plt.plot(x1, y1, 'ro')
                plt.text(x1 + 0.01, y1 + 0.01, nam)

    fig.patch.set_visible(False)
    plt.gca().axis('off')
    name_plot = '{}_azimuthcover.png'.format(phase) if phase else 'Tensor.png'
    plt.savefig(name_plot)
    plt.close()
    return
    
    
def _plot_waveforms(files, components, type_str, tensor_info,
                    start_margin=10, test=False, forward=False):
    """We plot the observed and synthetic data for a set of stations and
    channels.
    """
    files = [file for file in files if file['component'] in components]
    files = sorted(files, key=lambda k: k['azimuth'])
    azimuth = [file['azimuth'] for file in files]
    fig = plt.figure(figsize=(13, 9))
    numrows_phase = len(files) // 4 + 1
    gs = gridspec.GridSpec(max(4, numrows_phase), 4)
    for file in files:
        dt = file['dt']
        nstart = file['start_signal']
        margin = int(start_margin / dt) if nstart > int(start_margin / dt) else 0 
        obs = np.array(file['observed'])
        if nstart >= 0:
            obs = np.concatenate((np.zeros(nstart), obs))
        if type_str in ['tele_body', 'strong_motion', 'cgps']:
            obs = obs[nstart - margin:]
        if type_str == 'cgps':
            obs = obs - obs[10]
        if type_str == 'surf_tele':
            if nstart >= 0: obs = obs[nstart:]
            else: obs = np.concatenate((np.zeros(-nstart), obs))

        obs = obs if not forward else 0 * obs
        syn = np.array(file['synthetic'])
        syn = np.concatenate((np.zeros(margin), syn))
        length = min(len(obs), len(syn), file['duration'])
        length = min(length, int(950 / dt))\
            if not type_str in ['surf_tele'] else length
        obs = np.array(obs[:length])
        syn = np.array(syn[:length])
        dt = file['dt']
        az = file['azimuth']
        dist = file['distance']
        name = file['name']
        time = np.arange(-margin, length - margin) * dt\
            if not type_str=='dart' else np.arange(-margin, length - margin)
        jj = azimuth.index(az)
        weight = file['trace_weight']
        alpha = 1 if weight > 0 else 0.1
        ax = fig.add_subplot(gs[jj % numrows_phase, jj // numrows_phase])
        ax.plot(time, obs, 'k', linewidth=0.8, alpha=alpha)
        ax.plot(time, syn, 'r', linewidth=0.8, alpha=alpha)
#        if 'BHZ' in components or 'SH' in components:
        min_val = min(np.min(obs), np.min(syn))
        max_val = max(np.max(obs), np.max(syn))
        ax.vlines(0, min_val, max_val)
        ax.text(
            0.1, 0.9, '{:0.1f}'.format(az), ha='center',
            va='center', transform=ax.transAxes)
        ax.text(
            0.1, 0.1, '{:0.1f}'.format(dist), ha='center',
            va='center', transform=ax.transAxes)
        ax.text(
            0.9, 0.9, name, ha='center', va='center', transform=ax.transAxes)
        ax.set_xlim([np.min(time), np.max(time)])        
        ax.set_ylim([min_val, max_val])
        ax.xaxis.set_major_locator(
            ticker.MaxNLocator(nbins=3, min_n_ticks=3))
        ax.yaxis.set_major_locator(
            ticker.MaxNLocator(nbins=3, min_n_ticks=3))

    if type_str == 'cgps':
        if 'LXZ' in components: plot_name = 'LXZ_cgps_waves.png'
        if 'LXN' in components: plot_name = 'LXN_cgps_waves.png'
        if 'LXE' in components: plot_name = 'LXE_cgps_waves.png'

    if type_str == 'strong_motion':
        if 'HNZ' in components: plot_name = 'HNZ_strong_motion_waves.png'
        if 'HNN' in components: plot_name = 'HNN_strong_motion_waves.png'
        if 'HNE' in components: plot_name = 'HNE_strong_motion_waves.png'

    if type_str == 'tele_body':
        if 'P' in components: plot_name = 'P_body_waves.png'
        if 'SH' in components: plot_name = 'SH_body_waves.png'

    if type_str == 'surf_tele':
        if 'P' in components: plot_name = 'Rayleigh_surf_waves.png'
        if 'SH' in components: plot_name = 'Love_surf_waves.png'
    
    if type_str == 'dart':
        plot_name = 'Dart_waves.png'

    plt.savefig(plot_name, bbox_inches='tight')
    plt.close()
    return


#def _plot_runup(field=None):
#    """
#    """
#    inputfile = os.path.join(os.getcwd(),'Solucion.txt')
#    resolucion = 60 * 2 # Cantidad de puntos por grado (60 = 1 minuto resolucion)
#
#    #########################################
#    # Cambiar de formato los archivos
#
#    with open(inputfile) as f:
#        content = f.readlines()
#    content = [x.strip() for x in content]
#    print(content[1].split())
#    delta_x, delta_y = [val for val in content[1].split() if 'km' in val]
#    delta_x = delta_x.split('=')[1] if '=' in delta_x else delta_x
#    delta_x = delta_x[:-2]
#    print(delta_x)
#    delta_y = delta_y.split('=')[1] if '=' in delta_y else delta_y
#    delta_y = delta_y[:-2]
#    print(delta_y)
#    dx = float(delta_x) * np.ones((len(content[10:-1])))
#    dy = float(delta_y) * np.ones((len(content[10:-1])))
#    f.close()
#
#    input_units = {"length":"km", "width":"km", "depth":"km", "slip":"cm"}
#
#    data = np.loadtxt(inputfile, skiprows=10)
#    lat = data[:,0]
#    lon = data[:,1]
#    depth = data[:,2]
#    slip = data[:,3]
#    rake = data[:,4]
#    strike = data[:,5]
#    dip = data[:,6]
#
#    dir=os.getcwd()
#
#    ##look for runup field measure
#    import glob
#    q=glob.glob('r*.txt')
#    print(q)
#
#    field = None
#    if len(q)==1:
#        field=np.loadtxt(q[0])
#        field = field[:,[0,2]]    
#
#    InputFile = "Temp.csv"
#
#    hdr = "Longitude,Latitude,Depth(km),Length,Width,Strike,Dip,Rake,Slip"
#    np.savetxt(
#        InputFile,
#        list(zip(lon, lat, depth, dx, dy, strike, dip, rake, slip, slip)),
#        delimiter=',', header=hdr, fmt="%.5f", comments="")
#
#    fault = dtopotools.CSVFault()
#    fault.read(InputFile, input_units=input_units)
#
#    # Mostrar los valores relevantes
#    print("The seismic moment is %g N-m" % fault.Mo())
#    print("The Moment magnitude is %g" % fault.Mw())
#    print("  (Assuming the rigidity mu of all subfaults is the default value %g Pa)"\
#	  % fault.subfaults[0].mu)
#
#
#    # Cear archivo dtopo necesario para GeoClaw
#    xlower = min(lon) - 1
#    xupper = max(lon) + 1
#    ylower = min(lat) - 1
#    yupper = max(lat) + 1
#    dx = 1. / resolucion
#    mx = int((xupper - xlower)/dx + 1)
#    xupper = xlower + (mx-1)*dx
#    my = int((yupper - ylower)/dx + 1)
#    yupper = ylower + (my-1)*dx
#    x = np.linspace(xlower,xupper,mx)
#    y = np.linspace(ylower,yupper,my)
#
#    dtopo = fault.create_dtopography(x,y,times=[1.], verbose=True)
#    os.remove(InputFile)
#
#    ### INPUT ###############################
#
#    outputfile = os.path.join(os.getcwd(), "runup_output.dat")
#
#    ### Cargar variables ####################
#
#    lon = dtopo.x
#    lat = dtopo.y
#
#    deformation = dtopo.dZ[0,:,:]
#    max_def = np.amax(deformation,axis=1)
#
#    ### Calculo del runup ###################
#    H = np.amax(max_def)
#    d = 4*1000.
#    alpha = np.sqrt(1 + H/d)
#    beta_1 = 1./180.*2*np.pi
#    beta_2 = 2./180.*2*np.pi
#    theta = np.pi/2
#    f = max_def/H
#    f[f < 0] = 0
#    ro_1 = 2.831 * H * np.power(H/d, 1./4.)\
#        * np.sqrt(alpha * np.power(np.tan(beta_1), -1))
#    ro_2 = 2.831 * H * np.power(H/d, 1./4.)\
#        * np.sqrt(alpha * np.power(np.tan(beta_2), -1))
#    runup_1 =  f * ro_1 * np.sqrt(np.sin(theta))
#    runup_2 =  f * ro_2 * np.sqrt(np.sin(theta))
#
#    ### Graficar deformacion y runup calculado
#
#    d = {'figsize':(13,6)}
#    fig, (ax0, ax1, ax2) = plt.subplots(1, 3, sharey=True, **d)
#
#    fault.plot_subfaults(slip_color=True, axes=ax0)
#
#    cax1 = ax1.imshow(
#        deformation, extent=[lon[1], lon[-1], lat[1], lat[-1]],
#        origin='lower', cmap="seismic", vmin=-max_def.max(), vmax=max_def.max())
#    plt.colorbar(cax1, ax=ax1, orientation='vertical')
#    ax1.set_title('Seafloor Deformation [m]')
#    max_lev = round(max_def.max())
#    levels = np.linspace(-max_lev, max_lev, 4*int(max_lev) + 1)
#    ax1.contour(deformation, levels, linewidths=0.5, colors='black',
#                extent=[lon[1], lon[-1], lat[1], lat[-1]])
#
#    ax2.plot(runup_1, lat, label=r'$\beta=1$')
#    ax2.plot(runup_2, lat, label=r'$\beta=2$')
#    ax2.set_title('Runup analitico [m]')
#    ax2.legend()
#    ax2.grid(linestyle='--')
#
#    if field is not None:
#        for l,r in zip(field[:,0],field[:,1]):
##            ax2.plot([0,r], [l, l], color='k', linewidth=2, alpha=0.3)
#            ax2.scatter(r, l, alpha=0.5)
#
#    plt.autoscale(enable=True, tight=True)
#
#    ### Guardar datos
#    print(outputfile)
#    np.savetxt(outputfile, np.transpose([lat, runup_1]), fmt='%.3f')
#
#    fig.savefig(os.path.join(dir,'analitic_runup.png'),bbox_inches='tight')
#    plt.close()
#    return


def __extent_plot(plane_info):
    """
    """
    stk_subfaults, dip_subfaults, delta_strike, delta_dip, hyp_stk, hyp_dip\
        = pl_mng.__unpack_plane_data(plane_info)
    return [-hyp_stk * delta_strike, (stk_subfaults - hyp_stk) * delta_strike,
            -hyp_dip * delta_dip, (dip_subfaults - hyp_dip) * delta_dip]

def __several_axes(data, segment, point_source_seg, ax, min_val=None, max_val=None,
                   autosize=True, cmap=slipcpt, depth_ax=True):
    """
    """
    stk_subfaults, dip_subfaults, delta_strike, delta_dip, hyp_stk, hyp_dip\
        = pl_mng.__unpack_plane_data(segment)
    min_dist = - (hyp_dip + 0.5) * delta_dip
    max_dist = (dip_subfaults - hyp_dip - 0.5) * delta_dip
    min_strike = - (hyp_stk + 0.5) * delta_strike
    max_strike = (stk_subfaults - hyp_stk - 0.5) * delta_strike
    dep = point_source_seg[:, :, :, :, 2]
    min_depth = dep[0, 0, 0, 0]
    max_depth = dep[-1, 0, -1, 0]
    if stk_subfaults * dip_subfaults == 1:
        dip = segment['dip']
        delta_z = delta_dip * np.sin(dip * np.pi / 180.0) / 2
        min_depth = min_depth - delta_z
        max_depth = max_depth + delta_z
    max_val = np.max(data.flatten()) if not max_val else max_val
    min_val = np.min(data.flatten()) if not min_val else min_val
    im = ax.imshow(data, cmap=cmap, origin='lower', vmin=min_val, vmax=max_val,
                   aspect='auto',
                   extent=[min_strike, max_strike, min_dist, max_dist])
    ax.set(adjustable='datalim')
    if autosize:
        ax.figure.set_size_inches(4 * stk_subfaults * delta_strike / dip_subfaults / delta_dip, 4)
    if depth_ax:
        ax2 = ax.twinx()
        ax2.set_xlim([min_strike, max_strike])
        ax2.set_ylim([min_depth, max_depth])
        ax2.set(adjustable='datalim')
        ax2.set_ylabel('Depth (km)',fontsize=16)
        ax2.invert_yaxis()
        if autosize:
            ax2.figure.set_size_inches(4 * stk_subfaults * delta_strike / dip_subfaults / delta_dip, 4)
    ax.invert_yaxis()
    return ax, im


def __rupture_process(time, slip, trup, tmid, tstop, rise_time,
                      point_sources, dt):
    """We give slip rate, rupture front, and accumulated slip at a certain
    time ``time``.
    """
    ny_total, nx_total = np.shape(slip)
    srate = np.zeros((ny_total, nx_total))
    cslip = np.zeros((ny_total, nx_total))
    broken = np.ones((ny_total, nx_total))
    
    for i in range(ny_total):
        for j in range(nx_total):
            convolve = rise_time[i, j, :]
            index = int(time / dt - trup[i, j] / dt)
            if (time < trup[i, j]):
                broken[i, j] = 0.
            elif index < len(convolve):
                srate[i, j] = convolve[index] * slip[i, j]#srmax[iy, ix]
            if trup[i, j] < time <= tmid[i, j]:
                cslip[i, j] = (time - trup[i, j]) * srate[i, j] / 2.
            if tmid[i, j] < time <= tstop[i, j]:
                cslip[i, j]\
                    = slip[i, j] - (tstop[i, j] - time) * srate[i, j] / 2.
            if (time > tstop[i, j]):
                cslip[i, j] = slip[i, j]
    return srate, cslip, broken
    
    
def __add_watermark(fig):
    """
    """
    fig.text(0.1, 0.1, 'CSN Automatic \nSolution', fontsize=50, color='gray',
             ha='left', va='bottom', alpha=0.5, wrap=True)
    return fig
                
                
if __name__ == '__main__':
    """
    """
    import management as mng
    import eventpage_downloads as dwnlds
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--folder", default=os.getcwd(),
                        help="folder where there are input files")
    parser.add_argument("-gcmt", "--gcmt_tensor",
                        help="location of GCMT moment tensor file")
    parser.add_argument("-ffms", "--ffm_solution", action="store_true",
                        help="plot FFM solution slip maps, rise time")
    parser.add_argument("-t", "--tele", action="store_true",
                        help="plot misfit of teleseismic data")
    parser.add_argument("-su", "--surface", action="store_true",
                        help="plot misfit of surface waves data")
    parser.add_argument("-st", "--strong", action="store_true",
                        help="plot strong motion stations and strong motion misfit")
    parser.add_argument("--cgps", action="store_true",
                        help="plot misfit of cGPS data")
    parser.add_argument("--gps", action="store_true", help="plot GPS data")
    parser.add_argument("-o","--option",choices=['autoscale','noscale'],
                        help="choose whether Rupture plot needs to be scaled")
    parser.add_argument("-mr","--mrtime", default='0', type=float,
                        help="choose cutoff time for Moment Rate plot")
    parser.add_argument("-shakemap","--shakemappolygon",action="store_true",
                        help="create shakemap_polygon.txt")
    parser.add_argument("-ev","--EventID", nargs='?', const='Not Provided', type=str,
                        help="Provide event ID")
    parser.add_argument("-d","--downloads",action="store_true",
                        help="create event page download files")
    args = parser.parse_args()
    os.chdir(args.folder)
    used_data = []
    used_data = used_data + ['strong_motion'] if args.strong else used_data
    used_data = used_data + ['cgps'] if args.cgps else used_data
    used_data = used_data + ['tele_body'] if args.tele else used_data
    used_data = used_data + ['surf_tele'] if args.surface else used_data
    default_dirs = mng.default_dirs()
    if args.gcmt_tensor:
        cmt_file = args.gcmt_tensor
        tensor_info = tensor.get_tensor(cmt_file=cmt_file)
    else:
        tensor_info = tensor.get_tensor()
    segments_data = json.load(open('segments_data.json'))
    segments = segments_data['segments']
    rise_time = segments_data['rise_time']
    point_sources = pf.point_sources_param(segments, tensor_info, rise_time)
    solution = get_outputs.read_solution_static_format(segments)
    if args.ffm_solution:
        if not os.path.isfile('velmodel_data.json'):
            vel_model = mv.select_velmodel(tensor_info, default_dirs)
        else:
            vel_model = json.load(open('velmodel_data.json'))
        shear = pf.shear_modulous(point_sources, velmodel=vel_model)
        if args.option == 'autoscale':
            autosize = True
        else:
            autosize = False
        if args.mrtime > 0:
            mr_time = args.mrtime
        else:
            mr_time = None
        if args.EventID:
            evID = args.EventID
        else:
            evID = None
        plot_ffm_sol(tensor_info, segments_data, point_sources, shear, solution,
                     vel_model, default_dirs, autosize=autosize, mr_time=mr_time, evID=evID)
        
    if args.shakemappolygon:
        if args.EventID:
            evID = args.EventID
        else:
            evID = None
        shakemap_polygon(segments, point_sources, solution, tensor_info, evID=evID)

    traces_info, stations_gps = [None, None]
    if args.gps:
        names, lats, lons, observed, synthetic, error\
                = get_outputs.retrieve_gps()
        stations_gps = zip(names, lats, lons, observed, synthetic, error)
    if args.strong:
        traces_info = json.load(open('strong_motion_waves.json'))
    #if args.strong or args.gps:
    if args.gps:
        solution = get_outputs.read_solution_static_format(segments)
        _PlotMap(tensor_info, segments, point_sources, solution, default_dirs,
                 files_str=traces_info, stations_gps=stations_gps)
        #input_model = load_ffm_model.load_ffm_model(
        #        segments_data, point_sources, option='Fault.time')
        #_PlotSlipDist_Compare(segments, point_sources, input_model, solution)
        #_PlotComparisonMap(tensor_info, segments, point_sources, input_model,
        #                   solution)
    if args.downloads:
       dwnlds.write_CMTSOLUTION_file()
       if args.EventID:
            evID = args.EventID
       else:
            evID = None
       dwnlds.write_Coulomb_file(eventID=evID)
       dwnlds.write_Okada_displacements()

    plot_misfit(used_data)#, forward=True)
