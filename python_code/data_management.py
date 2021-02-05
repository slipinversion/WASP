# -*- coding: utf-8 -*-
"""Module with routines to create json files for the different data sets to be
used in kinematic modelling, from the station and channel metadata stored in
the files with the waveforms-

.. warning::
        
    Currently, only waveform files in sac format are supported. 
"""


import numpy as np
from obspy.io.sac import SACTrace
from obspy import read
import glob
import os
import json
from obspy import read
from scipy.stats import norm
import argparse
import errno
import management as mng
from obspy.taup import TauPyModel


def _dict_trace(file, name_station, channel, azimuth, distance, dt, duration,
               n_start_obs, trace_weight, wavelet_weight, synthetic_trace,
               observed_trace=None, location=None, derivative=False):
    """
    """
    info = {
            'file': file,
            'dt': dt,
            'duration': duration,
            'wavelet_weight': wavelet_weight,
            'start_signal': n_start_obs,
            'trace_weight': trace_weight,
            'synthetic': synthetic_trace,
            'observed': observed_trace,
            'name': name_station,
            'component': channel,
            'azimuth': azimuth,
            'distance': distance,
            'location': location,
            'derivative': derivative
    }
    return info
    

def tele_body_traces(files, tensor_info, data_prop):
    """Write json dictionary with specified properties for teleseismic data
    
    :param files: list of waveform files in sac format
    :param tensor_info: dictionary with moment tensor information
    :param data_prop: dictionary with waveform properties
    :type files: list
    :type tensor_info: dict
    :type data_prop: dict
    
    .. warning::
        
        Make sure the filters of teleseismic data agree with the values in
        sampling_filter.json! 
    """
    if len(files) == 0:
        return
    origin_time = tensor_info['date_origin']
    headers = [SACTrace.read(file) for file in files]
    streams = [read(file) for file in files]
    dt = headers[0].delta
    dt = round(dt, 1)
    n0, n1 = data_prop['wavelet_scales']
    filter0 = data_prop['tele_filter']
    duration = duration_tele_waves(tensor_info, dt)
    wavelet_weight0, wavelet_weight1 = wavelets_body_waves(
            duration, filter0, dt, n0, n1)
    info_traces = []
    event_lat = tensor_info['lat']
    event_lon = tensor_info['lon']
    depth = tensor_info['depth']
    model = TauPyModel(model="ak135f_no_mud")
    
    for file, header, stream in zip(files, headers, streams):
        __failsafe(filter0, header)
        distance, azimuth, back_azimuth = mng._distazbaz(
                header.stla, header.stlo, event_lat, event_lon)
        arrivals = mng.theoretic_arrivals(model, distance / 111.11, depth)
        if header.kcmpnm == 'BHZ':
            arrival = arrivals['p_arrival'][0].time
            weight = 1.0
            wavelet_weight = wavelet_weight0
        elif header.kcmpnm == 'SH':
            arrival = arrivals['s_arrival'][0].time
            weight = 0.5
            wavelet_weight = wavelet_weight1
        starttime = stream[0].stats.starttime
        begin = starttime - origin_time
        n_start_obs = int((arrival - begin) / dt + 0.5)
        info = _dict_trace(
                file, header.kstnm, header.kcmpnm, azimuth, distance / 111.11,
                dt, duration, n_start_obs, weight, wavelet_weight, [],   ##!!!!!!!!!!
                location=[header.stla, header.stlo], derivative=False)
        info_traces.append(info)
    with open('tele_waves.json','w') as f:
         json.dump(
                 info_traces, f, sort_keys=True, indent=4,
                 separators=(',', ': '), ensure_ascii=False)
    return info_traces
    
    
def tele_surf_traces(files, tensor_info, data_prop):
    """Write json dictionary with specified properties for surface wave data
    
    :param files: list of waveform files in sac format
    :param tensor_info: dictionary with moment tensor information
    :param data_prop: dictionary with waveform properties
    :type files: list
    :type tensor_info: dict
    :type data_prop: dict
    """
    if len(files) == 0:
        return
    origin_time = tensor_info['date_origin']
    n0, n1 = data_prop['wavelet_scales']
    wavelet_weight = wavelets_surf_tele(n0, n1)
    info_traces = []
    headers = [SACTrace.read(file) for file in files]
    streams = [read(file) for file in files]
    event_lat = tensor_info['lat']
    event_lon = tensor_info['lon']
    for file, header, stream in zip(files, headers, streams):
        npts = header.npts
        starttime = stream[0].stats.starttime
        begin = origin_time - starttime
        n_start = int(begin / 4.0)
        if n_start < 0:
            continue
        if npts < n_start:
            continue
        if npts + n_start < 0:
            continue
        length = 900
        if 900 >= (npts - n_start):
            length = npts - n_start if n_start > 0 else npts
        distance, azimuth, back_azimuth = mng._distazbaz(
                header.stla, header.stlo, event_lat, event_lon)
        info = _dict_trace(
                file, header.kstnm, header.kcmpnm, azimuth, distance / 111.11,
                4.0, length, n_start, 1.0, wavelet_weight, [],
                location=[header.stla, header.stlo])
        info_traces.append(info)
    with open('surf_waves.json','w') as f:
         json.dump(
                 info_traces, f, sort_keys=True, indent=4,
                 separators=(',', ': '), ensure_ascii=False)
    return info_traces
    
    
def strong_motion_traces(files, tensor_info, data_prop):
    """Write json dictionary with specified properties for strong motion data
    
    :param files: list of waveform files in sac format
    :param tensor_info: dictionary with moment tensor information
    :param data_prop: dictionary with waveform properties
    :type files: list
    :type tensor_info: dict
    :type data_prop: dict
    
    .. warning::
        
        Make sure the filters of strong motion data agree with the values in
        sampling_filter.json!
    """
    if len(files) == 0:
        return
    event_lat = tensor_info['lat']
    event_lon = tensor_info['lon']
    depth = tensor_info['depth']
    origin_time = tensor_info['date_origin']
    headers = [SACTrace.read(file) for file in files]
    dt_strong = headers[0].delta
    dt_strong = round(dt_strong, 1)
    values = [mng._distazbaz(header.stla, header.stlo, event_lat, event_lon)\
        for header in headers]
    distances = [value[0] for value in values]
    zipped = zip(distances, headers)
    arrivals = [np.sqrt(dist**2 + depth**2) / 5 + header.b for dist, header in zipped]
    filter0 = data_prop['strong_filter']
    seismic_moment = tensor_info['moment_mag']
#    outliers = strong_outliers(files, tensor_info)
    if seismic_moment < 2 * 10 ** 26:
        outliers = strong_outliers2(files, distances, tensor_info)
    else:
        outliers = strong_outliers(files, tensor_info)
    duration = duration_strong_motion(
            distances, arrivals, tensor_info, dt_strong)
    n0, n1 = data_prop['wavelet_scales']
    wavelet_weight = wavelets_strong_motion(
            duration, filter0, dt_strong, n0, n1)
    black_list = {
            'PB02': [
                    'HNE',
                    'HNN',
                    'HLE',
                    'HLN'
                    ],
            'PX02': [
                    'HNE',
                    'HLE'
                    ]
    }
    
    info_traces = []
    outlier_traces = []
    headers = [SACTrace.read(file) for file in files]
    streams = [read(file) for file in files]
    weights = [1.0 for st, file in zip(streams, files)]
    weights = [0 if header.kstnm in black_list\
        and header.kcmpnm in black_list[header.kstnm] else weight\
        for weight, header in zip(weights, headers)]
    
    zipped = zip(files, headers, weights, streams, arrivals)
    stat_list = ['KUSD', 'GMLD', 'CESE', 'DIDI', 'YKAV', 'FOCM', 'AKS', 'DATC',
                 'GOMA', 'KRBN']
        
    for file, header, weight, stream, arrival in zipped:
        # if not header.kstnm in stat_list:
        #     continue
        # weight = near_field_weight(tensor_info, header.stla, header.stlo)
        # weight = 0 if header.kstnm in black_list\
        # and header.kcmpnm in black_list[header.kstnm] else weight
        start = origin_time - stream[0].stats.starttime
        distance, azimuth, back_azimuth = mng._distazbaz(
                header.stla, header.stlo, event_lat, event_lon)
        info = _dict_trace(
                file, header.kstnm, header.kcmpnm, azimuth, distance / 111.11,
                dt_strong, duration, int(start / dt_strong), weight,
                wavelet_weight, [], location=[header.stla, header.stlo])
        info_traces.append(info)
    with open('strong_motion_waves.json','w') as f:
         json.dump(
                 info_traces, f, sort_keys=True, indent=4,
                 separators=(',', ': '), ensure_ascii=False)
    with open('outlier_strong_motion_waves.json','w') as f:
         json.dump(
                 outlier_traces, f, sort_keys=True, indent=4,
                 separators=(',', ': '), ensure_ascii=False)
    return info_traces


def cgps_traces(files, tensor_info, data_prop):
    """Write json dictionary with specified properties for cGPS data
    
    :param files: list of waveform files in sac format
    :param tensor_info: dictionary with moment tensor information
    :param data_prop: dictionary with waveform properties
    :type files: list
    :type tensor_info: dict
    :type data_prop: dict
    
    .. warning::
        
        Make sure the filters of cGPS data agree with the values in
        sampling_filter.json!
    """
    if len(files) == 0:
        return
    event_lat = tensor_info['lat']
    event_lon = tensor_info['lon']
    depth = tensor_info['depth']
    origin_time = tensor_info['date_origin']
    headers = [SACTrace.read(file) for file in files]
    dt_cgps = headers[0].delta
    dt_cgps = round(dt_cgps, 1)
    values = [mng._distazbaz(header.stla, header.stlo, event_lat, event_lon)\
        for header in headers]
    distances = [value[0] for value in values]
    zipped = zip(distances, headers)
    arrivals = [np.sqrt(dist**2 + depth**2) / 5 + header.b for dist, header in zipped]
    duration = duration_strong_motion(distances, arrivals, tensor_info, dt_cgps)
    filter0 = data_prop['strong_filter']
    n0, n1 = data_prop['wavelet_scales']
    wavelet_weight = wavelets_strong_motion(
            duration, filter0, dt_cgps, n0, n1, cgps=True)
    info_traces = []
    vertical = ['LXZ', 'LHZ', 'LYZ']
    headers = [SACTrace.read(file) for file in files]
    streams = [read(file) for file in files]
    channels = [header.kcmpnm for header in headers]
    weights = [1.2 if not channel in vertical else 0.6 for file, channel in\
               zip(files, channels)]

    for file, header, stream, weight in zip(files, headers, streams, weights):
        start = origin_time - stream[0].stats.starttime
        distance, azimuth, back_azimuth = mng._distazbaz(
                header.stla, header.stlo, event_lat, event_lon)
        info = _dict_trace(
                file, header.kstnm, header.kcmpnm, azimuth, distance / 111.11,
                dt_cgps, duration, int(start // dt_cgps), weight,
                wavelet_weight, [], location=[header.stla, header.stlo])
        info_traces.append(info)
    with open('cgps_waves.json','w') as f:
         json.dump(
                 info_traces, f, sort_keys=True, indent=4,
                 separators=(',', ': '), ensure_ascii=False)
    return info_traces


def static_data(tensor_info, unit='cm'):
    """Write json dictionary for static GPS data
    
    :param tensor_info: dictionary with moment tensor information
    :param unit: units of static data
    :type tensor_info: dict
    :type unit: string, optional
    """
    info_traces = []
    if not os.path.isfile('cgps_waves.json') and not os.path.isfile('gps_data'):
        return

    if os.path.isfile('cgps_waves.json') and not os.path.isfile('gps_data'):
        cgps_data = json.load(open('cgps_waves.json'))

        names = [file['name'] for file in cgps_data]
        names = list(set(names))
    
        for name in names:
            observed = [0, 0, 0]
            weights = [0, 0, 0]
            error = [0, 0, 0]
            for file in cgps_data:
                name2 = file['name']
                if not name2 == name:
                    continue
                comp = file['component']
                lat, lon = file['location']
                duration = file['duration']
                stream = read(file['file'])
                trace = stream[0].data
                # if len(trace) < 200:
                    # continue
                baseline0 = np.mean(trace[:20])
                baseline1 = np.mean(trace[-25:-5])
                # baseline1 = np.mean(trace[-100:])
                baseline = baseline1 - baseline0
                variance = np.std(trace[-25:-5])
                # variance = np.std(trace[-100:])
                weight = 1.0 if not comp[-1] == 'Z' else 0.5#max(0.1, 2 * norm.cdf(abs(baseline) / variance) - 1)
                index = 0 if comp[-1] == 'Z' else 1 if comp[-1] == 'N' else 2
                observed[index] = str(baseline)
                weights[index] = str(weight)
                error[index] = str(variance)
            info = _dict_trace(
                    [], name, [], [], [], [], 1, 10, weights,
                    [], [], observed_trace=observed, location=[lat, lon])
            info['data_error'] = error
            info_traces.append(info)

    if os.path.isfile('gps_data'):
        if unit == 'cm':
            factor = 1
        elif unit == 'm':
            factor = 100
        elif unit == 'mm':
            factor = 1 / 10
        with open('gps_data', 'r') as infile:
            lines = [line.split() for line in infile]
    
        for line in lines[2:200]:
            name = line[0]
            lon = float(line[1])
            lat = float(line[2])
            observed = [0, 0, 0]
            observed[2] = float(line[3]) * factor
            observed[1] = float(line[4]) * factor
            observed[0] = float(line[5]) * factor
            weights = [0, 0, 0]
            error = [0, 0, 0]
            variance = float(line[6]) * factor
            error[2] = variance
            weights[2] = 1.0#max(0.1, 2 * norm.cdf(abs(observed[2]) / variance) - 1)
            variance = float(line[7]) * factor
            error[1] = variance
            weights[1] = 1.0#max(0.1, 2 * norm.cdf(abs(observed[1]) / variance) - 1)
            variance = float(line[8]) * factor
            error[0] = variance
            weights[0] = 0.5#max(0.1, 2 * norm.cdf(abs(observed[0]) / variance) - 1)
            observed = [str(obs) for obs in observed]
            weights = [str(w) for w in weights]
            info = _dict_trace(
                    [], name, [], [], [], [], 1, 10, weights,
                    [], [], observed_trace=observed, location=[lat, lon])
            info['data_error'] = error
            info_traces.append(info) 
        
    with open('static_data.json','w') as f:
         json.dump(
                 info_traces, f, sort_keys=True, indent=4,
                 separators=(',', ': '), ensure_ascii=False)
    return info_traces


def strong_outliers(files, tensor_info):
    """Routine which seeks to detect data with poorly removed strong motion
    baselines.
    
    :param files: list of waveform files in sac format
    :param tensor_info: dictionary with moment tensor information
    :type files: list
    :type tensor_info: dict
    """
    ratios = []
    inv_ratios = []
    for file in files:
        stream = read(file)
        ratio = np.std(stream[0].data)#lf_ratio(stream, time_shift)
        ratios = ratios + [ratio]
        inv_ratios = inv_ratios + [max(1 / ratio, 10 ** -10)]
#
# modified z_score
#
    outliers = []
    mad_ratios = np.median(abs(ratios - np.median(ratios)))
    cutoff = 6.785# if tensor_info['moment_mag'] >= 8 * 10 ** 25 else 1
    for file, ratio in zip(files, ratios):
        mad = 0.6785 * np.abs(ratio - np.median(ratios)) / mad_ratios
        outliers = outliers + [file] if mad > cutoff else outliers
#    mad_ratios2 = max(
#        np.median(abs(inv_ratios - np.median(inv_ratios))), 10 ** -10)
#    for file, inv_ratio in zip(files, inv_ratios):
#        mad = 0.6785 * np.abs(inv_ratio - np.median(inv_ratios)) / mad_ratios2
#        outliers = outliers + [file] if mad > 6.7 else outliers
    return outliers


def strong_outliers2(files, distances, tensor_info):
    """Routine which seeks to detect data with poorly removed strong motion
    baselines.
    
    :param files: list of waveform files in sac format
    :param distances: list of distances from stations to the event hypocenter
    :param tensor_info: dictionary with moment tensor information
    :type files: list
    :type distances: list
    :type tensor_info: dict
    """
    ratios = []
    inv_ratios = []
    depth = tensor_info['depth']
    for file, dist in zip(files, distances):
        dist2 = np.sqrt(dist ** 2 + depth ** 2)
        stream = read(file)
        ratio = np.std(stream[0].data) * dist2
        ratios = ratios + [ratio]
        inv_ratios = inv_ratios + [max(1 / ratio, 10 ** -10)]
#
# modified z_score
#
    outliers = []
    mad_ratios = max(np.median(abs(ratios - np.median(ratios))), 10 ** -10)
    cutoff = 6.785# if tensor_info['moment_mag'] >= 8 * 10 ** 25 else 1
    for file, ratio in zip(files, ratios):
        mad = 0.6785 * np.abs(ratio - np.median(ratios)) / mad_ratios
        outliers = outliers + [file] if mad > cutoff else outliers
#    mad_ratios2 = max(
#        np.median(abs(inv_ratios - np.median(inv_ratios))), 10 ** -10)
#    for file, inv_ratio in zip(files, inv_ratios):
#        mad = 0.6785 * np.abs(inv_ratio - np.median(inv_ratios)) / mad_ratios2
#        outliers = outliers + [file] if mad > 6.7 else outliers
    return outliers



def __failsafe(filtro, header, cgps=False):
    """
    """
    low_freq = filtro['low_freq']
    high_freq = filtro['high_freq']
    low_freq2 = header.t8
    high_freq2 = header.t9
    if not cgps:
        if not abs(low_freq - low_freq2) + abs(high_freq - high_freq2) < 0.0001:
            print([low_freq, high_freq])
            print([low_freq2, high_freq2])
            raise RuntimeError('Selected filter doesn\'t match filter applied to data')
    else:
        if not abs(high_freq - high_freq2) < 0.0001:
            print([high_freq, high_freq2])
            raise RuntimeError('Selected filter doesn\'t match filter applied to data')


def filling_data_dicts(tensor_info, data_type, data_prop, data_folder):
    """Routine to fill JSON dictionaries containing data properties, for all
    data types selected.
    
    :param tensor_info: dictionary with moment tensor information
    :param data_type: list with data types to use in modelling.
    :param data_prop: dictionary with moment tensor information
    :param data_folder: string with folder where data is located.
    :type tensor_info: dict
    :type data_type: list
    :type data_prop: dict
    :type data_folder: string

        
    .. rubric:: Example:
    
    >>> data_prop = {
            "sampling": {
                "dt_strong": 0.4,
                "dt_tele": 0.2
            },
            "strong_filter": {
                "high_freq": 0.125,
                "low_freq": 0.01
            },
            "tele_filter": {
                "freq0": 0.002,
                "freq3": 1.2,
                "high_freq": 1.0,
                "low_freq": 0.004
            },
            "wavelet_scales": [
                1,
                8
            ]
        }
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
    >>> data_folder = '/path/to/data_folder'
    >>> filling_data_dicts(tensor_info, data_type, data_prop, data_folder)
    """
    folder = os.path.abspath(os.getcwd())
    if 'tele_body' in data_type:
        os.chdir(data_folder)
        tele_traces = get_traces_files('tele_body')
        os.chdir(folder)
        if not os.path.isfile('tele_waves.json'):
            tele_body_traces(tele_traces, tensor_info, data_prop)
    if 'surf_tele' in data_type:
        os.chdir(data_folder)
        surf_traces = get_traces_files('surf_tele')
        os.chdir(folder)
        if not os.path.isfile('surf_waves.json'):
            tele_surf_traces(surf_traces, tensor_info, data_prop)
    if 'strong_motion' in data_type:
        os.chdir(data_folder)
        strong_traces = get_traces_files('strong_motion')
        os.chdir(folder)
        if not os.path.isfile('strong_motion_waves.json'):
            strong_motion_traces(strong_traces, tensor_info, data_prop)
    if 'cgps' in data_type:
        os.chdir(data_folder)
        cgps_data = get_traces_files('cgps')
        os.chdir(folder)
        if not os.path.isfile('cgps_waves.json'):
            cgps_traces(cgps_data, tensor_info, data_prop)
    if 'gps' in data_type:
        static_data(tensor_info, unit='mm')


def get_traces_files(data_type):
    """Get list with waveform files (in sac format) for stations and
    channels selected for modelling
    
    :param data_type: list with data types to use in modelling.
    :type data_type: list
    """
    if data_type == 'tele_body':
        p_traces_files = glob.glob(os.path.join('P', '*tmp'))
        sh_traces_files = glob.glob(os.path.join('SH', '*tmp'))
        traces_files = p_traces_files + sh_traces_files
    if data_type == 'surf_tele':
        p_traces_files = glob.glob(os.path.join('LONG', '*.BHZ.*tmp'))
        sh_traces_files = glob.glob(os.path.join('LONG', '*.SH.*tmp'))
        traces_files = p_traces_files + sh_traces_files
    if data_type == 'strong_motion':
        # traces_files = glob.glob(os.path.join('STR', '*.H[LN][ENZ]*.VEL.tmp'))
        traces_files = glob.glob(os.path.join('STR', '*.VEL.tmp'))
        traces_files = traces_files + glob.glob(os.path.join('STR', '*.BN[ENZ]*.VEL.tmp'))
    if data_type == 'cgps':
        traces_files = glob.glob(os.path.join('cGPS', '*.L[HXY]*tmp'))
    traces_files = [os.path.abspath(file) for file in traces_files]
    return traces_files


def wavelets_body_waves(duration, filtro_tele, dt_tele, n_begin, n_end):
    """Automatic determination of weight of wavelet scales
    
    :param duration: length of signal to be considered for modelling
    :param filtro_tele: filtering properties for body waves
    :param dt_tele: sampling interval for body waves
    :param n_begin: minimum wavelet scale
    :param n_end: maximum wavelet scale
    :type duration: float
    :type filtro_tele: float
    :type dt_tele: float
    :type n_begin: int
    :type n_end: int
    """
    low_freq = filtro_tele['low_freq']
    min_wavelet = int(np.log2(3 * 2**8 * dt_tele * low_freq)) + 1
    min_index = int(duration / dt_tele)
    min_wavelet = max(min_wavelet, 10 - int(np.log2(min_index)))
    p_wavelet_weight = ['0'] + ['1'] * 7
    p_wavelet_weight[:min_wavelet] = ['0'] * min_wavelet
    p_wavelet_weight = ' '.join(p_wavelet_weight[n_begin - 1:n_end])
    p_wavelet_weight = ''.join([p_wavelet_weight, '\n'])
#
# the maximum frequency to be modelled for SH waves is 0.5 Hz
#
    s_wavelet_weight = ['0'] + ['1'] * 6 + ['0'] if dt_tele <= 0.2\
        else ['1'] * 8
    s_wavelet_weight[:min_wavelet] = ['0'] * min_wavelet
    s_wavelet_weight = ' '.join(s_wavelet_weight[n_begin - 1:n_end])
    s_wavelet_weight = ''.join([s_wavelet_weight, '\n'])
    return p_wavelet_weight, s_wavelet_weight
    
    
def wavelets_surf_tele(n_begin, n_end):
    """Automatic determination of weight of wavelet scales

    :param n_begin: minimum wavelet scale
    :param n_end: maximum wavelet scale
    :type n_begin: int
    :type n_end: int
    """
    wavelet_weight = ['0'] * 3 + ['2', '2', '2'] + ['0'] * 2
    wavelet_weight = ' '.join(wavelet_weight[n_begin - 1:n_end])
    wavelet_weight = ''.join([wavelet_weight, '\n'])
    return wavelet_weight


def __wavelets_dart(n_begin, n_end):
    """Automatic determination of weight of wavelet scales
    
    :param n_begin: minimum wavelet scale
    :param n_end: maximum wavelet scale
    :type n_begin: int
    :type n_end: int
    """
    wavelet_weight = ['0'] * 2 + ['2'] * 6# + ['0'] * 2
    wavelet_weight = ' '.join(wavelet_weight[n_begin - 1:n_end])
    wavelet_weight = ''.join([wavelet_weight, '\n'])
    return wavelet_weight


def wavelets_strong_motion(duration, filtro_strong, dt_strong, n_begin, n_end,
                           cgps=False):
    """Automatic determination of weight of wavelet scales
    
    :param duration: length of signal to be considered for modelling
    :param filtro_strong: filtering properties for strong motion data
    :param dt_strong: sampling interval for strong motion data
    :param n_begin: minimum wavelet scale
    :param n_end: maximum wavelet scale
    :param cgps: whether wavelet coefficients are for cGPS or strong motion data
    :type duration: float
    :type filtro_tele: float
    :type dt_tele: float
    :type n_begin: int
    :type n_end: int
    :type cgps: bool, optional
    """
    low_freq = filtro_strong['low_freq']
    high_freq = filtro_strong['high_freq']
    min_wavelet = int(np.log2(3 * 2**8 * dt_strong * low_freq)) + 1
    if cgps:
        min_wavelet = 3#4
    min_index = int(duration / dt_strong)
    min_wavelet = max(min_wavelet, 10 - int(np.log2(min_index)))
    max_wavelet = max(int(np.log2(3 * 2**10 * dt_strong * high_freq)), 1)
#
# largest frequency we can model with wavelets
#
    wavelet_weight = ['2'] * 10
    wavelet_weight[:min_wavelet] = ['0'] * min_wavelet
    wavelet_weight = wavelet_weight[:n_end]
    wavelet_weight = wavelet_weight[:max_wavelet] + ['0'] * (n_end - max_wavelet)
    wavelet_weight = ' '.join(wavelet_weight[n_begin - 1:n_end])
    wavelet_weight = ''.join([wavelet_weight, '\n'])
    return wavelet_weight
    
    
def duration_tele_waves(tensor_info, dt_tele):
    """Source duration plus depth phases
    
    :param tensor_info: dictionary with properties of strong motion data
    :param dt_tele: sampling interval for body waves
    :type tensor_info: dict
    :type dt_tele: float
    """
    depth = tensor_info['depth']
    time_shift = tensor_info['time_shift']
    if depth < 33.0:
        rdur = 2 * (time_shift + depth / 3.5)
    else:
        rdur = 2 * (time_shift + 33.0 / 3.5 + (depth - 33) / 5.0)
    max_length = 180 if depth < 300.0 else 300
    tele_syn_len = min(int(rdur * 1.5 / dt_tele), int(max_length / dt_tele))
    tele_syn_len = max(int(60 / dt_tele), tele_syn_len)
    return tele_syn_len
    
    
def duration_strong_motion(distances, arrivals, tensor_info, dt_strong):
    """Maximum duration estimation for strong motion records
    
    :param distances: distances from stations to hypocenter
    :param arrivals: estimation of first arrival from hypocenter to station
    :param dt_strong: sampling interval for strong motion data
    :param tensor_info: dictionary with moment tensor properties
    :type distances: array
    :type arrivals: array
    :type tensor_info: dict
    :type dt_strong: float
    """
    depth = float(tensor_info['depth'])
    time_shift = float(tensor_info['time_shift'])
    max_dist = max(distances)
    max_len = 350
    if dt_strong == 0.2:
        max_len = 200
    elif dt_strong == 0.4:
        max_len = 350
    elif dt_strong >= 0.8:
        max_len = 600
    max_arrival = max([arrival for arrival in arrivals])
    syn_len = 4 * time_shift + 7 * depth / 50
    syn_len = min(
            int((min(syn_len, max_len) + max_arrival) / dt_strong), 950)
    return syn_len


def duration_dart(distances, arrivals, tensor_info, dt_dart):
    """Maximum duration estimation for DART records
    
    :param distances: distances from stations to hypocenter
    :param arrivals: estimation of first arrival from hypocenter to station
    :param dt_dart: sampling interval for dart data
    :param tensor_info: dictionary with moment tensor properties
    :type distances: array
    :type arrivals: array
    :type tensor_info: dict
    :type dt_dart: float
    """
    depth = float(tensor_info['depth'])
    time_shift = float(tensor_info['time_shift'])
    max_dist = max(distances)
    max_arrival = 0
   # max_arrival = max([arrival for arrival in arrivals])
    max_len = 300
    syn_len = 4 * time_shift + 15 * max_dist / 111.11 + 7 * depth / 50
    syn_len = min(
            int((min(syn_len, max_len) + max_arrival) / dt_dart), 950)
    return syn_len
    
    
def __is_number(string):
    try:
        float(string)
        return True
    except ValueError:
        return False


if __name__ == '__main__':
    import seismic_tensor as tensor

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--folder", default=os.getcwd(),
                        help="folder where there are input files")
    parser.add_argument("-d", "--data_folder", default=os.getcwd(),
                        help="folder with waveform data")
    parser.add_argument("-gcmt", "--gcmt_tensor",
                        help="location of GCMT moment tensor file")
    parser.add_argument("-t", "--tele", action="store_true",
                        help="create JSON for teleseismic body waves")
    parser.add_argument("-su", "--surface", action="store_true",
                        help="create JSON for surface waves")
    parser.add_argument("-st", "--strong", action="store_true",
                        help="create JSON for strong motion data")
    parser.add_argument("--cgps", action="store_true",
                        help="create JSON for cGPS data")
    parser.add_argument("--gps", action="store_true",
                        help="create JSON for static GPS data")
    args = parser.parse_args()
    data_folder = args.data_folder
    if not os.path.isfile('sampling_filter.json'):
        raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT),
                'sampling_filter.json')
    data_prop = json.load(open('sampling_filter.json'))
    os.chdir(args.folder)
    if args.gcmt_tensor:
        cmt_file = args.gcmt_tensor
        tensor_info = tensor.get_tensor(cmt_file=cmt_file)
    else:
        tensor_info = tensor.get_tensor()
    data_type = []
    data_type = data_type + ['gps'] if args.gps else data_type
    data_type = data_type + ['strong_motion'] if args.strong else data_type
    data_type = data_type + ['cgps'] if args.cgps else data_type
    data_type = data_type + ['tele_body'] if args.tele else data_type
    data_type = data_type + ['surf_tele'] if args.surface else data_type
    filling_data_dicts(tensor_info, data_type, data_prop, data_folder)
