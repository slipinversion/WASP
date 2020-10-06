# -*- coding: utf-8 -*-
"""Module with automatic settings for sampling, length and frequency filtering
for teleseismic and near field data
"""


import numpy as np
import json


###################################
# automatic settings
###################################
    
    
def nyquist_frequency(dt):
    """
    """
    return 1 / 2 / dt
    

def sampling(tensor_info, dt_cgps):
    """Automatic setting for data sampling
    
    :param tensor_info: dictionary with moment tensor information
    :type tensor_info: dict
    """
    moment_mag = tensor_info['moment_mag']
    time_shift = tensor_info['time_shift'] 
    depth = tensor_info['depth']
    dt_tele = 0.2    
    if moment_mag > 10 ** 29:
        dt_tele = 0.4
    if 200 < depth < 400:
        dt_tele = 0.4
    elif depth >= 400:
        dt_tele = 0.5
#
# I had some issues modelling deep earthquakes with smaller sampling.
#
    dt_strong = 0.2
    if time_shift <= 40:
        dt_strong = 0.2
    elif time_shift <= 80:
        dt_strong = 0.4
    else:
        dt_strong = 0.8
    dt_strong = dt_strong if depth < 500 else 0.8
    return {'dt_tele': dt_tele, 'dt_strong': dt_strong, 'dt_cgps': dt_cgps}


def filtro_tele(tensor_info):
    """Automatic settings for filter of teleseismic body waves
    
    :param tensor_info: dictionary with moment tensor information
    :type tensor_info: dict
    """
    time_shift = tensor_info['time_shift']
    depth = tensor_info['depth']
    freq0 = 0.003
    freq2 = 1.0
    freq3 = 1.2
    if time_shift < 10:
        freq1 = 0.01
    elif time_shift < 30:
        freq1 = 0.006
    elif time_shift < 80:
        freq0 = 0.002
        freq1 = 0.004
    if time_shift >= 80 or depth >= 200:
        freq0 = 0.001
        freq1 = 0.002
        freq2 = 0.8
        freq3 = 0.9
    filtro = {
            'freq0': freq0, 
            'low_freq': freq1, 
            'high_freq': freq2, 
            'freq3': freq3        
            }
    return filtro
    
    
def filtro_strong(tensor_info):
    """Automatic settings for filter of strong motion data.
    
    :param tensor_info: dictionary with moment tensor information
    :type tensor_info: dict
    """
    time_shift = tensor_info['time_shift']

    if time_shift <= 10:
        min_freq = 0.02
    elif time_shift < 25:
        min_freq = 0.01
    elif time_shift < 50:
        min_freq = 0.005
    elif time_shift < 100:
        min_freq = 0.002
    else:
        min_freq = 0.001
    min_freq = max(min_freq, 1 / 100)
    
    max_freq = 0.125# if tensor_info['time_shift'] > 10 else 0.25
    max_freq = max_freq if tensor_info['depth'] < 300 else 0.05
    filtro = {
            'low_freq': min_freq,
            'high_freq': max_freq
    }
    return filtro


def properties_json(tensor_info, dt_cgps):
    """We set automatic properties for waveform data to be used in modelling
    and write those to a JSON file.
    
    :type tensor_info: dict
    :param tensor_info: dictionary with moment tensor information
    """
    dict1 = sampling(tensor_info, dt_cgps)
    dict2 = filtro_tele(tensor_info)
    dict3 = filtro_strong(tensor_info)
    scales = wavelet_scales(tensor_info)
    dict4 = {
            'sampling': dict1,
            'tele_filter': dict2,
            'strong_filter': dict3,
            'wavelet_scales': scales
    }
    with open('sampling_filter.json', 'w') as f:
         json.dump(
                 dict4, f, sort_keys=True, indent=4, separators=(',', ': '),
                 ensure_ascii=False)
    return dict4
    
    
def wavelet_scales(tensor_info):
    """Scales of wavelet transform to be used in modelling
    
    :param tensor_info: dictionary with moment tensor information
    :type tensor_info: dict
    """
    n_begin = 1
    n_end = 8
#    muestreo = sampling(tensor_info)
#    dt_tele = muestreo['dt_tele']
#    n_end = 7
#    if tensor_info['time_shift'] < 15 or dt_tele >= 0.4:
#        n_end = 8
    return n_begin, n_end
    
    
#def wavelets_body_waves(duration, filtro_tele, dt_tele, n_begin, n_end):
#    """Automatic determination of weight of wavelet scales
#    
#    :param duration: length of signal to be considered for modelling
#    :param filtro_tele: filtering properties for body waves
#    :param dt_tele: sampling interval for body waves
#    :param n_begin: minimum wavelet scale
#    :param n_end: maximum wavelet scale
#    :type duration: float
#    :type filtro_tele: float
#    :type dt_tele: float
#    :type n_begin: int
#    :type n_end: int
#    """
#    low_freq = filtro_tele['low_freq']
#    min_wavelet = int(np.log2(3 * 2**8 * dt_tele * low_freq)) + 1
#    min_index = int(duration / dt_tele)
#    min_wavelet = max(min_wavelet, 10 - int(np.log2(min_index)))
#    p_wavelet_weight = ['0'] + ['1'] * 7
#    p_wavelet_weight[:min_wavelet] = ['0'] * min_wavelet
#    p_wavelet_weight = ' '.join(p_wavelet_weight[n_begin - 1:n_end])
#    p_wavelet_weight = ''.join([p_wavelet_weight, '\n'])
##
## the maximum frequency to be modelled for SH waves is 0.5 Hz
##
#    s_wavelet_weight = ['0'] + ['1'] * 6 + ['0'] if dt_tele <= 0.2\
#        else ['1'] * 8
#    s_wavelet_weight[:min_wavelet] = ['0'] * min_wavelet
#    s_wavelet_weight = ' '.join(s_wavelet_weight[n_begin - 1:n_end])
#    s_wavelet_weight = ''.join([s_wavelet_weight, '\n'])
#    return p_wavelet_weight, s_wavelet_weight
#    
#    
#def wavelets_surf_tele(n_begin, n_end):
#    """Automatic determination of weight of wavelet scales
#
#    :param n_begin: minimum wavelet scale
#    :param n_end: maximum wavelet scale
#    :type n_begin: int
#    :type n_end: int
#    """
#    wavelet_weight = ['0'] * 3 + ['2', '4', '2'] + ['0'] * 2
#    wavelet_weight = ' '.join(wavelet_weight[n_begin - 1:n_end])
#    wavelet_weight = ''.join([wavelet_weight, '\n'])
#    return wavelet_weight
#
#
#def __wavelets_dart(n_begin, n_end):
#    """Automatic determination of weight of wavelet scales
#    
#    :param n_begin: minimum wavelet scale
#    :param n_end: maximum wavelet scale
#    :type n_begin: int
#    :type n_end: int
#    """
#    wavelet_weight = ['0'] * 2 + ['2'] * 6# + ['0'] * 2
#    wavelet_weight = ' '.join(wavelet_weight[n_begin - 1:n_end])
#    wavelet_weight = ''.join([wavelet_weight, '\n'])
#    return wavelet_weight
#
#
#def wavelets_strong_motion(duration, filtro_strong, dt_strong, n_begin, n_end):
#    """Automatic determination of weight of wavelet scales
#    
#    :param duration: length of signal to be considered for modelling
#    :param filtro_strong: filtering properties for strong motion data
#    :param dt_strong: sampling interval for strong motion data
#    :param n_begin: minimum wavelet scale
#    :param n_end: maximum wavelet scale
#    :type duration: float
#    :type filtro_tele: float
#    :type dt_tele: float
#    :type n_begin: int
#    :type n_end: int
#    """
#    low_freq = filtro_strong['low_freq']
#    high_freq = filtro_strong['high_freq']
#    min_wavelet = int(np.log2(3 * 2**8 * dt_strong * low_freq)) + 1
#    if abs(dt_strong - 1.0) < 0.01:
#        min_wavelet = 4
#    min_index = int(duration / dt_strong)
#    min_wavelet = max(min_wavelet, 10 - int(np.log2(min_index)))
#    max_wavelet = max(int(np.log2(3 * 2**10 * dt_strong * high_freq)), 1)
##
## largest frequency we can model with wavelets
##
#    wavelet_weight = ['2'] * 10
#    wavelet_weight[:min_wavelet] = ['0'] * min_wavelet
#    wavelet_weight = wavelet_weight[:n_end]
#    wavelet_weight = wavelet_weight[:max_wavelet] + ['0'] * (n_end - max_wavelet)
#    wavelet_weight = ' '.join(wavelet_weight[n_begin - 1:n_end])
#    wavelet_weight = ''.join([wavelet_weight, '\n'])
#    return wavelet_weight
#    
#    
#def duration_tele_waves(tensor_info, dt_tele):
#    """Source duration plus depth phases
#    
#    :param tensor_info: dictionary with properties of strong motion data
#    :param dt_tele: sampling interval for body waves
#    :type tensor_info: dict
#    :type dt_tele: float
#    """
#    depth = tensor_info['depth']
#    time_shift = tensor_info['time_shift']
#    if depth < 33.0:
#        rdur = 2 * (time_shift + depth / 3.5)
#    else:
#        rdur = 2 * (time_shift + 33.0 / 3.5 + (depth - 33) / 5.0)
#    tele_syn_len = min(int(rdur * 1.5 / dt_tele), 950)
#    tele_syn_len = max(int(60 / dt_tele), tele_syn_len)
#    return tele_syn_len
#    
#    
#def duration_strong_motion(distances, arrivals, tensor_info, dt_strong):
#    """Maximum duration estimation for strong motion records
#    
#    :param distances: distances from stations to hypocenter
#    :param arrivals: estimation of first arrival from hypocenter to station
#    :param dt_strong: sampling interval for strong motion data
#    :param tensor_info: dictionary with moment tensor properties
#    :type distances: array
#    :type arrivals: array
#    :type tensor_info: dict
#    :type dt_strong: float
#    """
#    depth = float(tensor_info['depth'])
#    time_shift = float(tensor_info['time_shift'])
#    max_dist = max(distances)
#    if dt_strong == 0.2:
#        max_len = 200
#    elif dt_strong == 0.4:
#        max_len = 350
#    elif dt_strong >= 0.8:
#        max_len = 600
#    max_arrival = max([arrival for arrival in arrivals])
#    syn_len = 4 * time_shift + 15 * max_dist / 111.11 + 7 * depth / 50
#    syn_len = min(
#            int((min(syn_len, max_len) + max_arrival) / dt_strong), 950)
#    return syn_len
#
#
#def duration_dart(distances, arrivals, tensor_info, dt_dart):
#    """Maximum duration estimation for DART records
#    
#    :param distances: distances from stations to hypocenter
#    :param arrivals: estimation of first arrival from hypocenter to station
#    :param dt_dart: sampling interval for dart data
#    :param tensor_info: dictionary with moment tensor properties
#    :type distances: array
#    :type arrivals: array
#    :type tensor_info: dict
#    :type dt_dart: float
#    """
#    depth = float(tensor_info['depth'])
#    time_shift = float(tensor_info['time_shift'])
#    max_dist = max(distances)
#    max_arrival = 0
##    max_arrival = max([arrival for arrival in arrivals])
#    max_len = 300
#    syn_len = 4 * time_shift + 15 * max_dist / 111.11 + 7 * depth / 50
#    syn_len = min(
#            int((min(syn_len, max_len) + max_arrival) / dt_dart), 950)
#    return syn_len


if __name__ == '__main__':
    import argparse
    import seismic_tensor as tensor
    import os
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--folder", default=os.getcwd(),
                        help="folder where there are input files")
    parser.add_argument("-gcmt", "--gcmt_tensor",
                        help="location of GCMT moment tensor file")
    args = parser.parse_args()
    os.chdir(args.folder)
    if args.gcmt_tensor:
        cmt_file = args.gcmt_tensor
        tensor_info = tensor.get_tensor(cmt_file=cmt_file)
    else:
        tensor_info = tensor.get_tensor() 
    properties_json(tensor_info)
