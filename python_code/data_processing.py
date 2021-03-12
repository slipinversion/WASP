#!/usr/bin/python
# -*- coding: utf-8 -*-
"""Module with routines to perform data processing (bandpass filtering,
downsampling, removal of instrumental response) to use such data in finite
fault modelling
"""

import subprocess
from shutil import copy2
from multiprocessing import Pool#, cpu_count
from obspy import read
from obspy.io.sac import SACTrace
from obspy.taup import TauPyModel
from obspy.core.utcdatetime import UTCDateTime
from obspy.geodetics import locations2degrees, kilometers2degrees
import numpy as np
import os
import glob
import logging
import management as mng
import modulo_logs as ml
import wang_baseline_removal_v1 as wang1
import time
import seismic_tensor as tensor
from scipy.integrate import cumtrapz
from scipy.signal import butter, filtfilt
import errno
import json
import remove_response as response


###################################
# teleseismic body waves
###################################


def select_process_tele_body(tele_files0, tensor_info, data_prop):
    """Module for automatic selection and processing of teleseismic body waves.
    
    :param tensor_info: dictionary with moment tensor information
    :param data_prop: dictionary with waveform properties
    :param tele_files0: files with body waves to be selected and processed
    :type tensor_info: dict
    :type data_prop: dict
    :type tele_files0: list
    .. rubric:: Example:
    
    >>> from obspy.core.utcdatetime import UTCDateTime
    >>> import glob
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
            'centroid_depth': 35,
            'timedelta': 30 * 60
        }
    >>> tele_files0 = glob.glob('*sac')
    >>> select_process_tele_body(tele_files0, tensor_info, data_prop)
    .. note::
        
        Currently, this code allows only to process files in sac format, where the 
        instrumental response is in a SACPZ file. Other data formats are not
        supported.
    """
    time1 = time.time()
    depth = tensor_info['depth']
    timedelta = tensor_info['timedelta']
    tele_files1 = __pre_select_tele(tele_files0, tensor_info, timedelta)
    if not tele_files1:
        return
    if not os.path.isdir('logs'):
        os.mkdir('logs')
    if not os.path.isdir('P'):
        os.mkdir('P')
    if not os.path.isdir('SH'):
        os.mkdir('SH')
    logger1 = ml.create_log('body_tele_processing',
                            os.path.join('logs', 'body_tele_processing.log'))
    logger1 = ml.add_console_handler(logger1)
    logger1.info('Start selection of teleseismic body data')

    data_folder = os.path.abspath(os.getcwd())
    
    logger1.info('Get theoretical arrivals of P and S waves')
    model = TauPyModel(model="ak135f_no_mud")
    __picker(tensor_info, tele_files1, depth, model)
    logger1.info('time spent getting theoretic arrivals: {}'.format(
        time.time() - time1))
    logger1.info('Remove response for selected teleseismic traces')
    response_files = glob.glob('SACPZ*') + glob.glob('SAC_PZs*')
    __remove_response_body(tele_files1, response_files, tensor_info, data_prop,
                           logger=logger1)
    logger1.info('time until remove response: {}'.format(time.time() - time1))
    logger1.info('Rotate horizontal components body waves')
#
    os.chdir(os.path.join(data_folder, 'SH'))
    horizontal_sacs = glob.glob('*DIS')
    __rotation(horizontal_sacs, logger=logger1)
    
    logger1.info('Process body waves')
    os.chdir(data_folder)
    os.chdir(os.path.join(data_folder, 'P'))
    p_files = glob.glob('*BHZ*DIS')
    process_body_waves(p_files, 'P', data_prop, tensor_info)
    os.chdir(os.path.join(data_folder, 'SH'))
    s_files = glob.glob('*SH*DIS')
    process_body_waves(s_files, 'SH', data_prop, tensor_info)
    
    logger1.info('Body waves have been succesfully processed')
    os.chdir(data_folder)
    logger1.info('time until create files: {}'.format(time.time() - time1))
    files = glob.glob(os.path.join('P', '*BHZ*DIS.tmp0'))
    select_tele_stations(files, 'P', tensor_info)
    files = glob.glob(os.path.join('SH', '*SH*DIS.tmp0'))
    select_tele_stations(files, 'SH', tensor_info)
    logger1.info('total time spent: {}'.format(time.time() - time1))
    ml.close_log(logger1)
    return
    

def __pre_select_tele(tele_files0, tensor_info, timedelta):
    """Module which pre-selects teleseismic waveforms for processing
    """
    date_origin = tensor_info['date_origin']
    sacheaders = [SACTrace.read(sac) for sac in tele_files0]
    tele_sacs = [
        sac for sac, sacheader in zip(tele_files0, sacheaders) if
        (sacheader.stla and sacheader.stlo)]
    sacheaders = [SACTrace.read(sac) for sac in tele_sacs]
    
    for sacheader, sac in zip(sacheaders, tele_sacs):
        sacheader.reftime = UTCDateTime(date_origin)
        sacheader.write(sac, byteorder='little')
    
    used_tele_sacs = []
    for sac, sacheader in zip(tele_sacs, sacheaders):
        if not __select_distance(tensor_info, sacheader.stla, sacheader.stlo,
                                 max_distance=89, min_distance=31):
            continue
        if not __select_time(sac, tensor_info):
            continue
        if not sacheader.kcmpnm in sac:
            continue
        if not sacheader.khole in ['00', '', None]:
            continue
        used_tele_sacs = used_tele_sacs + [sac]
    return used_tele_sacs


def __picker(tensor_info, used_tele_sacs, depth, model):
    """Obspy taupy picker.
    """
    event_lat = tensor_info['lat']
    event_lon = tensor_info['lon']
    sacheaders = [SACTrace.read(sac) for sac in used_tele_sacs]
    for sacheader, sac in zip(sacheaders, used_tele_sacs):
        dist = locations2degrees(
                event_lat, event_lon, sacheader.stla, sacheader.stlo)
        arrivals = mng.theoretic_arrivals(model, dist, depth)
        sacheader.t1 = arrivals['p_arrival'][0].time
        sacheader.t2 = arrivals['pp_arrival'][0].time
        sacheader.t3 = arrivals['p_slowness']
        sacheader.t5 = arrivals['s_arrival'][0].time
        sacheader.t6 = arrivals['ss_arrival'][0].time
        sacheader.t7 = arrivals['s_slowness']
        sacheader.evla = event_lat
        sacheader.evlo = event_lon
        sacheader.write(sac, byteorder = 'little')
    return


def __remove_response_body(used_tele_sacs, response_files, tensor_info,
                           data_prop, logger=None):
    """We remove instrumental response for data, and replace it for a common
    instrumental response
    """
    string = [
            'ZEROS             2\n',
            '0.               0.\n',
            '0.               0.\n',
            'POLES             4\n',
            '-6.17E-03  6.17E-03\n',
            '-6.17E-03 -6.17E-03\n',
            '-39.18        49.12\n',
            '-39.18       -49.12\n',
            'CONSTANT         3948'
    ]
    string = ''.join(string)
    with open('response_file', 'w') as resp_file:
        resp_file.write(string)
    filtro = data_prop['tele_filter']
    depth = tensor_info['depth']
    time_shift = tensor_info['time_shift']
    t_min = 120.0
    t_max = 400.0
    if depth < 300:
        t_max = 240.0 if time_shift < 80 else 400.0 
    freq0 = filtro['freq0']
    freq1 = filtro['low_freq']
    freq2 = filtro['high_freq']
    freq3 = filtro['freq3']

    sacheaders = [SACTrace.read(sac) for sac in used_tele_sacs]
    full_signal1 = lambda arrival, begin, end:\
        arrival - t_min >= begin and arrival + t_max <= end
    new_name = lambda station, channel, network:\
        os.path.join('{}.{}.{}.DIS'.format(station, channel, network))

    sac_files = []
    old_responses = []
    names = []

    for sacheader, sac in zip(sacheaders, used_tele_sacs):
        channel = sacheader.kcmpnm
        phase = 'P' if channel == 'BHZ' else 'SH'
        arrival = sacheader.t1 if channel == 'BHZ' else sacheader.t5
        begin = sacheader.b
        end = sacheader.e
        full_signal = full_signal1(arrival, begin, end)
        name = sacheader.kstnm
        network = sacheader.knetwk
        str_name = os.path.join(phase, new_name(name, channel, network))
        loc_code = sacheader.khole# if sacheader.khole else '--'
        # pz_files = [response for response in response_files\
        #             if name in response and channel in response\
        #             and loc_code in response]
        pz_files0 = [response for response in response_files\
                     if name in response and channel in response]
        if loc_code == '00':
            pz_files = [response for response in pz_files0 if '00' in response]
        if loc_code in [None, '']:
            pz_files = [response for response in pz_files0 if '--' in response]
            pz_files = pz_files +\
                [response for response in pz_files0 if '__' in response]
        if not pz_files:
            continue
        pzfile_val = next(iter(pz_files))#
        if not os.path.isfile(pzfile_val):
            continue

        if full_signal:
            sac_files = sac_files + [sac]
            old_responses = old_responses + [pzfile_val]
            names = names + [str_name]

    pre_processing = 'rtr \ntaper \nhp c {} n 2 p 2 \n'.format(freq1)
    out, err = response.replace_response_sac(
        sac_files, old_responses, names, new_response='response_file',
        pre_processing=pre_processing, add_response=True, freq0=freq0,
        freq1=freq1, freq2=freq2, freq3=freq3)

    if logger:
        logger.debug(out.decode('utf-8'))
        if err: logger.warning(err.decode('utf-8'))
    return


def __rotation(used_tele_sacs, logger=None):
    """We rotate horizontal body waves to match transverse and radial
    directions to the event.
    """
    input_sac = ''
    string = lambda x, y, z, w: '{}.{}.{}.{}'.format(x, y, z, w)
    
    headers = [SACTrace.read(sac) for sac in used_tele_sacs]
    stations = [header.kstnm for header in headers]
    stations = list(set(stations))

    for sta in stations:
        file1 = glob.glob('{}.BH[N1]*DIS'.format(sta))
        file2 = glob.glob('{}.BH[E2]*DIS'.format(sta))
        if file1 and file2:
            file1 = file1[0]
            head1 = SACTrace.read(file1)
            file2 = file2[0]
            head2 = SACTrace.read(file2)
            sta = head1.kstnm
            net = head1.knetwk
            transverse = string(sta, 'SH', net, 'DIS')
            radial = string(sta, 'RD', net, 'DIS')
            cmpaz1 = head1.cmpaz
            cmpaz2 = head2.cmpaz
            if np.remainder(abs(cmpaz1 - cmpaz2), 90.0) == 0:
                east0 = head2.b
                east1 = head2.e
                north0 = head1.b
                north1 = head1.e
                begin = min(east0, north0) + 1
                end = max(east1, north1) - 1
                input_sac = '{}\ncut {} {} \n'.format(input_sac, begin, end)\
                            + 'r {} {} \n'.format(file1, file2)\
                            + 'rot \n write {} {}'.format(radial, transverse)

    input_sac = '{}\nquit\n'.format(input_sac)
    out, err = mng.run_sac(input_sac)
    
    if logger:
        logger.debug(out.decode('utf-8'))
        if err: logger.warning(err.decode('utf-8'))
    return


def process_body_waves(tele_files, phase, data_prop, tensor_info):
    r"""We high pass filter and downsample teleseismic data. 
    
    :param phase: string indicating whether phase is P or SH
    :param tensor_info: dictionary with moment tensor information
    :param data_prop: dictionary with waveform properties
    :param tele_files: body wave files to be resampled and filtered
    :type phase: string
    :type tensor_info: dict
    :type data_prop: dict
    :type tele_files: list
    """
    dt = data_prop['sampling']['dt_tele']
    date_origin = tensor_info['date_origin']
    filtro = data_prop['tele_filter']
    low_freq = filtro['low_freq']
    high_freq = filtro['high_freq']
    t_min = 120.0
    t_max = 360.0
    
    sacheaders = [SACTrace.read(sac) for sac in tele_files]       
    streams = [read(sac) for sac in tele_files]
    for sac, sacheader in zip(tele_files, sacheaders):
        sacheader.t8 = low_freq
        sacheader.t9 = high_freq
        if phase == 'SH':
            sacheader.kcmpnm = 'SH'
        sacheader.write(sac, byteorder = 'little')
    if phase == 'SH':
        for sac, sacheader in zip(tele_files, sacheaders):
            sacheader.kcmpnm = 'SH'
            sacheader.write(sac, byteorder = 'little')
    
    sacheaders = [SACTrace.read(sac) for sac in tele_files]
    streams = [read(sac) for sac in tele_files]
    high_freq = 1 if phase == 'P' else 0.5
    high_freq = min(high_freq, 1 / 2.5 / dt)
    for sacheader, stream, sac in zip(sacheaders, streams, tele_files):
        arrival = sacheader.t1 if phase == 'P' else sacheader.t5
        start = date_origin + arrival - t_min
        end = date_origin + arrival + t_max
        tr = stream[0]
        tr.trim(starttime=start, endtime=end)
        nyq = 0.5 / sacheader.delta
        high_freq2 = high_freq / nyq
        b, a = butter(2, high_freq2, btype='lowpass')
        filt_data = filtfilt(b, a, tr.data)
        tr.data = filt_data
        factor = dt * tr.stats.sampling_rate
        if abs(int(round(factor)) - factor) > 10 ** - 2: continue
        tr.decimate(int(round(factor)), no_filter=True)
        tr.data = tr.data - tr.data[0]#tr.data)
        if phase == 'SH': tr.detrend(type='linear')
        tr.data = 10 ** 6 * tr.data
        if np.max(np.abs(tr.data)) > 5 * 10 ** 5: continue
        stream[0] = tr
        stream.write('{}.tmp0'.format(sac), format='SAC', byteorder = 0)
    return
    
    
##########################
# We proceed to select teleseismic stations
##########################


def select_tele_stations(files, phase, tensor_info):
    """We select body and/or surface waves to use in finite fault modelling.
    
    :param files: list of used waveforms in sac format
    :param phase: string indicating whether wave is P or SH or Love or Rayleigh
    :type tensor_info: list
    :type phase: string
    """
    event_lat = tensor_info['lat']
    event_lon = tensor_info['lon']
    weight = 1.0 if not phase == 'SH' else 0.5
    if phase in ['P', 'Rayleigh']:
        min_snr = 6.0 if phase == 'P' else 4.0
        window = 100 if phase == 'P' else 1500
        jump = 2
        total = __used_stations(jump, files, tensor_info)
        if total > 50:
            jump = 4
            total = __used_stations(jump, files, tensor_info)
        if total > 50:
            jump = 8
    if phase in ['SH', 'Love']:
        min_snr = 3.0#4.0
        window = 200 if phase == 'SH' else 1500
        jump = 8
        total = __used_stations(jump, files, tensor_info)
        if total > 30:
            jump = 10
            total = __used_stations(jump, files, tensor_info)
        if total > 30:
            jump = 12

    sacheaders = [SACTrace.read(sac) for sac in files]
    signal2noise = [__s2nr(sacfile, phase, window) for sacfile in files]
    min_signal2noise = len(
        [ratio for ratio in signal2noise if ratio > min_snr])
    if min_signal2noise < 4: min_snr = min_snr / 2.0
    score = lambda x, y, z, y0, y1, z0:\
        x - (y - y0)*(y - y1)/((y - y0)**2.0 + (y - y1)**2.0)\
        + (z - z0)/(20 - z0)
    
# Limiting the number of data. Choose one station for every degree in
# azimuth

    phase = 'LONG' if phase in ['Rayleigh', 'Love'] else phase        
    with open(os.path.join(phase, 'selected_stations_list'), 'w') as outfile:
        for az0 in range(0, 360, jump):
            az1 = az0 + jump
            best = ''
            best_sta = ''
            best_score = -1.0
            for sacheader, snr, sac in zip(sacheaders, signal2noise, files):
                dis, az, baz = mng._distazbaz(
                    sacheader.stla, sacheader.stlo, event_lat, event_lon)
                # az = sacheader.az
                if az0 > az or az >= az1 or snr <= min_snr: continue
                value = 1 if kilometers2degrees(dis) >= 45 else 0
                value = value if phase in ['P', 'SH'] else 1.0
                value = score(value, az, snr, az0, az1, min_snr)
                if value > best_score:
                    best = '{}\nsnr: {}\nbody wave weight: {}\n'\
                    'surf wave weight: 1.0\n'.format(sac, snr, weight)
                    best_score = value
                    best_sta = sac

            if best_sta:
                outfile.write(best)
                copy2(best_sta, best_sta[:-1])

    return


def __s2nr(sacfile, phase, signal_length):
    r"""Signal to noise ratio for data.
    """
    sacheader = SACTrace.read(sacfile)
    stream = read(sacfile)
    trace = stream[0].data
    begin = sacheader.b
    dt = sacheader.delta
    arrival = sacheader.t1 if phase == 'P' else sacheader.t5
    error_window = 5.0 if phase == 'P' else 10.0
    if phase in ['Rayleigh', 'Love']:
        arrival = sacheader.t1
        error_window = 100.0
    arrival = int((arrival - begin) / dt)
    length = int(signal_length / dt + 0.1)
    error_window = int(error_window / dt)
    length = min(length, arrival - error_window)
    begin_signal = arrival + error_window - 1
    signal = trace[begin_signal:begin_signal + length]
    begin_noise = arrival - error_window - length + 1
    noise = trace[begin_noise:begin_noise + length]
    signal_std = np.std(signal)
    noise_std = np.std(noise)
    return min(20, signal_std / noise_std)


def __used_stations(jump, sacfiles, tensor_info):
    sacheaders = [SACTrace.read(sac) for sac in sacfiles]
    event_lat = tensor_info['lat']
    event_lon = tensor_info['lon']
    values = [mng._distazbaz(header.stla, header.stlo, event_lat, event_lon)\
              for header in sacheaders]
    azimuths = [az for dis, az, baz in values]
    total = 0
    for az0 in range(0, 360, jump):
        list_ = [az for az in azimuths if az0 <= az < az0 + jump]
        if list_:
            total = total + 1
    return total


##################################
# teleseismic surface waves
##################################


def select_process_surf_tele(tele_files0, tensor_info):
    """Module for selecting and processing surface wave data.
    
    :param tensor_info: dictionary with moment tensor information
    :param tele_files0: files with surface wave to be selected and processed
    :type tensor_info: dict
    :type tele_files0: list
    
    .. rubric:: Example:
    
    >>> from obspy.core.utcdatetime import UTCDateTime
    >>> import glob
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
            'centroid_depth': 35,
            'timedelta': 30 * 60
        }
    >>> tele_files0 = glob.glob('*sac')
    >>> select_process_surf_tele(tele_files0, tensor_info)
    .. note::
        
        Currently, this code allows only to process files in sac format, where the 
        instrumental response is in a SACPZ file. Other data formats are not
        supported.
    """
    timedelta = tensor_info['timedelta']
    depth = tensor_info['depth']
    tele_files1 = __pre_select_tele(tele_files0, tensor_info, timedelta)
    if not tele_files1:
        return
    if not os.path.isdir('logs'):
        os.mkdir('logs')
    if not os.path.isdir('LONG'):
        os.mkdir('LONG')
    logger1 = ml.create_log('surf_tele_processing',
                            os.path.join('logs', 'surf_tele_processing.log'))
    logger1 = ml.add_console_handler(logger1)
    logger1.info('Start processing of long period surface waves')

    data_folder = os.path.abspath(os.getcwd())
    
    os.chdir(data_folder)
    logger1.info('Get theoretical arrivals of P and S waves')
    model = TauPyModel(model="ak135f_no_mud")
    __picker(tensor_info, tele_files1, depth, model)
    logger1.info('Remove instrumental response for surface waves')
    response_files = glob.glob('SACPZ*') + glob.glob('SAC_PZs*')
    __remove_response_surf(tele_files1, response_files, logger=logger1)
    logger1.info('Rotate horizontal components for surface waves')
    os.chdir(os.path.join(data_folder, 'LONG'))
    horizontal_sacs = glob.glob('*BH[EN12]*DIS')
    __rotation(horizontal_sacs)
    logger1.info('Get final surface waves')
    os.chdir(os.path.join(data_folder, 'LONG'))
    surf_files = glob.glob('*.BHZ.*DIS') + glob.glob('*.SH.*DIS')
    __create_long_period(surf_files)
    
    logger1.info('Long period surface waves have been succesfully processed')
    os.chdir(data_folder)
    ml.close_log(logger1)
    files = glob.glob(os.path.join('LONG', '*BHZ*DIS.tmp0'))
    select_tele_stations(files, 'Rayleigh', tensor_info)
    files = glob.glob(os.path.join('LONG', '*SH*DIS.tmp0'))
    select_tele_stations(files, 'Love', tensor_info)
    return
    
    
def __remove_response_surf(used_tele_sacs, response_files, logger=None):
    """We remove instrumental response for surface wave data.
    """
    sacheaders = [SACTrace.read(sac) for sac in used_tele_sacs]
    new_name = lambda station, channel, network:\
        os.path.join('{}.{}.{}.DIS'.format(station, channel, network))

    sac_files = []
    old_responses = []
    names = []

    for sacheader, sac in zip(sacheaders, used_tele_sacs):
        name = sacheader.kstnm
        channel = sacheader.kcmpnm
        loc_code = sacheader.khole# if sacheader.khole else '--'
        # pz_files = [response for response in response_files\
        #             if name in response and channel in response\
        #             and loc_code in response]
        pz_files0 = [response for response in response_files\
                     if name in response and channel in response]
        if loc_code == '00':
            pz_files = [response for response in pz_files0 if '00' in response]
        if loc_code in [None, '']:
            pz_files = [response for response in pz_files0 if '--' in response]
            pz_files = pz_files +\
                [response for response in pz_files0 if '__' in response]
        # print(name, channel, loc_code)
        # print(pz_files)
        if not pz_files:
            continue
        pzfile_val = next(iter(pz_files))#
        netwk = sacheader.knetwk
        str_name = os.path.join('LONG', new_name(name, channel, netwk))
        if not os.path.isfile(pzfile_val):
            continue

        sac_files = sac_files + [sac]
        old_responses = old_responses + [pzfile_val]
        names = names + [str_name]

    pre_processing = 'dif \ntaper \nint \n'
    out, err = response.replace_response_sac(
        sac_files, old_responses, names, pre_processing=pre_processing,
        freq0=0.003, freq1=0.004, freq2=0.006, freq3=0.007)

    if logger:
        logger.debug(out.decode('utf-8'))
        if err: logger.warning(err.decode('utf-8'))
    return
    
    
def __create_long_period(surf_files):
    """We downsample long period data.
    """
    sacheaders = [SACTrace.read(sac) for sac in surf_files]
    for sac, sacheader in zip(surf_files, sacheaders):
        sacheader.kcmpnm = 'SH' if sac.split('.')[1] == 'SH'\
            else sacheader.kcmpnm
        sacheader.write(sac, byteorder = 'little')
        
    for sac in surf_files:
        st = read(sac)
        tr = st[0]
        factor = 4 * tr.stats.sampling_rate
        tr.decimate(int(round(factor)), no_filter=True)
        tr.data = 1000 * tr.data
        if np.max(np.abs(tr.data)) < 10 ** -3: continue
        st[0] = tr
        st.write('{}.tmp0'.format(sac), format='SAC', byteorder = 0)
    return
    

###################################
# cgps
###################################


def select_process_cgps(cgps_files, tensor_info, data_prop):
    """Routine for selecting and processing cgps data
    
    :param tensor_info: dictionary with moment tensor information
    :param data_prop: dictionary with waveform properties
    :param cgps_files: files with cGPS data to be selected and processed
    :type tensor_info: dict
    :type data_prop: dict
    :type cgps_files: list
    .. rubric:: Example:
    
    >>> from obspy.core.utcdatetime import UTCDateTime
    >>> import glob
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
            'centroid_depth': 35,
            'timedelta': 30 * 60
        }
    >>> cgps_files0 = glob.glob('*sac')
    >>> select_process_cgps(cgps_files0, tensor_info, data_prop)
    .. note::
        
        Currently, this code allows only to process files in sac format.
        Other data formats are not supported.
    """
    if not os.path.isdir('logs'):
        os.mkdir('logs')
    logger1 = ml.create_log('cgps_processing',
                            os.path.join('logs', 'cgps_processing.log'))
    logger1 = ml.add_console_handler(logger1)
    logger1.info('Process cGPS data')
    data_folder = os.path.abspath(os.getcwd())
    separator = '\n*************************************************\n'
    logger1.info('Select cGPS traces')
    __select_cgps_files(cgps_files, tensor_info)
    logger1.info('Process selected cGPS traces')
    os.chdir(os.path.join(os.getcwd(), 'cGPS'))
    cgps_files1 = glob.glob('*.L[HXY]*.sac') + glob.glob('*.L[HXY]*SAC')
    __change_start(cgps_files1, tensor_info)
    new_process_cgps(tensor_info, cgps_files1, data_prop, logger=logger1)
    os.chdir(data_folder)
    logger1.info('cGPS waves have been succesfully processed')
    ml.close_log(logger1)
    return


def __select_distance(tensor_info, station_lat, station_lon,
                      max_distance=180, min_distance=0, use_centroid=False):
    """
    """
    event_lat = tensor_info['lat']
    event_lon = tensor_info['lon']
    distance = locations2degrees(event_lat, event_lon, station_lat, station_lon)
    if use_centroid:
        centroid_lat = tensor_info['centroid_lat']
        centroid_lon = tensor_info['centroid_lon']
        distance2 = locations2degrees(
            centroid_lat, centroid_lon, station_lat, station_lon)
        distance = min(distance, distance2)
    return min_distance <= distance < max_distance


def __select_time(sac, tensor_info):
    """
    """
    datetime = tensor_info['datetime']
    origin_time = UTCDateTime(datetime)
    st = read(sac)
    return st[0].stats.starttime - origin_time < 20 * 60

    
def __select_cgps_files(cgps_files, tensor_info):
    """We select cgps data. We must make sure our cgps data is not too bad
    """
    event_lat = tensor_info['lat']
    event_lon = tensor_info['lon']
    depth = tensor_info['depth']
    time_shift = tensor_info['time_shift']
    centroid_lat = tensor_info['centroid_lat']
    streams = [read(file) for file in cgps_files]
    event_lon = tensor_info['lon']
    date_origin = tensor_info['date_origin']
    depth = tensor_info['depth']
    distance = 2 if time_shift < 50 else 4

    for i, (stream, sac) in enumerate(zip(streams, cgps_files)):
        station_lat = stream[0].stats.sac.stla
        station_lon = stream[0].stats.sac.stlo
        if station_lat == None or station_lon == None:
            continue
        syn_len = 2 * time_shift + 7 * depth / 50
        start = date_origin
        starttime = stream[0].stats.starttime
        delta = min(start - starttime, 20)
        stream.trim(starttime=start - delta)
        start, end = [stream[0].stats.starttime, stream[0].stats.endtime]
        if end - start < syn_len:
            continue
        if not __select_distance(tensor_info, station_lat, station_lon,
                                 max_distance=distance, use_centroid=True):
            continue
        #print(stream[0].stats.station)
        indexes = np.isfinite(stream[0].data)        
        if np.max(stream[0].data[indexes]) == np.min(stream[0].data[indexes]):
            continue
        
        stream[0].data = stream[0].data\
            - np.sum(stream[0].data[np.isfinite(stream[0].data)])\
            / np.sum(np.isfinite(stream[0].data))
        
        if np.max(np.abs(stream[0].data)) > 10 ** 4:
            continue

        if not np.sum(np.isnan(stream[0].data)) == 0:
            station_lat = stream[0].stats.sac.stla
            station_lon = stream[0].stats.sac.stlo
            delta = stream[0].stats.delta
            distance = locations2degrees(
                event_lat, event_lon, station_lat, station_lon)
            distance = 111.11 * distance
            distance2 = np.sqrt(distance ** 2 + depth ** 2)
    
            index0 = distance2 / 5 / delta
            indexes2 = [i for i, value in enumerate(stream[0].data)\
                if not np.isfinite(value)]
            indexes3 = [i for i in indexes2 if index0 <= i < index0 + 100]
            gaps = __get_gaps(stream[0].data)
            if not max(indexes2) < index0\
            and len(indexes3) >= 5 and max([len(gap) for gap in gaps]) >= 4:
                continue
            else:
                for gap in gaps[:-1]:
                    begin = gap[0] - 1
                    end = gap[-1] + 1
                    stream[0].data = __linear_fill(stream[0].data, begin, end)

        station = stream[0].stats.station
        channel = stream[0].stats.channel
        name = os.path.join('cGPS', '{}.{}.sac'.format(station, channel))
        stream.write(name, format='SAC', byteorder=0)


def __change_start(stations_str, tensor_info):
    """Routine for modifying start time of cGPS data
    """
    sacheaders = [SACTrace.read(sac) for sac in stations_str]
    for sac, sacheader in zip(stations_str, sacheaders):
        st = read(sac)
        begin = sacheader.b
        reftime = sacheader.reftime
        st[0].stats.starttime = UTCDateTime(reftime) + begin
        st[0].trim(starttime=tensor_info['date_origin'] - 110,
            endtime=tensor_info['date_origin'] + 1100)# - 60 * 10)
        st[0].data = st[0].data - np.mean(st[0].data[:10])
        st.write(sac, format='SAC', byteorder=0)


def new_process_cgps(tensor_info, stations_str, data_prop, logger=None):
    """Routine for processing (filtering and filling header) selected cGPS data.
    
    :param tensor_info: dictionary with moment tensor information
    :param stations_str: list with waveforms (in sac format) to use
    :param data_prop: dictionary with waveform properties
    :param logger: where to log results of processing data
    :type tensor_info: dict
    :type stations_str: list
    :type data_prop: dict
    :type logger: logging, optional
    """
    start = tensor_info['date_origin']
    filtro_strong = data_prop['strong_filter']
    lat_ep = tensor_info['lat']
    lon_ep = tensor_info['lon']
    sampling = data_prop['sampling']
    dt_cgps = sampling['dt_cgps']       
       
    high_freq = min(filtro_strong['high_freq'], 0.5)

    streams = [read(sac) for sac in stations_str]
    
    for st, sac in zip(streams, stations_str):
        if dt_cgps > st[0].stats.delta:
            ratio = int(dt_cgps / st[0].stats.delta)
            st[0].decimate(ratio)
        if st[0].stats.delta > dt_cgps:
            st[0].resample(sampling_rate=1/dt_cgps)
        sacheader = SACTrace.from_obspy_trace(st[0])
        sacheader.t9 = high_freq
        sacheader.write(sac, byteorder='little')

    sacheaders = [SACTrace.read(sac) for sac in stations_str]        
    
    for sac, sacheader in zip(stations_str, sacheaders):
        begin = sacheader.b
        start_time = UTCDateTime(sacheader.reftime) + begin
        sacheader.o = start - start_time
        sacheader.write(sac, byteorder='little')
    
    _filter_decimate(stations_str, filtro_strong, dt_cgps, corners=4, passes=2,
                     decimate=False, filter='lowpass', logger=logger)
    return


def __linear_fill(data, begin, end):
    """
    """
    slope = (data[begin] - data[end]) / (begin - end) if begin >= 0 else 0
    intercept = data[begin] - slope * begin if begin >= 0 else data[end]
    for index in range(begin + 1, end):
        data[index] = slope * index + intercept
    return data


def __get_gaps(data):
    """
    """
    indexes = [i for i, value in enumerate(data) if not np.isfinite(value)]
    if not indexes:
        return None
    start_gaps = [index for index in indexes if index - 1 not in indexes]
    end_gaps = [index for index in indexes if index + 1 not in indexes]
    gaps = []
    for start, end in zip(start_gaps, end_gaps):
        gaps = gaps + [[index for index in range(start, end + 1)]]
    return gaps


###################################
# custom strong motion
###################################


def _filter_decimate(sac_files, filtro, dt, corners=4, passes=2,
                     decimate=True, logger=None, filter='bandpass'):
    """
    """
    low_freq = filtro['low_freq']
    high_freq = filtro['high_freq']
    power = lambda m, n: max([d for d in range(10) if m % (n**d) == 0])
    input_sac = ''
    for sac in sac_files:
        sacheader = SACTrace.read(sac)
        delta = sacheader.delta
        ratio = int(round(dt / delta))
        power5 = power(ratio, 5)
        power2 = power(ratio, 2)
        power3 = power(ratio, 3)
        if ratio > 2**power2 * 3**power3 * 5**power5: continue
        input_sac = '{}\nread {} \n'.format(input_sac, sac)
        if filter == 'bandpass':
            input_sac = input_sac + 'bp c {} {} n {} p {} \n'.format(
                        low_freq, high_freq, corners, passes)
        elif filter == 'lowpass':
            input_sac = input_sac + 'lp c {} n {} p {} \n'.format(
                        high_freq, corners, passes)
        if decimate:
            input_sac = input_sac\
                    + 'decimate 5\n' * power5\
                    + 'decimate 3\n' * power3\
                    + 'decimate 2\n' * power2
        input_sac = input_sac + 'mul 100\nwrite {}.tmp \n'.format(sac)

    input_sac = '{}\nquit\n'.format(input_sac)
    out, err = mng.run_sac(input_sac)
    if logger:
        logger.debug(out.decode('utf-8'))
        if err: logger.warning(err.decode('utf-8'))
    return


#####################
# automatic strong motion
#####################



def select_process_strong(strong_files0, tensor_info, data_prop, custom=False):
    """Module for selecting and processing strong motion data
    
    :param tensor_info: dictionary with moment tensor information
    :param data_prop: dictionary with waveform properties
    :param strong_files0: files with strong motions to be selected and processed
    :type tensor_info: dict
    :type data_prop: dict
    :type strong_files0: list
    .. rubric:: Example:
    
    >>> from obspy.core.utcdatetime import UTCDateTime
    >>> import glob
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
            'centroid_depth': 35,
            'timedelta': 30 * 60
        }
    >>> strong_files0 = glob.glob('*sac')
    >>> select_process_strong(strong_files0, tensor_info, data_prop)
    .. note::
        Currently, this code allows only to process files in sac format, where the 
        instrumental response is in a SACPZ file. Other data formats are not
        supported.
    """
    data_folder = os.path.abspath(os.getcwd())
    if not os.path.isdir('logs'):
        os.mkdir('logs')
    if not os.path.isdir('STR'):
        os.mkdir('STR')
    logger1 = ml.create_log('strong_motion_processing',
                            os.path.join('logs', 'strong_motion_processing.log'))
    logger1 = ml.add_console_handler(logger1)
    logger1.info('Process strong motion data')
    logger1.info('Remove response for selected strong motion traces')
    response_files = glob.glob('SACPZ*') + glob.glob('SAC_PZs*')
    for file in response_files:
        __convert_response_acc(file)
    __remove_response_str(strong_files0, response_files, logger=logger1)
    logger1.info('Select strong motion traces')
    os.chdir(os.path.join(data_folder, 'STR'))
    strong_files1 = glob.glob('*ACC')
    strong_files2 = __select_str_files(strong_files1, tensor_info)
    logger1.info('Update duration of traces')
    _trim_accel(strong_files2, tensor_info)
    logger1.info('Process selected strong motion traces')
    process_strong_motion(strong_files2, tensor_info, data_prop, logger=logger1)
    os.chdir(data_folder)
    strong_files3 = glob.glob(os.path.join('STR', '*tmp'))
    select_strong_stations(strong_files3)
    logger1.info('Strong motion waves have been succesfully processed')
    ml.close_log(logger1)
    return


def __convert_response_acc(resp_file):
    """
    """
    with open(resp_file, 'r') as infile:
        lines = [line for line in infile]
    indexes = [i for i, line in enumerate(lines) if 'ZEROS' in line.split()]
    index0 = 0
    lines2 = []
    for index in indexes:
        lines2 = lines2 + lines[index0:index]
        string, nzeroes = lines[index].split()
        if nzeroes == '2':
            lines2 = lines2 + ['ZEROS\t 0\n']
            index0 = index + 3
        if nzeroes == '3':
            lines2 = lines2 + ['ZEROS\t 0\n']
            index0 = index + 4
        lines2 = lines2 + lines[index0:]
    with open(resp_file, 'w') as outfile:
        for line in lines2:
            outfile.write(line)
    return resp_file
    
    
def __select_str_files(strong_files, tensor_info):
    r"""We pre-select strong motion data
    """
    time_shift = tensor_info['time_shift']
    centroid_lat = tensor_info['centroid_lat']
    centroid_lon = tensor_info['centroid_lon']
    depth = tensor_info['depth']
    strong_files2 = []
    distance = 2 if time_shift < 50 else 4
    
    for sac in strong_files:
        select_stat = True
        sacfile = SACTrace.read(sac)
        station_lat = sacfile.stla
        station_lon = sacfile.stlo
      
        st = read(sac)
        syn_len = 20 + 4 * time_shift + 7 * depth / 50
        start, end = [st[0].stats.starttime, st[0].stats.endtime]
        if end - start < syn_len:
            continue
        if not __select_distance(tensor_info, station_lat, station_lon,
                                 max_distance=distance, use_centroid=True):
            continue
        strong_files2 = strong_files2 + [sac]
    
    return strong_files2


def __remove_response_str(stations_str, response_files, logger=None):
    """Remove instrumental response.
    """
    sacheaders = [SACTrace.read(sac) for sac in stations_str]

    channel_name = lambda x: '{}N{}'.format(x[0], x[2])
    new_name = lambda station, channel, network:\
        'STR.{}.{}.{}.ACC'.format(station, channel_name(channel), network)
    # channel_name = lambda x: '{}H{}'.format(x[0], x[2])
    # new_name = lambda station, channel, network:\
    #     'STR.{}.{}.{}.VEL'.format(station, channel_name(channel), network)

    sac_files = []
    old_responses = []
    names = []
    pre_processing = 'rtr\n'

    for sacheader, sac in zip(sacheaders, stations_str):
        name = sacheader.kstnm
        channel = sacheader.kcmpnm
        network = sacheader.knetwk
        pz_files = [response for response in response_files\
                    if name in response and channel in response]
        if not pz_files:
            continue
        pzfile2 = next(iter(pz_files))
        if not os.path.isfile(pzfile2):
            continue
        filename = os.path.join('STR', new_name(name, channel, network))

        sac_files = sac_files + [sac]
        old_responses = old_responses + [pzfile2]
        names = names + [filename]

    out, err = response.replace_response_sac(
        sac_files, old_responses, names, pre_processing=pre_processing)

    if logger:
        logger.debug(out.decode('utf-8'))
        if err: logger.warning(err.decode('utf-8'))
    return


def _trim_accel(strong_files, tensor_info):
    """
    """
    start = tensor_info['date_origin']
    streams = [read(sac) for sac in strong_files]
    for sac, stream in zip(strong_files, streams):
        starttime = stream[0].stats.starttime
        diff = start - starttime
        # if diff > 20:
        stream[0].trim(starttime = start - 20, pad=True,
            fill_value=0)#stream[0].data[0])
        header = SACTrace.from_obspy_trace(stream[0])
        stream.write(sac, format='SAC', byteorder=0)


def process_strong_motion(strong_files, tensor_info, data_prop, logger=None):
    r"""We integrate strong motion to velocity and proceed to downsample it.
    
    :param tensor_info: dictionary with moment tensor information
    :param data_prop: dictionary with waveform properties
    :param strong_files: files with strong motion to be resampled and integrated to velocity
    :param logger: where to log results of strong motion resampling and integration
    :type tensor_info: dict
    :type data_prop: dict
    :type strong_files: list
    :type logger: logging, optional
    """
    muestreo = data_prop['sampling']
    dt_strong = muestreo['dt_strong']
    filtro_strong = data_prop['strong_filter']
    
    low_freq = filtro_strong['low_freq']       
    high_freq = filtro_strong['high_freq']
    
    sacheaders = [SACTrace.read(sac) for sac in strong_files]
    
    for sacheader, sac in zip(sacheaders, strong_files):
        sacheader.t8 = low_freq
        sacheader.t9 = high_freq
        sacheader.write(sac, byteorder='little')
    
    print('Baseline removal procedure')
    cpus = 6#cpu_count()
    lists = [strong_files[i:i + cpus] for i in range(0, len(strong_files), cpus)]
    pool = Pool(cpus)
    results = [
        pool.apply_async(__worker, args=(sta_list, ))\
        for sta_list in lists]
    [result.get() for result in results]
    pool.close()
    pool.join()
    strong_files2 = glob.glob('*VEL')
    
    _filter_decimate(strong_files2, filtro_strong, dt_strong)
    return


def small_events(select_str_data):
    """A simple method for removing strong motion baselines. Based on Boore et
    al. [2002]_
    
    :param select_st_data: list with waveforms (in sac format) to use
    :type select_st_data: list
    
    This method consists on padding the record with zeros on the extremes and
    applying a low pass filter.
    
    .. [2002] Boore DM: On pads and filters: processing strong-motion data.
        Bulletin of the Seismological Society of America 2005,95(2):745-750.
        10.1785/0120040160
    """
    order = 4
    for sac in select_str_data:
        stream = read(sac)
        delta = stream[0].stats.delta
        t_pad = 1.5 * order / 0.02 / 2
        t_pad = int(t_pad / delta)
        fs = stream[0].stats['sampling_rate']
        if _delete_criteria(stream[0].data):
            continue
        data = stream[0].data
        data = np.pad(data, t_pad, 'constant')
        data = cumtrapz(data, dx=delta, initial=0)
        st_vel = stream.copy()
        nyq = 0.5 * fs
        low = 0.02 / nyq
        b, a = butter(4, low, btype='highpass')
#
# apply filter
#
        data = filtfilt(b, a, data)
        st_vel[0].data = data[t_pad:-t_pad]
        st_vel.write('{}.VEL'.format(sac[:-4]), format='SAC')
    return


def __worker(select_str_data):
    """Helper routine. Method for baseline removal based on Wang et al.
    """
    for sac in select_str_data:
        stream = read(sac)
        if _delete_criteria(stream[0].data):
            continue
        st_vel = wang1.wang_process(sac)
        if not st_vel:
            continue
        st_vel.write('{}.VEL'.format(sac[:-4]), format='SAC')
    return


def _delete_criteria(data):
    """
    """
    data = data - np.mean(data)
    if np.max(np.abs(data)) < 10 ** -6:
        return True
    return False


##########################
# We proceed to select strong motion stations
##########################


def select_strong_stations(files):
    """We select strong motion data to use in finite fault modelling
    
    :param files: list of waveform files (in sac format) to select
    :type files: list
    """
    if len(files) < 150:
        return files
    streams = [read(sac) for sac in files]
    sacheaders = [SACTrace.read(sac) for sac in files]
    names = [header.kstnm for header in sacheaders]
    azimuth0 = np.array([header.az for header in sacheaders])
    az0 = np.amin(azimuth0)
    az1 = np.amax(azimuth0)
    jump = int((az1 - az0) / 60) + 1
    names = list(set(names))
    select_stations = []
    
# Limiting the number of data. Choose at most one station (3 channels)
# for every degree in azimuth. Preferably use data with higher PGV
     
    for az0 in range(0, 360, jump):
        az1 = az0 + jump
        best_pgv = 0.0
        for name in names:
            indexes = [i for i, header in enumerate(sacheaders)\
                       if name==header.kstnm]
            azimuth = next(header.az for header in sacheaders\
                           if name==header.kstnm)
            if not az0 <= azimuth < az1:
                continue
            streams2 = [st for i, st in enumerate(streams) if i in indexes]
            pgvs = np.array([np.max(np.abs(st[0].data)) for st in streams2])
            pgv = np.sqrt(np.sum(pgvs ** 2))
            if pgv > best_pgv:
                best_pgv = pgv
                chosen_sta = [sac for i, sac in enumerate(files)\
                              if i in indexes]
        
        if best_pgv > 0.0:
            select_stations = select_stations + chosen_sta
    
    for file in files:
        if not file in select_stations:
            os.remove(file)
    return select_stations
    
    
if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--folder", default=os.getcwd(),
                        help="folder where there are input files")
    parser.add_argument("-gcmt", "--gcmt_tensor",
                        help="location of GCMT moment tensor file")
    parser.add_argument("-t", "--tele", action="store_true",
                        help="process teleseismic data")
    parser.add_argument("-su", "--surface", action="store_true",
                        help="process surface waves data")
    parser.add_argument("-st", "--strong", action="store_true",
                        help="process strong motion data")
    parser.add_argument("--cgps", action="store_true",
                        help="process cGPS data")
    parser.add_argument("-dt", "--sampling_delta",
                        help="process cGPS data")
    parser.add_argument("--pick", action="store_true",
                        help="pick of teleseismic waves")
    args = parser.parse_args()
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
    if args.tele:      
        tele_files = glob.glob('*BH*SAC') + glob.glob('*BH*sac')
        select_process_tele_body(tele_files, tensor_info, data_prop)
    if args.surface:
        tele_files = glob.glob('*BH*SAC') + glob.glob('*BH*sac')
        select_process_surf_tele(tele_files, tensor_info)
    if args.strong:
        strong_files = glob.glob('*.HN*SAC') + glob.glob('*.HL*SAC')\
                    + glob.glob('*.HN*sac') + glob.glob('*.HL*sac')\
                    + glob.glob('*.AH?.*') + glob.glob('*.AH?.*')\
                    + glob.glob('*_HN*sac') + glob.glob('*_HL*sac')
        # strong_motion_processing_custom(
        #     strong_files, tensor_info, data_prop, type='acc')
        select_process_strong(strong_files, tensor_info, data_prop)
    if args.cgps:
        cgps_files = glob.glob('*.L[HXY]*SAC') + glob.glob('*.L[HXY]*sac')
        select_process_cgps(cgps_files, tensor_info, data_prop)
    
