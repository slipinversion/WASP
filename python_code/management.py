# -*- coding: utf-8 -*-
"""Routine for performing administrative tasks, such as changing folders, or
moving to a different folder.
"""


import os
import numpy as np
import json
from obspy.geodetics import locations2degrees, degrees2kilometers
from obspy.core.utcdatetime import UTCDateTime
import subprocess
from read_config import read_config


def theoretic_arrivals(model, dist, depth):
    """
    """
    p_arrival = model.get_travel_times(
            source_depth_in_km=depth, distance_in_degree=dist,
            phase_list=['P'])
    pp_arrival = model.get_travel_times(
            source_depth_in_km=depth, distance_in_degree=dist,
            phase_list=['PP'])
    s_arrival = model.get_travel_times(
            source_depth_in_km=depth, distance_in_degree=dist,
            phase_list=['S'])
    ss_arrival = model.get_travel_times(
            source_depth_in_km=depth, distance_in_degree=dist,
            phase_list=['SS'])
    p_arrivals = model.get_travel_times(
            source_depth_in_km=depth, distance_in_degree=dist + 0.1,
            phase_list=['P'])
    p_slowness = (p_arrivals[0].time - p_arrival[0].time) / 0.1
    s_arrivals = model.get_travel_times(
            source_depth_in_km=depth, distance_in_degree=dist + 0.1,
            phase_list=['S'])
    s_slowness = (s_arrivals[0].time - s_arrival[0].time) / 0.1
    arrivals = {
        'p_arrival': p_arrival,
        'pp_arrival': pp_arrival,
        's_arrival': s_arrival,
        'ss_arrival': ss_arrival,
        'p_slowness': p_slowness,
        's_slowness': s_slowness
    }
    return arrivals
    
    
def update_data(tensor_info, data_type=None):
    """
    
    :param data_type: list of data types to be used in modelling.
    :param tensor_info: dictionary with moment tensor information
    :type data_type: list, optional
    :type tensor_info: dict
    """
    if not data_type:
        data_type = [
                'tele_body',
                'surf_tele',
                'strong_motion',
                'cgps',
                'gps'
        ]
#    data_type = inversion_data(tensor_info) if not data_type else data_type
    data_type2 = []
    if 'tele_body' in data_type:
        data_type2 = data_type2 + ['tele_body']\
            if os.path.isfile('tele_waves.json') else data_type2
    if 'surf_tele' in data_type:
        data_type2 = data_type2 + ['surf_tele']\
            if os.path.isfile('surf_waves.json') else data_type2
    if 'strong_motion' in data_type:
        data_type2 = data_type2 + ['strong_motion']\
            if os.path.isfile('strong_motion_waves.json') else data_type2
    if 'cgps' in data_type:
        data_type2 = data_type2 + ['cgps']\
            if os.path.isfile('cgps_waves.json') else data_type2
    if 'gps' in data_type:
        data_type2 = data_type2 + ['gps']\
            if os.path.isfile('static_data.json') else data_type2
#    os.chdir(current_dir)
    return set(data_type2) & set(data_type) 
    
    
def default_dirs():
    """Environment variables.
    """
    config = read_config()
    paths = config['PATHS']
    info = paths['info']
    compute_near_gf = paths['compute_near_gf']
    get_near_gf = paths['get_near_gf']
    modelling = paths['modelling']
    default_dirs = {
            'root_dir': paths['code_path'],
            'long_gf_bank': paths['surf_gf_bank'],
            'crust_codes': os.path.join(info, 'CNtype2.txt'),
            'models_codes': os.path.join(info, 'CNtype2_key.txt'),
            'litho_model': os.path.join(info, 'LITHO1.0.nc'),
            'gf_bank': paths['surf_gf_bank'],
            'strong_motion_gf_bank2': os.path.join(
                    compute_near_gf, 'green_bank_openmp_f95'),
            'strong_motion_gf': os.path.join(get_near_gf, 'get_strong_motion'),
            'cgps_gf_bank': os.path.join(get_near_gf, 'cgps'),
            'gps_gf': os.path.join(compute_near_gf, 'gf_static_f95'),
            'tele_gf': os.path.join(modelling, 'green_tele'),
            'finite_fault': os.path.join(modelling, 'run_modelling'),
            'forward': os.path.join(modelling, 'run_forward'),
            'trench_graphics': os.path.join(
                    paths['cartopy_files'], 'PB2002_plates'),
            'sac_exec': paths['sac_exec']
    }
    return default_dirs


def start_time_id(tensor_info):
    """We get a string, by concatenating the data for the origin time of
    an earthquake
    
    :param tensor_info: dictionary with moment tensor information
    :type tensor_info: dict
    """
    datetime = tensor_info['datetime']
    chunk1, chunk2 = datetime.split('T')
    year, month, day = chunk1.split('-')[:3]
    hour, minute, second = chunk2.split(':')

    time_origin = '{}{}{}{}{}{}'.format(year, month, day, hour, minute, second)
    return time_origin


def _distazbaz(station_lat, station_lon, event_lat, event_lon):
    """We compute the distance, azimuth and back_azimuth, between two pairs
    lat lon.
    """
    degrees2rad = np.pi / 180.0
    arco_circulo = locations2degrees(
            event_lat, event_lon, station_lat, station_lon)
    distance = degrees2kilometers(arco_circulo)
    azimuth = np.arctan2(
            np.cos(station_lat * degrees2rad) * np.cos(event_lat * degrees2rad)\
            * np.sin((station_lon - event_lon) * degrees2rad),
            np.sin(station_lat * degrees2rad) - np.cos(arco_circulo * degrees2rad)\
            * np.sin(event_lat * degrees2rad))
    azimuth = azimuth / degrees2rad
    azimuth = np.remainder(azimuth, 360)

    sin_comp_baz = np.sin(
            (station_lon - event_lon) * degrees2rad)\
            * np.sin((90 - event_lat) * degrees2rad)\
            / np.sin(arco_circulo * degrees2rad)
    back_azimuth = 2 * np.pi - np.arcsin(sin_comp_baz)
    back_azimuth = back_azimuth / degrees2rad
    back_azimuth = np.remainder(back_azimuth, 360)
    return distance, azimuth, back_azimuth


def __bool_finite_fault(tensor_info):
    """
    """
    moment = tensor_info['moment_mag']
    lat = tensor_info['lat']
    lon = tensor_info['lon']
    finite_fault = False if moment < 8 * 10**25 else True
    if 30 < lat < 40 and -125 < lon < -110:
        finite_fault = False if moment < 10**25 else True
    return finite_fault


def correct_response_file(tensor_info, pzfile):
    """Routine for selecting instrumental response only in period of earthquake
    """
    date_origin = tensor_info['date_origin']
    with open(pzfile, 'r') as infile:
        lines = [line.split() for line in infile]
    start_times = [line[-1] for line in lines if 'START' in line]
    end_times = [line[-1] for line in lines if 'END' in line]
    end_times[0] = '2030-12-31T23:59:59'
    print(start_times)
    print(end_times)
    print('\n\n\n\n')
    network_lines = [index for index, line in enumerate(lines) if 'NETWORK' in line]
    constant_lines = [index for index, line in enumerate(lines) if 'CONSTANT' in line]
    
    start_times2 = [UTCDateTime(time) for time in start_times]
    end_times2 = [UTCDateTime(time) for time in end_times]
    
    zipped1 = zip(start_times2, end_times2)
    
    index = next(i for i, (start, end) in enumerate(zipped1)\
                 if start <= date_origin <= end)
    index = max(0, index - 1)
    index1 = network_lines[index] - 1
    index2 = constant_lines[index] + 1
    
    with open(pzfile, 'w') as outfile:
        for line in lines[index1:index2]:
            line2 = ' '.join(line)
            outfile.write('{}\n'.format(line2))


def run_sac(command):
    """Routine to execute sac with a given input command from python
    """
    dirs = default_dirs()
    sac_exec = dirs['sac_exec']
    p = subprocess.Popen(
        sac_exec, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
    out, err = p.communicate(input=command.encode('utf-8'))
    p.terminate()
    return out, err


if __name__ == '__main__':
    import glob
    import argparse
    import seismic_tensor as tensor

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--folder", default=os.getcwd(),
                        help="folder where there are input files")
    parser.add_argument("-gcmt", "--gcmt_tensor",
                        help="location of GCMT moment tensor file")
    
    args = parser.parse_args()
    os.chdir(args.folder)
    files = glob.glob('SACPZ*')
    tensor_info = tensor.get_tensor(cmt_file=args.gcmt_tensor)
    for pzfile in files:
        correct_response_file(tensor_info, pzfile)
        
    

