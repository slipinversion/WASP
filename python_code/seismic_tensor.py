# -*- coding: utf-8 -*-
"""
Module for retrieving data from moment tensor files in CMTSOLUTION format, and
retrieve the nodal planes from such a moment tensor.
"""


import numpy as np
from obspy.core.utcdatetime import UTCDateTime
import json
import os


def get_tensor(cmt_file=None, timedelta=None):
    """From a cmt tensor file we get tensor information in a dictionary.
    
    :param cmt_file: location of text file to be used
    :param timedelta: difference between O.T. and run time
    :type cmt_file: string, optional
    :type timedelta: float, optional
    """
    if cmt_file:
        tensor_info = read_gcmt_file(cmt_file)
        tensor_info = modify_tensor(tensor_info)
        tensor_info['timedelta']\
            = UTCDateTime.utcnow() - tensor_info['date_origin']
        if timedelta:
            tensor_info['timedelta'] = int(timedelta) * 60
    else:
        if not os.path.isfile('tensor_info.json'):
            raise RuntimeError('No file named tensor_info.json located in '\
                               'folder {}'.format(os.getcwd()))
        tensor_info = json.load(open('tensor_info.json'))
        tensor_info = modify_tensor(tensor_info)
    return tensor_info


def write_tensor(tensor_info):
    """We write a JSON with the information of the tensor
    
    :param tensor_info: dictionary with moment tensor information
    :type tensor_info: dict
    """
    tensor_info2 = tensor_info.copy()
    with open('tensor_info.json','w') as f:
        try:
            del tensor_info2['date_origin']
        except KeyError:
            pass
        json.dump(
            tensor_info2, f, sort_keys=True, indent=4, separators=(',', ': '),
            ensure_ascii=False)
    return


def read_gcmt_file(cmt_file):
    """We read the moment tensor information from a GCMT moment tensor file,
    and store its content in a python dictionary
    
    :param cmt_file: location of text file to be used
    :type cmt_file: string
    """
    with open(cmt_file, 'r') as input_file:
        lines = [line.split() for line in input_file]

    if 'Q' in lines[0][0]:
        new_line = lines[0]
        code, year = lines[0][0].split('Q')
        lines[0] = [code, year] + new_line[1:]
    if 'W' in lines[0][0]:
        new_line = lines[0]
        code, year = lines[0][0].split('W')
        lines[0] = [code, year] + new_line[1:]
        
    print(lines[0])
    
    year, month, day, hour, minute, second0, lat, lon, hyp_depth\
        = lines[0][1:10]
    second0 = float(second0)
    if second0 - int(second0) > 0:
        second, msec = str(second0).split('.')
    else:
        second = str(second0).split('.')[0]
        msec = '000'
    month = str(month).zfill(2)
    day = str(day).zfill(2)
    hour = str(hour).zfill(2)
    minute = str(minute).zfill(2)
    second = str(second).zfill(2)
    if len(msec) == 1:
        msec = str(100 * int(msec))
    else:
        msec = str(10 * int(msec))
    msec = msec.zfill(3)
    datetime = '{}-{}-{}T{}:{}:{}'.format(
        year, month, day, hour, minute, second)
    print(lines[2])
    time_shift = max(float(lines[2][2]), 5)
    print(lines[3])
    half_duration = float(lines[3][2])

    centroid_lat = float(lines[4][1])
    centroid_lon = float(lines[5][1])
    centroid_depth = float(lines[6][1])

    mrr = float(lines[7][1])
    mtt = float(lines[8][1])
    mpp = float(lines[9][1])
    mrt = float(lines[10][1])
    mrp = float(lines[11][1])
    mtp = float(lines[12][1])

    input_dict = {
        'mrr': mrr,
        'mtt': mtt,
        'mpp': mpp,
        'mrt': mrt,
        'mrp': mrp,
        'mtp': mtp,
        'datetime': datetime,
        'lat': lat,
        'lon': lon,
        'depth': hyp_depth,
        'time_shift': time_shift,
        'half_duration': half_duration,
        'centroid_lat': centroid_lat,
        'centroid_lon': centroid_lon,
        'centroid_depth': centroid_depth
    }

    return input_dict
    
    
def modify_tensor(tensor_info):
    """Add some extra information to the tensor dictionary.

    :param tensor_info: dictionary with moment tensor information
    :type tensor_info: dict
    """
    datetime = tensor_info['datetime']
    tensor_info = {
        key: float(value) for (key, value) in tensor_info.items()\
        if __is_number(value)}
    
    chunk1, chunk2 = datetime.split('T')
    year, month, day = chunk1.split('-')[:3]
    hour, minute, second = chunk2.split(':')

    tensor_info['date_origin'] = UTCDateTime(
        int(year), int(month), int(day), int(hour), int(minute),
        int(second))
    tensor_info['datetime'] = datetime

    mzz = tensor_info['mrr']
    mxx = tensor_info['mtt']
    myy = tensor_info['mpp']
    mxz = tensor_info['mrt']
    myz = - tensor_info['mrp']
    mxy = - tensor_info['mtp']
    moment_tensor = np.zeros((3, 3))
    moment_tensor[0, :] = np.array([mxx, mxy, mxz])
    moment_tensor[1, :] = np.array([mxy, myy, myz])
    moment_tensor[2, :] = np.array([mxz, myz, mzz])

    tensor_info['moment_mag'] = __get_moment_mag(moment_tensor)
    return tensor_info
    
    
def planes_from_tensor(tensor_info, dip_sens=0, stk_sens=0):
    """From moment tensor, we generate both nodal planes, in other words,
    two tuples, :math:`(\phi, \delta, \lambda)`.
    
    :param tensor_info: dictionary with moment tensor information
    :param dip_sens: modify dip angle by specified amount
    :param stk_sens: modify strike angle by specified amount
    :type tensor_info: dict
    :type dip_sens: float, optional
    :type stk_sens: float, optional
    """
    mzz = tensor_info['mrr']
    mxx = tensor_info['mtt']
    myy = tensor_info['mpp']
    mxz = tensor_info['mrt']
    myz = - tensor_info['mrp']
    mxy = - tensor_info['mtp']
    moment_tensor = np.zeros((3, 3))
    moment_tensor[0, :] = np.array([mxx, mxy, mxz])
    moment_tensor[1, :] = np.array([mxy, myy, myz])
    moment_tensor[2, :] = np.array([mxz, myz, mzz])
        
    mt_eigenvalue, mt_eigenvectors = np.linalg.eigh(moment_tensor)
    imax = np.argmax(mt_eigenvalue)
    imin = np.argmin(mt_eigenvalue)
    slip_vector = (mt_eigenvectors[:, imax] + mt_eigenvectors[:, imin])
    fault_normal = (mt_eigenvectors[:, imax] - mt_eigenvectors[:, imin])
    if fault_normal[2] > 0:
        fault_normal = - fault_normal
        slip_vector = - slip_vector
    strike1, dip1, rake1 =\
        __strike_dip_rake_from_ln(slip_vector, fault_normal)
    if slip_vector[2] > 0:
        fault_normal = - fault_normal
        slip_vector = - slip_vector
    strike2, dip2, rake2 =\
        __strike_dip_rake_from_ln(fault_normal, slip_vector)
    strike1 = np.remainder(strike1, 360) + stk_sens
    strike2 = np.remainder(strike2, 360) + stk_sens
    dip1 = dip1 + dip_sens
    dip2 = dip2 + dip_sens
    rake1 = np.remainder(rake1, 360)
    rake2 = np.remainder(rake2, 360)
    if rake1 > 180:
        rake1 = rake1 - 360
    if rake2 > 180:
        rake2 = rake2 - 360
    np1_plane_info = {'strike': strike1, 'dip': dip1, 'rake': rake1}
    ffm_np1_data = {'plane_info': np1_plane_info}
    np2_plane_info = {'strike': strike2, 'dip': dip2, 'rake': rake2}
    ffm_np2_data = {'plane_info': np2_plane_info}
    return ffm_np1_data, ffm_np2_data


def __strike_dip_rake_from_ln(slip_vector, fault_normal):
    """We compute strike, dip and rake from the slip vector and the vector
    which is normal to the fault plane, following Udias Fig. 16.19

    :return (strike, dip, rake): strike, dip and rake
    :type (strike, dip, rake): tuple of floats

    .. warning::
        This routine is copied from the same routine in ``instaseis`` library.
        The only modification is here we work in xyz coordinates, as opposed
        to mrt as in said routine.
    
    """
    l_norm = slip_vector / np.linalg.norm(slip_vector, 2)
    n_norm = fault_normal / np.linalg.norm(fault_normal, 2)
    delta = np.arccos(-n_norm[2])
    phi = np.arctan2(-n_norm[0], n_norm[1])

    # needs two different formulas, because the first is unstable for dip = 0
    # and the second for dip = 90
    if delta > 0.1:
        lambd = np.arctan2(
            -l_norm[2], np.sin(delta)
            * (l_norm[0] * np.cos(phi) + l_norm[1] * np.sin(phi)))
    else:
        lambd = np.arctan2(
            (l_norm[0] * np.sin(phi) - l_norm[1] * np.cos(phi)),
            np.cos(delta)
            * (l_norm[0] * np.cos(phi) + l_norm[1] * np.sin(phi)))

    strike = np.rad2deg(phi)
    dip = np.rad2deg(delta)
    rake = np.rad2deg(lambd)

    return strike, dip, rake


def __get_moment_mag(moment_tensor):
    """Given the moment tensor, we compute the seismic moment, following
    Herrmann (1989).
    """
    mt_eigenvalue = np.linalg.eigh(moment_tensor)[0]
    imax = np.argmax(mt_eigenvalue)
    amax = mt_eigenvalue[imax]
    imin = np.argmin(mt_eigenvalue)
    amin = mt_eigenvalue[imin]
    moment_mag = (np.abs(amax) + np.abs(amin)) / 2
    return moment_mag
    
    
def __is_number(string):
    """Rutina auxiliar.
    """
    try:
        float(string)
        return True
    except ValueError:
        return False


if __name__ == '__main__':
    import sys
    cmt_file = sys.argv[1]
    tensor_info = read_gcmt_file(cmt_file)
    tensor_info = modify_tensor(tensor_info)
    moment_mag = tensor_info['moment_mag']
    moment_mag = 2 * np.log10(moment_mag) / 3 - 10.7
    print(moment_mag)
    write_tensor(tensor_info)
