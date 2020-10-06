#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Module to find and store modelling parameters, such as number of iterations
of simmulated annealing and its cooling rate, weight of regularization
constraints, and boundaries of model space.
"""


import json
import management as mng
import plane_management as pl_mng
import argparse
import os
import seismic_tensor as tensor


def modelling_prop(tensor_info, data_type=None, moment_mag=None):
    r"""We write a file which is a very important input for performing
    FFM modelling using simmulated annealing.

    :param moment_mag: seismic moment to be used in modelling
    :param data_type: list of data types to be used in modelling.
    :param tensor_info: dictionary with moment tensor information
    :type moment_mag: float, optional
    :type data_type: list, optional
    :type tensor_info: dict
    """
    moment_mag = tensor_info['moment_mag'] if not moment_mag else moment_mag
    time_shift = tensor_info['time_shift']
    if not data_type:
        data_type = mng.update_data(tensor_info)
    data_type2 = set(data_type) - {'gps'}

    syn_len = int(2.5 * time_shift)
    factor = 1 / len(data_type2)
    moment_weight = 0.01 * factor
    slip_weight = 0.1 * factor
    time_weight = 0.1 * factor
    dictionary = {
            'initial_temperature': 0.01,
            'seismic_moment': moment_mag,
            'moment_weight': moment_weight,
            'slip_weight': slip_weight,
            'time_weight': time_weight,
            'max_source_dur': syn_len,
            'iterations': 250,
            'cooling_rate': 0.971
    }

    with open('annealing_prop.json','w') as f:
         json.dump(
                 dictionary, f, sort_keys=True, indent=4,
                 separators=(',', ': '), ensure_ascii=False)

    moment_mag = tensor_info['moment_mag']
    segments, rise_time = pl_mng.__get_planes_json()
    if moment_mag < 8 * 10**25:
        dmax = 500.0
    elif moment_mag < 10 ** 28:
        dmax = 1500.0
    elif moment_mag < 7 * 10**28:
        dmax = 2000.0
    elif moment_mag < 4 * 10**29:
        dmax = 4000.0
    else:
        dmax = 8000.0
    step = 20
    nstep = min(int(dmax / step) + 1, 51)

    dictionary2 = {
            'min_slip': 0,
            'max_slip': dmax,
            'slip_step': nstep
    }

    segments2 = []
    for segment in segments:
        rake = segment['rake']
        rake_min = rake - 20
        rake_max = rake + 20
        rstep = 21
        if rake_min < 0:
            rake_min = rake_min + 360
            rake_max = rake_max + 360
        dictionary3 = {
                'rake_min': rake_min,
                'rake_max': rake_max,
                'rake_step': rstep
        }
        dictionary3.update(dictionary2)
        segments2 = segments2 + [dictionary3]

    with open('model_space.json','w') as f:
         json.dump(
                 segments2, f, sort_keys=True, indent=4,
                 separators=(',', ': '), ensure_ascii=False)
    return dictionary, segments2


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--folder", default=os.getcwd(),
                        help="folder where there are input files")
    parser.add_argument("-gcmt", "--gcmt_tensor",
                        help="location of GCMT moment tensor file")
    parser.add_argument("-t", "--tele", action="store_true",
                        help="compute files with teleseismic data")
    parser.add_argument("-su", "--surface", action="store_true",
                        help="compute files with surface waves data")
    parser.add_argument("-st", "--strong", action="store_true",
                        help="compute files with strong motion data")
    parser.add_argument("--cgps", action="store_true",
                        help="compute files with cGPS data")
    args = parser.parse_args()

    os.chdir(args.folder)
    if args.gcmt_tensor:
        cmt_file = args.gcmt_tensor
        tensor_info = tensor.get_tensor(cmt_file=cmt_file)
    else:
        tensor_info = tensor.get_tensor()
    data_type = []
    if args.tele:
        data_type = data_type + ['tele_body']
    if args.surface:
        data_type = data_type + ['surf_tele']
    if args.strong:
        data_type = data_type + ['strong_motion']
    if args.cgps:
        data_type = data_type + ['cgps']
    modelling_prop(tensor_info, data_type=data_type)
