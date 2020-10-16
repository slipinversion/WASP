#!/usr/bin/python
# -*- coding: utf-8 -*-
"""Routine for doing a time shift of observed data.
"""


import numpy as np
import os
import json
import glob
import get_outputs
from obspy import read
import matplotlib.pyplot as plt
import seismic_tensor as tensor


def shift_match(data_type, plot=False, method='full'):
    """We shift synthetic data to maximize cross-correlation with observed data.

    :param data_type: list of data types to be used in modelling.
    :param plot: whether to plot results of waveform shift
    :param method: method for cross-correlation
    :type data_type: list
    :type tensor_info: bool, optional
    :type method: string, optional
    """
    # print(data_type)
    if data_type == 'tele':
        json_file = 'tele_waves.json'
    if data_type == 'strong':
        json_file = 'strong_motion_waves.json'
    if data_type == 'cgps':
        json_file = 'cgps_waves.json'
    if data_type == 'surf':
        json_file = 'surf_waves.json'
    files = json.load(open(json_file))
    plot_folder = 'tele_shift' if data_type == 'tele' else 'strong_shift'\
        if data_type == 'strong' else 'surf_shift' if data_type == 'surf'\
        else 'cgps_shift'
    if not os.path.isdir(plot_folder):
        os.mkdir(plot_folder)
    synthetics_file = 'synm.tele' if data_type == 'tele' else 'synm.str'\
        if data_type == 'strong' else 'synm.str_low' if data_type == 'surf'\
        else 'synm.cgps'

    dt = float(files[0]['dt'])
    files = get_outputs.get_data_dict(files, syn_file=synthetics_file)
    for file in files:
        derivative = False if not 'derivative' in file else file['derivative']
        # if file['component'] in ['P', 'BHZ']:
        #     file['synthetic'] = []
        #     file['observed'] = []
        #     continue
        print(file['name'], file['component'])
        synt_tr = file['synthetic']
        stream = read(file['file'])
        nshift = int(5 / dt) if data_type == 'tele' else int(5 / dt)\
            if data_type == 'strong' else int(12 / dt) if data_type == 'surf'\
            else int(4 / dt)
        obser_tr = stream[0].data
        if derivative:
            obser_tr = np.diff(obser_tr)
        length = int(float(file['duration'])) if method == 'full' else int(25 / dt)
        start = int(file['start_signal'])
        print(length)
        tr_shift = _shift(obser_tr, synt_tr, nshift, length, start, plot=plot)
        file['start_signal'] = start + tr_shift

        if plot:
            name = file['name']
            component = file['component']
            synthetic = synt_tr[:length]
            fig, axes = plt.subplots(2, 1)
            fig.suptitle('{} {}'.format(name, component))
            if start >= 0:
                observed0 = np.array(
                    [val for i, val in enumerate(obser_tr[start:])\
                    if i < length])
            else:
                observed0 = np.array(
                    [val for i, val in enumerate(obser_tr) if i < length])
                observed0 = np.concatenate(([0] * -start, observed0))
            print(len(observed0), len(synthetic))
            axes[0].plot(observed0)
            axes[0].plot(synthetic, 'r')
            axes[0].set_title('Before Shift')
            start2 = start + tr_shift
            if start2 >= 0:
                observed1 = np.array(
                    [val for i, val in enumerate(obser_tr[start2:])\
                    if i < length])
            else:
                observed1 = np.array(
                    [val for i, val in enumerate(obser_tr) if i < length])
                observed1 = np.concatenate(([0] * -start2, observed1))
            axes[1].plot(observed1)
            axes[1].plot(synthetic, 'r')
            axes[1].set_title('After Shift')
            name_file = os.path.join(
                plot_folder, '{}_{}.png'.format(name, component))
            plt.savefig(name_file)
            plt.close(fig)
        file['synthetic'] = []
        file['observed'] = []

#    if do_shift:
    with open(json_file,'w') as f:
         json.dump(
             files, f, sort_keys=True, indent=4,
             separators=(',', ': '), ensure_ascii=False)
    return


def _shift(obser_tr, syn_tr, nshift, length, start_pos, plot=False):
    """Routine for shifting an observed waveform to maximize the
    cross-correlation to a synthetic waveform.
    """
    err_max = 0
    j_min = 0
#    length = min(length, len(syn_tr), 150)
    synthetic = syn_tr[:length]
    for j in range(-nshift, nshift + 1):
        # observed = np.array(
        #     [val for i, val in enumerate(obser_tr[start_pos + j:])\
        #      if i < length])
        start2 = start_pos + j
        if start2 < 0:
            observed = np.array(
                [val for i, val in enumerate(obser_tr) if i < length + start2])
            exy = np.sum(observed * synthetic[-start2:])
        else:
            observed = np.array(
                [val for i, val in enumerate(obser_tr[start_pos + j:])\
                if i < length])
            synthetic2 = synthetic[:len(observed)]
            # print(len(observed), len(synthetic2))
            exy = np.sum(observed * synthetic2)
        err = 2 * exy
#        err = np.max(np.abs(c_obs - c_syn))
#        err = err / max(np.max(np.abs(c_syn), np.max(np.abs(c_obs)))#1.0 - 2.0 * exy / (exx + eyy)

        if err_max <= err:
            err_max = err
            j_min = j
    print('Final err max:', err_max)
    return j_min


def _print_arrival(tensor_info):
    """
    """
    date_origin = tensor_info['date_origin']
    other_files = glob.glob('../data/*BHZ*sac')
    files = json.load(open('tele_waves.json'))
    plot_folder = 'tele_arrival'
    if not os.path.isdir(plot_folder):
        os.mkdir(plot_folder)
    for file in files:
        start = file['start_signal']
        sac = file['file']
        stream = read(sac)
        trace = stream[0]
        stat = trace.stats.station
        chan = trace.stats.channel
        if not chan == 'BHZ':
            continue
        other_file = [v for v in other_files if stat in v]
        other_stream = read(other_file[0])
        trace2 = other_stream[0]
        if not 't1' in trace2.stats.sac:
            continue
        arrival = trace2.stats.sac['t1']
        start = date_origin + arrival - 15
        end = date_origin + arrival + 15
        trace2.trim(starttime=start, endtime=end)
        fig = plt.figure()
        data = trace2.data
        time = np.linspace(-15, 15, len(data))
        plt.plot(time, data)
        plt.title(stat)
        plt.axvline(x=0, color='r')
#        plt.xticks(time)
        fig.savefig(os.path.join(plot_folder, '{}_pick'.format(stat)))
        del fig
    return


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--folder", default=os.getcwd(),
                        help="folder where there are input files")
    parser.add_argument('-o', '--option', choices=['match', 'manual'],
                        help='choose whether to shift by cross-correlation\
                        or plot to pick manually')
    parser.add_argument("-t", "--type", choices=['tele', 'strong', 'cgps', 'surf'],
                        help="type of data to apply shift")
    parser.add_argument("-p", "--plot", action="store_true",
                        help="plot or not results of shift")
    parser.add_argument("-gcmt", "--gcmt_tensor",
                        help="location of GCMT moment tensor file")
    args = parser.parse_args()
    os.chdir(args.folder)
    if args.gcmt_tensor:
        cmt_file = args.gcmt_tensor
        tensor_info = tensor.get_tensor(cmt_file=cmt_file)
    if args.option == 'match':
        plot = False if not args.plot else True
        shift_match(args.type, plot=plot)#, method='start')
    else:
        _print_arrival(tensor_info)

