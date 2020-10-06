# -*- coding: utf-8 -*-
"""Management of fault planes files and solution model
"""


import numpy as np
import json
import errno
import os


def __get_planes_json():
    """
    """
    if not os.path.isfile('segments_data.json'):
        raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT), 'segments_data.json')
    data = json.load(open('segments_data.json'))
    rise_time = data['rise_time']
    segments = data['segments']
    return segments, rise_time


def __read_planes_info():
    """
    """
    
    segments, rise_time = __get_planes_json()
    n_segments = len(segments)
    
    with open('Fault.time', 'r') as subfaults_file:
        lines = [line.split() for line in subfaults_file]
    
    nx_ps = int(lines[1][3])
    ny_ps = int(lines[1][4])
    
    with open('Fault.pos', 'r') as point_sources_file:
        lines2 = [line.split() for line in point_sources_file]
#
# read Fault.pos later
#
    point_sources = [[]] * n_segments
    line = 0
    for i_segment, segment in enumerate(segments):
        n_sub_y = segment['n_sub_y']
        n_sub_x = segment['n_sub_x']
        matrix = np.zeros((n_sub_y, n_sub_x, ny_ps, nx_ps, 7))
        line = line + 1
        for iys in range(n_sub_y):
            for ixs in range(n_sub_x):
                for j in range(ny_ps):
                    for i in range(nx_ps):
                        array = [float(val) for val in lines2[line]]
                        matrix[iys, ixs, j, i, :] = np.array(array)
                        line = line + 1
        point_sources[i_segment] = matrix
    return segments, rise_time, point_sources


def __unpack_plane_data(plane_info):
    """Auxiliary routine. We extract information from a fault segment.
    """
    n_sub_stk = plane_info['n_sub_x']
    n_sub_dip = plane_info['n_sub_y']
    delta_x = plane_info['delta_x']
    delta_y = plane_info['delta_y']
    hyp_stk = plane_info['hyp_stk'] - 1
    hyp_dip = plane_info['hyp_dip'] - 1
    return n_sub_stk, n_sub_dip, delta_x, delta_y, hyp_stk, hyp_dip


if __name__ == '__main__':
    import os
#    segments, rise_time, point_sources = __read_planes_info()
#    print(segments)
    os.chdir('/home/pk/Inversion_Chen_Ji/Inversiones/20140401234647/paper_challa/pk.5/NP2/srcmod_solutions/')
    tensor_info, solution = read_solution_fsp_format('hayes')
    print(solution)