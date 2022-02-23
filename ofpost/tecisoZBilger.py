#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This is a converter from 1D line timeseries data to Tecplot 2D data format and
get flame position (defined with maximum Qdot).
"""

import numpy as np
import sys, os
import shutil
# import argparse


def is_digit(string):
    try:
        float(string)
    except:
        return False
    else:
        return True

def first(iterable, condition = lambda x: True):
    try:
        return next(i for i, x in enumerate(iterable) if condition(x))
    except StopIteration:
        return None

def item_file_remove(item_path):
    for item in item_path:
        try:
            os.remove(item)
        except FileNotFoundError as e:
            pass
        else:
            pass

def interpolated_intercepts(x, y1, y2):
    """Find the intercepts of two curves, given by the same x data"""
    def intercept(point1, point2, point3, point4):
        """
        find the intersection between two lines
        the first line is defined by the line between point1 and point2
        the first line is defined by the line between point3 and point4
        each point is an (x,y) tuple.
        Returns: the intercept, in (x,y) format
        """    
        def line(p1, p2):
            A = (p1[1] - p2[1])
            B = (p2[0] - p1[0])
            C = (p1[0]*p2[1] - p2[0]*p1[1])
            return A, B, -C
        def intersection(L1, L2):
            D  = L1[0] * L2[1] - L1[1] * L2[0]
            Dx = L1[2] * L2[1] - L1[1] * L2[2]
            Dy = L1[0] * L2[2] - L1[2] * L2[0]
            if D == 0:
                raise ZeroDivisionError 
            x = Dx / D
            y = Dy / D
            return x, y
        L1 = line([point1[0],point1[1]], [point2[0],point2[1]])
        L2 = line([point3[0],point3[1]], [point4[0],point4[1]])
        R = intersection(L1, L2)
        return R

    idxs = np.argwhere(np.diff(np.sign(y1 - y2)) != 0)
    xcs = []
    ycs = []
    for idx in idxs:
        xc, yc = intercept((x[idx], y1[idx]),((x[idx+1], y1[idx+1])), ((x[idx], y2[idx])), ((x[idx+1], y2[idx+1])))
        xcs.append(xc)
        ycs.append(yc)
    return np.array(xcs), np.array(ycs)

def tecplotMerge(path_root, time_start, time_end, isoSurface, isoValue):
    try:
        time_series = list(filter(is_digit, os.listdir(path_root)))
        time_series.sort(key=lambda item:eval(item))
    except:
        print("Failed to get time series ")
        return
    else:
        pass

    # find time start and end index
    start_index = first(time_series, condition=lambda x:eval(x)>=time_start)
    if start_index == None:
        print("No data for time after {0}.".format(time_start))
        return
    if (time_end == -1):
        end_index = len(time_series)
    else:
        end_index = first(time_series, condition=lambda x:eval(x)>=time_end) + 1
    if (end_index < start_index):
        print("End time is smaller than start time.")
        return
    
    item_dataXYZ = ['X', 'Y', 'Z']
    # item_dataS = ['p', 'Qdot', 'T', 'rho', 'H', 'O', 'OH', 'HO2', 'H2O2', 'H2', 'O2', 'H2O', \
    #             'X_H', 'X_H2', 'X_H2O', 'X_H2O2', 'X_HO2', 'X_N2', 'X_O', 'X_O2', 'X_OH']
    # item_dataS = ['H', 'H2', 'H2O', 'H2O2', 'HO2', 'O', 'O2', 'OH', 'Qdot', 'T', \
    #             'X_H', 'X_H2', 'X_H2O', 'X_H2O2', 'X_HO2', 'X_N2', 'X_O', 'X_O2', 'X_OH', 'p', 'rho']
    
    item_dataS = ['H', 'H2', 'H2O', 'H2O2', 'HO2', 'O', 'O2', 'OH', 'Qdot', 'T', 'p', 'rho', 'ZBilger']
    item_dataV = ['U']
    item_dataVS = [item for itemV in item_dataV for item in [itemV + '_x', itemV + '_y', itemV + '_z']]
    
    item_dataSdS = ['K1', 'K2', 'K3', 'K4', 'cur1', 'cur2', 'igradU', 'nngradU', \
        'sd_conv', 'sd_corr', 'sd_diff', 'sd_rr', 'sd', 'sd_unsteady']
    item_dataSdV = ['W', 'grad', 'n']
    item_dataSdSIso = [item + '_sd' + isoSurface for item in item_dataSdS]
    item_dataSdVIso = [item + '_sd' + isoSurface for item in item_dataSdV]
    item_dataSdVSIso = [item for itemV in item_dataSdVIso \
        for item in [itemV + 'x', itemV + 'y', itemV + 'z']]
        
    item_data_fileS = item_dataS + item_dataSdSIso # scalar file
    item_data_fileV = item_dataV + item_dataSdVIso # vector file

    item_data = item_dataXYZ + item_data_fileS + item_dataVS + item_dataSdVSIso

    iso_file_ext = '_iso' + 'ZBilger' + isoValue + '.raw'
    item_data_pathS = [item_name + iso_file_ext \
        for index, item_name in enumerate(item_data_fileS)]
    item_data_pathV = [item_name + iso_file_ext \
        for index, item_name in enumerate(item_data_fileV)]
    

    index_ZBilger = first(item_data, condition=lambda x:x=='ZBilger')

    # start iterable
    dataf = np.zeros((end_index - start_index, len(item_data)+1))
    yf_ = np.zeros((end_index - start_index, 1))
    iter_t = 0
    iter_e = 0
    for time_index in range(start_index, end_index):
        iter_t += 1
        time_now = eval(time_series[time_index])
        if(iter_t%50 == 0):
            print("Working on time {0}.".format(time_now))
        
        data_pathN = os.path.join(path_root, time_series[time_index])
        scalar_pathN = [os.path.join(data_pathN, item) \
            for index, item in enumerate(item_data_pathS)]
        vector_pathN = [os.path.join(data_pathN, item) \
            for index, item in enumerate(item_data_pathV)]

        # Generate out put for current time step
        try:
            index_ZBilger_path = first(item_data_fileS, condition=lambda x:x=='ZBilger')
            # print(scalar_pathN[index_ZBilger_path])
            with open(scalar_pathN[index_ZBilger_path]) as f:
                ZBilgerN = np.loadtxt(f, delimiter = " ", skiprows=2)
        except BaseException as e:
            print("Failed to load ZBilger at time {0}.".format(time_now))
            continue
        else:
            pass
        try:
            datalenN = ZBilgerN.shape[0]
            dataN = np.zeros((datalenN, len(item_data)))
            dataN[:, 0:3] = ZBilgerN[:, 0:3]
            for index, item_scalar_pathN in enumerate(scalar_pathN):
                with open(item_scalar_pathN) as f:
                    scalar_dataN = np.loadtxt(f, delimiter = " ", skiprows=2)
                dataN[:, 3+index] = scalar_dataN[:, 3]
            for index, item_vector_pathN in enumerate(vector_pathN):
                with open(item_vector_pathN) as f:
                    vector_dataN = np.loadtxt(f, delimiter = " ", skiprows=2)
                data_indexS = 3 + len(scalar_pathN) + index *3
                dataN[:, data_indexS:data_indexS+3] = vector_dataN[:, 3:6]
        except BaseException as e:
            print("Failed to load data at time {0}.".format(time_now))
            continue
        else:
            pass
        # postprocess for dataN
        index = np.logical_and(np.abs(np.diff(dataN[:, 0])) < 1E-8, np.abs(np.diff(dataN[:, 1])) < 1E-8)
        indexS0 = np.append(True, ~index)
        indexS1 = np.append(~index, True)
        dataNN = (dataN[indexS0, :] + dataN[indexS1, :])/2.0
        dataN = dataNN

        # write to file 
        tecdataN_pathN = os.path.join(data_pathN, 'data' + '_iso' + 'ZBilger' + isoValue + '.dat')
        with open(tecdataN_pathN, 'w') as f:
            file_header = 'TITLE="{0}_{1}_{2}"\r\n'.format(time_now, 'ZBilger', isoValue)
            variables = 'VARIABLES="{0}"\r\n'.format('" "'.join(item_data))
            zone_header = 'ZONE T="{0}_{1}_{2}" I={3} J={4}\r\n  \
                ZONETYPE = Ordered, DATAPACKING = BLOCK\r\n'.format(time_now, 'ZBilger', isoValue, \
                dataN.shape[0], 1)
            f.writelines([file_header, variables, zone_header])
        with open(tecdataN_pathN, 'ab') as f:
            for i in range(dataN.shape[1]):
                np.savetxt(f, dataN[:, i], delimiter=',', newline='\n')
    
        # postprocess for ZBilger
        index_Qdotm = np.argmax(dataN[:, first(item_data_fileS, condition=lambda x:x=='Qdot')])
        dataf[iter_e, 0] = time_now
        dataf[iter_e, 1:] = dataN[index_Qdotm, :]
        iter_e += 1
        
    dataf = np.nan_to_num(dataf)    
    dataf = dataf[0:iter_e, :]
    tecflamedata = os.path.join(path_root, "Qdotmax-S{0}-E{1}-iso{2}_{3}.dat".format(\
        time_series[start_index], time_series[end_index-1], 'ZBilger', isoValue))
    with open(tecflamedata, 'w') as f:
        file_header = 'TITLE="Qdotmax-S{0}-E{1}-iso{2}_{3}"\r\n'.format(\
            time_series[start_index], time_series[end_index-1], 'ZBilger', isoValue)
        variables = 'VARIABLES="TIME" "{0}"\r\n'.format('" "'.join(item_data))
        zone_header = 'ZONE T="Qdotmax-iso{0}_{1}" I={3} J={4}\r\n  \
            ZONETYPE = Ordered, DATAPACKING = BLOCK\r\n'.format( \
            'ZBilger', isoValue, '', dataf.shape[0], 1)
        f.writelines([file_header, variables, zone_header])
    with open(tecflamedata, 'ab') as f:
        for i in range(dataf.shape[1]):
            np.savetxt(f, dataf[:, i], delimiter=',', newline='\n')

if __name__ == "__main__":
    # parse = argparse.ArgumentParser()
    # parse.add_argument('--timeStart', action=store, nargs=1, type=float)
    # args = parse.parse_args()
    
    if '--timeS' in sys.argv:
        time_string_index = sys.argv.index('--timeS') + 1
        try:
            time_string_start = sys.argv[time_string_index]
            time_start = float(time_string_start)
        except BaseException as e:
            print('Error in input start time')
            exit()
        else:
            pass
    else:
        print("Please reinput start time.")
        exit()

    if '--timeE' in sys.argv:
        time_string_index = sys.argv.index('--timeE') + 1
        try:
            time_string_end = sys.argv[time_string_index]
            time_end = float(time_string_end)
        except BaseException as e:
            print('Error in input end time')
            exit()
        else:
            pass
    else:
        time_end = -1

    if '--case' in sys.argv:
        path_root_index = sys.argv.index('--case') + 1
        try:
            path_root = sys.argv[path_root_index]    
        except IndexError:
            path_root = '.'
        else:
            pass
    else:
        path_root = '.'    
    path_root = os.path.abspath(path_root)
    if (os.access(path_root, os.F_OK) == False or os.path.isdir(path_root) == False):
        print('Path not exists.')
        exit()
        
    if '--isoSurface' in sys.argv:
        isoSurface_index = sys.argv.index('--isoSurface') + 1
        try:
            isoSurface = sys.argv[isoSurface_index]
            # print(type(isoValue))
        except BaseException as e:
            print('Error in input isoSurface.')
            exit()
        else:
            pass
    else:
        exit()


    if '--isoValue' in sys.argv:
        isoValue_index = sys.argv.index('--isoValue') + 1
        try:
            isoValue = sys.argv[isoValue_index]
        except BaseException as e:
            print('Error in input isoValue.')
            exit()
        else:
            pass
    else:
        exit()
    

    # Print summary of current input.
    print('*******************************')
    print("Merge timeseries data in {0} from time {1} to time {2} for iso {3} = {4}.".format(
        path_root, time_start, time_end, isoSurface, isoValue))
    tecplotMerge(path_root, time_start, time_end, isoSurface, isoValue)

    
