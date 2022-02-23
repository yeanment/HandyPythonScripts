#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This is a converter from 1D line timeseries data to Tecplot 2D data format and
get flame position (defined with maximum Qdot).
"""

import numpy as np
import scipy.signal
import sys, os
import shutil
import argparse
import utils

def tecplotMerge(path_root, time_start, time_end, item, isoSdSurface, isoSdSurfaceValue, \
    writeCombine=False, averageWdith=64):
    try:
        time_series = list(filter(utils.is_digit, os.listdir(os.path.join(path_root, 'lineSample'))))
        time_series.sort(key=lambda item:eval(item))
    except:
        print("Failed to get time series ")
        return

    # find time start and end index
    start_index = utils.first(time_series, condition=lambda x:eval(x)>=time_start)
    if start_index == None:
        print("No data for time after {0}.".format(time_start))
        return
    if (time_end == -1):
        end_index = len(time_series)
    else:
        end_index = utils.first(time_series, condition=lambda x:eval(x)>=time_end) + 1
    if (end_index < start_index):
        print("End time is smaller than start time.")
        return
    

    item_dataXYZ = ['X', 'Y', 'Z']
    item_dataS = ['H', 'H2', 'H2O', 'H2O2', 'HO2', 'N2', 'O', 'O2', 'OH', 'Qdot', 'T', 'p', 'rho']
    # item_dataS = ['p', 'Qdot', 'T', 'rho', 'H', 'O', 'OH', 'HO2', 'H2O2', 'H2', 'O2', 'H2O', 'N2']
    item_dataV = ['U']
    item_dataVS = [item for itemV in item_dataV for item in [itemV + '_x', itemV + '_y', itemV + '_z']]
    
    item_dataSdS = ['K1', 'K2', 'K3', 'K4', 'cur1', 'cur2', 'igradU', 'nngradU', \
                'sd_conv', 'sd_corr', 'sd_diff', 'sd_rr', 'sd', 'sd_unsteady']
    item_dataSdV = ['W', 'grad', 'n']
    item_dataSdVS = [item for itemV in item_dataSdV for item in [itemV + 'x', itemV + 'y', itemV + 'z']]
    item_dataSdSIso = [item + '_sd' + isoSdSurface for item in item_dataSdS]
    item_dataSdVIso = [item + '_sd' + isoSdSurface for item in item_dataSdV]
    item_dataSdVSIso = [item + '_sd' + isoSdSurface for item in item_dataSdVS]
    item_dataSdIso = item_dataSdSIso + item_dataSdVSIso

    # obtain the item data name
    item_data = item_dataXYZ + item_dataS + item_dataVS + item_dataSdIso

    isoSdSurface_index = utils.first(item_data, condition=lambda x:x==isoSdSurface)
    
    item_data_path = [os.path.join(path_root, item_name) for index, item_name in enumerate(item_data)]
    time_data_path = os.path.join(path_root, 'time')
    utils.item_file_remove(item_data_path)
    utils.item_file_remove([time_data_path])

    # start iterable
    lenSdValue = len(isoSdSurfaceValue)
    if lenSdValue == 0:
        # default Qdot max version
        # interpolate the local flow velocity
        dataf = np.zeros((1, end_index - start_index, len(item_data)+4))
        yf_ = np.zeros((1, end_index - start_index, 1))
        iter_e = np.zeros((1), dtype=np.int16)
        # lenSdValue = 1
        # isoSdSurfaceValue = [-1]
    else:
        dataf = np.zeros((lenSdValue, end_index - start_index, len(item_data)+4))
        yf_ = np.zeros((lenSdValue, end_index - start_index, 1))
        iter_e = np.zeros((lenSdValue), dtype=np.int16)

    for time_index in range(start_index, end_index):
        time_now = eval(time_series[time_index])
        if((time_index - start_index)%50 == 0):
            print("Working on time {0}.".format(time_now))
        
        path_data1 = os.path.join(os.path.join(path_root, 'lineSample'), time_series[time_index])
        path_data2 = os.path.join(os.path.join(path_root, 'lineSampleSd'+isoSdSurface), time_series[time_index])

        # file_list = os.listdir(path_data)
        file_scalar1 = "{0}_{1}.csv".format(item, '_'.join(item_dataS))
        file_vector1 = "{0}_{1}.csv".format(item, '_'.join(item_dataV))
        
        file_scalar2 = "{0}_{1}.csv".format(item, '_'.join(item_dataSdSIso))
        file_vector2 = "{0}_{1}.csv".format(item, '_'.join(item_dataSdVIso))

        scalar1_path = os.path.join(path_data1, file_scalar1)
        vector1_path = os.path.join(path_data1, file_vector1)
        scalar2_path = os.path.join(path_data2, file_scalar2)
        vector2_path = os.path.join(path_data2, file_vector2)

        if(not (os.path.exists(scalar1_path) and os.path.exists(scalar2_path) \
            and os.path.exists(vector1_path) and os.path.exists(vector2_path))):
            print("Error in get file.")
            continue

        with open( scalar1_path) as f:
            scalar1_time = np.loadtxt(f,delimiter = ",", skiprows = 1)
        with open( vector1_path) as f:
            vector1_time = np.loadtxt(f,delimiter = ",", skiprows = 1)

        with open( scalar2_path) as f:
            scalar2_time = np.loadtxt(f,delimiter = ",", skiprows = 1)
        with open( vector2_path) as f:
            vector2_time = np.loadtxt(f,delimiter = ",", skiprows = 1)

        data1_time = np.hstack((scalar1_time, vector1_time[:, 3:]))
        data2_time = np.hstack((scalar2_time[:,3:], vector2_time[:, 3:]))
        # print(data1_time.shape, data2_time.shape)
        # print(time_now)
        data_time = np.hstack((data1_time[:data2_time.shape[0], :], data2_time))
        
        nx = data_time.shape[0]
        
        # with open(time_data_path, 'a') as f:
        #     f.write("{0}*{1}\r\n".format(nx, time_now))
        # for i in range(len(item_data_path)):
        #     with open(item_data_path[i], 'a') as f:
        #         np.savetxt(f, data_time[:, i], delimiter=',', newline='\n')
                # f.write("\n\r")
        
        # postprocess 
        if lenSdValue == 0:
            iflame = np.argmax(data_time[:, utils.first(item_data, condition=lambda x:x=='Qdot')])
            dataf[0, iter_e[0], 0] = time_now
            dataf[0, iter_e[0], 1:len(item_data)+1] = data_time[iflame, 0:len(item_data)]
            iter_e[0] += 1
        else:
            dsN = np.sqrt(np.power(np.diff(data_time[:, 0]), 2.0) + np.power(np.diff(data_time[:, 1]), 2.0))
            sN = np.append(0, np.cumsum(dsN))
           
            # interpolating yf
            for indexIsoSdValue in range(lenSdValue):
                isoLine = np.zeros_like(data_time[:, isoSdSurface_index]) + isoSdSurfaceValue[indexIsoSdValue]
                # rLine = np.sqrt((data_time[:, 0] - data_time[0, 0]) * (data_time[:, 0] - data_time[0, 0]) \
                #     + (data_time[:, 1] - data_time[0, 1]) * (data_time[:, 1] - data_time[0, 1]))
                # # rLine = data_time[:, 1]
                try:
                    isoX, isoY = utils.interpolated_intercepts(sN, data_time[:, isoSdSurface_index], isoLine)
                except:
                    continue
                if isoY.size == 0:
                    continue
                sf = isoX[-1]
                
                if(iter_e[indexIsoSdValue] > 6):
                    sf_a = np.average(yf_[indexIsoSdValue, iter_e[indexIsoSdValue]-5:iter_e[indexIsoSdValue]])
                    sf = isoX[np.argmin(np.abs(sf - sf_a))]
                yf_[indexIsoSdValue, iter_e[indexIsoSdValue]] = sf

                qdot = data_time[:, utils.first(item_data, condition=lambda x:x=='Qdot')]
                qdot_index = np.where(qdot > 0.01)[0]
                if(len(qdot_index) < 2):
                    # consider all
                    for isoIndex in range(1, len(item_data)+1):
                        dataf[indexIsoSdValue, iter_e[indexIsoSdValue], isoIndex] = np.interp(sf, sN, data_time[:, isoIndex-1])
                else:
                    for isoIndex in range(1, len(item_data)+1):
                        if (qdot_index[-1] - qdot_index[0])%2 == 1:
                            window_length = min(25, qdot_index[-1] - qdot_index[0])
                        else:
                            window_length = min(25, qdot_index[-1] - qdot_index[0] - 1)
                        filterData = scipy.signal.savgol_filter(data_time[qdot_index[0]:qdot_index[-1], isoIndex-1], window_length, 3)
                        dataf[indexIsoSdValue, iter_e[indexIsoSdValue], isoIndex] = np.interp(sf, sN[qdot_index[0]:qdot_index[-1]], filterData)

                dataf[indexIsoSdValue, iter_e[indexIsoSdValue], 0] = time_now

                qdot = data_time[:, utils.first(item_data, condition=lambda x:x=='Qdot')]
                qdot_index = np.where(qdot > np.max(qdot)/2)[0]
                if(len(qdot_index) < 2):
                    dataf[indexIsoSdValue, iter_e[indexIsoSdValue], len(item_data)+1] = 0
                else:
                    dx_f = data_time[qdot_index[-1], 0] - data_time[qdot_index[0], 0]
                    dy_f = data_time[qdot_index[-1], 1] - data_time[qdot_index[0], 1]
                    dz_f = data_time[qdot_index[-1], 2] - data_time[qdot_index[0], 2]
                    dataf[indexIsoSdValue, iter_e[indexIsoSdValue], len(item_data)+1] = np.sqrt(dx_f*dx_f + dy_f*dy_f + dz_f*dz_f)
                # using maximum temperature gradient
                # using grad

                gradIsoIndex = utils.first(item_data, condition=lambda x:x=='gradx_sd' + isoSdSurface)
                nIsoIndex = utils.first(item_data, condition=lambda x:x=='nx_sd' + isoSdSurface)
                grad = data_time[:, gradIsoIndex]*data_time[:, nIsoIndex] + \
                    data_time[:, gradIsoIndex+1]*data_time[:, nIsoIndex+1] + \
                    data_time[:, gradIsoIndex+2]*data_time[:, nIsoIndex+2]
                dataf[indexIsoSdValue, iter_e[indexIsoSdValue], len(item_data)+2] = \
                    (np.max(data_time[:, isoSdSurface_index]) - np.min(data_time[:, isoSdSurface_index]))/np.max(np.abs(grad))
                iter_e[indexIsoSdValue] += 1

    dataf = np.nan_to_num(dataf)
    print(dataf.shape)
    # dataf = dataf[utils.first(dataf[:,2], condition=lambda x:x!=0)+1:iter_e, :]
    for indexIsoSdValue in range(0, len(iter_e)):
        datafN = dataf[indexIsoSdValue, 0:iter_e[indexIsoSdValue], :]
        print(datafN.shape)
        # Calculate dY/dt 
        index = np.logical_and(np.abs(np.diff(datafN[:, 1])) < 1E-20, np.abs(np.diff(datafN[:, 2])) < 1E-20)
        splits = np.where(~index)[0] + 1
        if splits.size != 0:
            dataNN = np.stack(list(\
                map(lambda item:np.mean(item, axis=0), np.split(datafN, splits, axis=0))), axis=0)
        else:
            dataNN = datafN
        print(dataNN.shape)
        # calculate dY/dt over 
        for indexS in range(averageWdith, dataNN.shape[0] - averageWdith - 1, 1):
            dataX = dataNN[indexS-averageWdith:indexS+averageWdith+1, 0]
            dataY = dataNN[indexS-averageWdith:indexS+averageWdith+1, 2]
            dataNN[indexS, len(item_data) + 3] = \
                (np.mean(dataX*dataY) - np.mean(dataY)*np.mean(dataX))/(np.mean(dataX*dataX) - np.mean(dataX)*np.mean(dataX))
        
        if lenSdValue == 0:
            isoSdSurface = 'Qdot'
            isoSdSurfaceValue = ['max']
            tecflamedata = os.path.join(path_root, "{0}_S{1}_E{2}_iso{3}_{4}_rfQdotm.dat".format(item, time_series[start_index], time_series[end_index-1], isoSdSurface, isoSdSurfaceValue[0]))
        else:
            tecflamedata = os.path.join(path_root, "{0}_S{1}_E{2}_iso{3}_{4}_rfQdotm.dat".format( \
                item, time_series[start_index], time_series[end_index-1], isoSdSurface, isoSdSurfaceValue[indexIsoSdValue]))
        try:
            with open(tecflamedata, 'w') as f:
                file_header = 'TITLE="{0}_{1}_{2}_rf"\r\n'.format(item, time_series[start_index], time_series[end_index-1])
                variables = 'VARIABLES="TIME" "{0}" "l_TQ" "l_GRAD" "dY/dt"\r\n'.format('" "'.join(item_data))
                zone_header = 'ZONE T="{0}_{1}_{2}_rf" I={3} J={4}\r\n  \
                            ZONETYPE = Ordered, DATAPACKING = BLOCK\r\n'.format(item, \
                            time_series[start_index], time_series[end_index-1], \
                            dataNN.shape[0], 1)
                f.writelines([file_header, variables, zone_header])
                for i in range(dataNN.shape[1]):
                    np.savetxt(f, dataNN[:, i], delimiter=',', newline='\n')
        except:
            continue
    # finish time iter
    # combine into one file
    # tecdata = os.path.join(path_root, "{0}_{1}_{2}.dat".format(item, time_series[start_index], time_series[end_index-1]))
    # with open(tecdata, 'w') as f:
    #     file_header = 'TITLE="{0}_{1}_{2}"\r\n'.format(item, time_series[start_index], time_series[end_index-1])
    #     variables = 'VARIABLES="TIME" "{0}"\r\n'.format('" "'.join(item_data))
    #     zone_header = 'ZONE T="{0}_{1}_{2}" I={3} J={4}\r\n  \
    #                 ZONETYPE = Ordered, DATAPACKING = BLOCK\r\n'.format(item, \
    #                 time_series[start_index], time_series[end_index-1], \
    #                 nx, iter_t)
    #     f.writelines([file_header, variables, zone_header])
    #     f.write(open(time_data_path, 'r').read())
        
    # for i in range(len(item_data_path)):
    #     with open(tecdata, 'a') as f:
    #         f.write(open(item_data_path[i], 'r').read())
    # utils.item_file_remove(item_data_path)
    # utils.item_file_remove([time_data_path])


if __name__ == "__main__":
    
    parse = argparse.ArgumentParser( \
        description='Process some input parameters for post-processing.')
    parse.add_argument('--timeStart', action='store', nargs=1, type=float, \
        default = [0.0], help='Starting time.')
    parse.add_argument('--timeEnd', action='store', nargs=1, type=float, \
        default = [-1.0], help='Ending time.')
    parse.add_argument('--case', action='store', nargs=1, type=utils.check_path, \
        default = ['.'], help='Root path.')
    parse.add_argument('--item', action='store', nargs=1, type=str, \
        default = [], help='Interpolated items for post processing.')
    parse.add_argument('--isoSdSurface', action='store', nargs=1, type=str, \
        default = [], help='Iso-surfaces for Sd evaluation.')
    parse.add_argument('--isoSdSurfaceValue', action='store', nargs='*', type=float, \
        help='Iso-surfaces values in Sd evaluation.')
    # parse.add_argument('--averageWidth', action='store', nargs=1, type=int, \
    #     default=[32], help='Determine whether write combined file for each step.')
    
    parse.add_argument('--writeTimeStepCombine', action='store', nargs=1, type=bool, \
        default=[False], help='Determine whether write combined file for each step.')
    
    args = parse.parse_args()
    print(args)
    path_root = os.path.abspath(args.case[0])

    # Print summary of current input.
    print('*******************************')
    print("Merge timeseries data in {0} from time {1} to time {2} for iso {3} = {4}.".format(
        path_root, args.timeStart, args.timeEnd, args.isoSdSurface, args.isoSdSurfaceValue))
    tecplotMerge(path_root, args.timeStart[0], args.timeEnd[0], args.item[0], \
        args.isoSdSurface[0], args.isoSdSurfaceValue, \
        args.writeTimeStepCombine[0])
    # , args.averageWidth[0])

    
