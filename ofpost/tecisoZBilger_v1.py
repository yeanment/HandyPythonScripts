#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This is a converter from 1D line timeseries data to Tecplot 2D data format and
get flame position (defined with maximum Qdot).
"""

import numpy as np
import scipy.interpolate, scipy.signal
import sys, os
import shutil
import argparse
import utils

def tecplotMerge(path_root, time_start, time_end, rfPosMaximal, isoSdSurface, isoSdValue, isoZBilgerValue, \
    savenpz=False, writeCombine=False, averageWdith=32, keepDataOption=0):
    
    # sort the time series
    try:
        time_series = list(filter(utils.is_digit, os.listdir(path_root)))
        time_series.sort(key=lambda item:eval(item))
    except:
        print("Failed to get time series.")
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

    # obtain the item data
    item_dataXYZ = ['X', 'Y', 'Z']
    # thermo and flow field variables
    item_dataS = ['H', 'H2', 'H2O', 'H2O2', 'HO2', 'O', 'O2', 'OH', 'Qdot', 'T', 'p', 'rho', 'ZBilger']
    item_dataV = ['U']
    item_dataVS = [item for itemV in item_dataV for item in [itemV + '_x', itemV + '_y', itemV + '_z']]
    # Sd related variables
    item_dataSdS = ['K1', 'K2', 'K3', 'K4', 'cur1', 'cur2', 'igradU', 'nngradU', \
                    'sd_conv', 'sd_corr', 'sd_diff', 'sd_rr', 'sd', 'sd_unsteady']
    item_dataSdV = ['W', 'grad', 'n']
    item_dataSdSIso = [item + '_sd' + isoSdSurface for item in item_dataSdS]
    item_dataSdVIso = [item + '_sd' + isoSdSurface for item in item_dataSdV]
    item_dataSdVSIso = [item for itemV in item_dataSdVIso \
                             for item in [itemV + 'x', itemV + 'y', itemV + 'z']]
    
    # obtain the item data name
    item_data_fileS = item_dataS + item_dataSdSIso # scalar file
    item_data_fileV = item_dataV + item_dataSdVIso # vector file
    item_data = item_dataXYZ + item_data_fileS + item_dataVS + item_dataSdVSIso

    # obtain the item data path
    iso_file_ext = '_iso' + 'ZBilger' + str(isoZBilgerValue) + '.raw'
    item_data_pathS = [item_name + iso_file_ext \
                       for index, item_name in enumerate(item_data_fileS)]
    item_data_pathV = [item_name + iso_file_ext \
                       for index, item_name in enumerate(item_data_fileV)]
    
    index_ZBilger = utils.first(item_data, condition=lambda x:x=='ZBilger')
    
    # start iterable
    if(rfPosMaximal):
        isoSdValue = [-1]
    isoSdValue.sort()
    lenSdValue = len(isoSdValue)

    dataf = np.zeros((lenSdValue, end_index - start_index, len(item_data)+1+3+2+2+2))
    iter_e = np.zeros((lenSdValue), dtype=np.int16)

    for time_index in range(start_index, end_index):
        time_now = eval(time_series[time_index])
        if((time_index - start_index)%50 == 0):
            print("Working on time {0}.".format(time_now))
        
        data_pathN = os.path.join(path_root, time_series[time_index])
        scalar_pathN = [os.path.join(data_pathN, item) \
                            for index, item in enumerate(item_data_pathS)]
        vector_pathN = [os.path.join(data_pathN, item) \
                            for index, item in enumerate(item_data_pathV)]

        # Generate out put for current time step
        try:
            index_ZBilger_path = utils.first(item_data_fileS, condition=lambda x:x=='ZBilger')
            # print(scalar_pathN[index_ZBilger_path])
            with open(scalar_pathN[index_ZBilger_path]) as f:
                ZBilgerN = np.loadtxt(f, delimiter = " ", skiprows=2)
        except BaseException as e:
            print("Failed to load ZBilger at time {0}.".format(time_now))
            continue
        
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
        
        # postprocess for dataN
        # x, y the same, average z
        if(keepDataOption):
            index = np.logical_and(np.abs(np.diff(dataN[:, 0])) < 1E-20, np.abs(np.diff(dataN[:, 1])) < 1E-20)
            splits = np.where(~index)[0] + 1
            if splits.size != 0:    
                dataN = np.stack(map(lambda item:np.mean(item, axis=0), \
                                filter(lambda item:item.shape[0] == keepDataOption, \
                                        np.split(dataN, splits, axis=0))), axis=0)
        sortIndex = np.argsort(dataN[:, 1])
        dataN = dataN[sortIndex, :]
        # determine the local tangential vector
        tck, u = scipy.interpolate.splprep(dataN[:, 0:3].transpose(), s=0)
        dx, dy, dz = scipy.interpolate.splev(u, tck,der=1)
        ds = np.sqrt(dx*dx + dy*dy + dz*dz)

        if(writeCombine):
            # write to file 
            tecdata_pathN = os.path.join(data_pathN, \
                                         'dataCombine-Time_' + time_series[time_index] + \
                                         '-iso' + 'ZBilger_{0}'.format(isoZBilgerValue) + '.dat')
            with open(tecdata_pathN, 'w') as f:
                file_header = 'TITLE="{0}_{1}_{2}"\r\n'.format(time_now, 'ZBilger', isoZBilgerValue)
                variables = 'VARIABLES="{0}" "n_x" "n_y" "n_z" \r\n'.format('" "'.join(item_data))
                zone_header = 'ZONE T="{0}_{1}_{2}" I={3} J={4}\r\n \
                               ZONETYPE = Ordered, DATAPACKING = BLOCK\r\n'.format( \
                               time_now, 'ZBilger', isoZBilgerValue, dataN.shape[0], 1)
                f.writelines([file_header, variables, zone_header])
            with open(tecdata_pathN, 'ab') as f:
                for i in range(dataN.shape[1]):
                    np.savetxt(f, dataN[:, i], delimiter=',', newline='\n')
                np.savetxt(f, dx/ds, delimiter=',', newline='\n')
                np.savetxt(f, dy/ds, delimiter=',', newline='\n')
                np.savetxt(f, dz/ds, delimiter=',', newline='\n')
        
        
        # postprocess for ZBilger
        # Special concern for -1 value
        if(rfPosMaximal):
            # find outmost local maximal
            QdotN = dataN[:, utils.first(item_data, condition=lambda x:x=='Qdot')]
            QdotNm = np.max(QdotN) # np.argmax
            curl = dataN[:, utils.first(item_data, condition=lambda x:x==('cur1_sd' + isoSdSurface))]
            index_QdotmList = np.where(\
                np.logical_and(np.diff(np.diff(curl)) < 0, QdotN[1:-1] > 0.5*QdotNm))[0]
            if index_QdotmList.size == 0:
                continue
            index_Qdotm = index_QdotmList[np.argmax(curl[index_QdotmList])]
            # index_Qdotm = index_QdotmList[-1]

            dataf[0, iter_e[0], 0] = time_now
            dataf[0, iter_e[0], 1:len(item_data)+1] = dataN[index_Qdotm, :]
            
            dataf[0, iter_e[0], len(item_data)+1] = dx[index_Qdotm]/ds[index_Qdotm]
            dataf[0, iter_e[0], len(item_data)+2] = dy[index_Qdotm]/ds[index_Qdotm]
            dataf[0, iter_e[0], len(item_data)+3] = dz[index_Qdotm]/ds[index_Qdotm]

            # determine the interpolated local velocity
            # indexColdS = utils.first(dataN[:, 1], condition=lambda x:x>0.03)
            # indexColdE = utils.first(dataN[:, 1], condition=lambda x:x>0.04)
            # A = np.vstack([dataN[indexColdS:indexColdE, 1], np.ones(indexColdE - indexColdS)]).T
            # m, c = np.linalg.lstsq(A, \
            #     dataN[indexColdS:indexColdE, utils.first(item_data, condition=lambda x:x=='U_y')])[0]
            # dataf[0, iter_e[0], len(item_data)+1] = m*dataf[0, iter_e[0], 2] + c
            
            iter_e[0] += 1
        else:
            # print("interpolating yf")
            isoSdSurface_index = utils.first(item_data, condition=lambda x:x==isoSdSurface)
            dsN = np.sqrt(np.power(np.diff(dataN[:, 0]), 2.0) + np.power(np.diff(dataN[:, 1]), 2.0))
            sN = np.append(0, np.cumsum(dsN))
           
            for indexIsoSdValue in range(lenSdValue):
                isoLine = np.zeros_like(dataN[:, isoSdSurface_index]) + isoSdValue[indexIsoSdValue]
                isoX, isoY = utils.interpolated_intercepts(sN, dataN[:, isoSdSurface_index], isoLine)
                if isoY.size == 0:
                    continue
                sf = isoX[-1]
                for isoIndex in range(1, len(item_data)+1):
                    dataf[indexIsoSdValue, iter_e[indexIsoSdValue], isoIndex] = np.interp(sf, sN, dataN[:, isoIndex-1])
                
                dataf[indexIsoSdValue, iter_e[indexIsoSdValue], len(item_data)+1] = np.interp(sf, sN, dx/ds)
                dataf[indexIsoSdValue, iter_e[indexIsoSdValue], len(item_data)+2] = np.interp(sf, sN, dy/ds)
                dataf[indexIsoSdValue, iter_e[indexIsoSdValue], len(item_data)+3] = np.interp(sf, sN, dz/ds)

                # determine the interpolated local velocity
                # indexColdS = utils.first(dataN[:, 1], condition=lambda x:x>0.03)
                # indexColdE = utils.first(dataN[:, 1], condition=lambda x:x>0.04)
                # A = np.vstack([dataN[indexColdS:indexColdE, 1], np.ones(indexColdE - indexColdS)]).T
                # m, c = np.linalg.lstsq(A, \
                #     dataN[indexColdS:indexColdE, utils.first(item_data, condition=lambda x:x=='U_y')])[0]
                # dataf[indexIsoSdValue, iter_e[indexIsoSdValue], len(item_data)+1] = m*dataf[indexIsoSdValue, iter_e[indexIsoSdValue], 2] + c
                dataf[indexIsoSdValue, iter_e[indexIsoSdValue], 0] = time_now
                iter_e[indexIsoSdValue] += 1
            
    dataf = np.nan_to_num(dataf)
    # calculate dY/dt
    for indexPost in range(0, len(iter_e)):
        # obtain the dY/dt by least square interpolation
        # postprocess for dataN
        datafN = dataf[indexPost, 0:iter_e[indexPost], :]
        # print(datafN.shape)
        
        # datafN = np.squeeze(datafN, axis=0)
        # index = np.logical_and(np.abs(np.diff(datafN[:, 1])) < 5E-6, np.abs(np.diff(datafN[:, 2])) < 5E-6)
        # splits = np.where(~index)[0] + 1
        # if splits.size != 0:
        #     dataNN = np.stack(list(\
        #         map(lambda item:np.mean(item, axis=0), np.split(datafN, splits, axis=0))), axis=0)
        # else:
        #     dataNN = datafN
        # print(dataNN.shape)
        dataNN = datafN
        # calculate dX/dt dY/dt over 
        dataNN[:, len(item_data) + 4] = utils.computedYdx(dataNN[:, 0], dataNN[:, 1], averageWidth=averageWdith)
        dataNN[:, len(item_data) + 5] = utils.computedYdx(dataNN[:, 0], dataNN[:, 2], averageWidth=averageWdith)

        # calculate filtered X and Y
        dataNN[:, len(item_data) + 6] = scipy.signal.savgol_filter(dataNN[:, 1], 21, 3)
        dataNN[:, len(item_data) + 7] = scipy.signal.savgol_filter(dataNN[:, 2], 21, 3)
        dataNN[:, len(item_data) + 8] = utils.computedYdx(dataNN[:, 0], dataNN[:, len(item_data) + 6], averageWidth=50)
        dataNN[:, len(item_data) + 9] = utils.computedYdx(dataNN[:, 0], dataNN[:, len(item_data) + 7], averageWidth=50)

        if (rfPosMaximal):
            tecflamedata = os.path.join(path_root, "S{0}-E{1}-iso{2}_{3}-rfCur1_max.dat".format(\
                time_series[start_index], time_series[end_index-1], 'ZBilger', isoZBilgerValue))
            isoSdSurface = 'cur1'
            isoSdValue = ['max']
        else:
            tecflamedata = os.path.join(path_root, "S{0}-E{1}-iso{2}_{3}-rfiso{4}_{5}.dat".format(\
                time_series[start_index], time_series[end_index-1], 'ZBilger', isoZBilgerValue, isoSdSurface, isoSdValue[indexPost]))
        
        if(savenpz):
            np.savez(tecflamedata[:-4], dataNN)
            
        with open(tecflamedata, 'w') as f:
            file_header = 'TITLE="S{0}-E{1}-iso{2}_{3}"\r\n'.format(\
                time_series[start_index], time_series[end_index-1], 'ZBilger', isoZBilgerValue)
            variables = 'VARIABLES="TIME" "{0}" "n_x" "n_y" "n_z" "dX/dt" "dY/dt" "X_" "Y_" "dX_/dt" "dY_/dt"\r\n'.format('" "'.join(item_data))
            zone_header = 'ZONE T="rfiso{0}_{1}" I={2} J={3}\r\n  \
                ZONETYPE = Ordered, DATAPACKING = BLOCK\r\n'.format( \
                isoSdSurface, isoSdValue[indexPost], dataNN.shape[0], 1)
            f.writelines([file_header, variables, zone_header])
        with open(tecflamedata, 'ab') as f:
            for i in range(dataNN.shape[1]):
                np.savetxt(f, dataNN[:, i], delimiter=',', newline='\n')

if __name__ == "__main__":
    parse = argparse.ArgumentParser( \
        description='Process some input parameters for post-processing.')
    parse.add_argument('--case', action='store', nargs=1, type=utils.check_path, \
        default = '.', help='Root path.')
    parse.add_argument('--timeStart', action='store', nargs=1, type=float, \
        default = [0.0], help='Starting time.')
    parse.add_argument('--timeEnd', action='store', nargs=1, type=float, \
        default = [-1.0], help='Ending time.')
    parse.add_argument('--rfPosMaximal', action='store', nargs=1, type=utils.str2bool, \
        default=[False], help='Evaluate rf from maximal isoSurface.')
    parse.add_argument('--isoSdSurface', action='store', nargs=1, type=str, \
        help='Iso-surfaces for Sd evaluation.')
    parse.add_argument('--isoSdSurfaceValue', action='store', nargs='*', type=float, \
        help='Iso-surfaces values in Sd evaluation.')
    parse.add_argument('--isoZBilgerValue', action='store', nargs=1, type=float, \
        default=[0.5], help='Iso-surfaces values of ZBilger.')
    
    parse.add_argument('--averageWidth', action='store', nargs=1, type=int, \
        default=[32], help='Determine whether write combined file for each step.')
    parse.add_argument('--savenpz', action='store', nargs=1, type=utils.str2bool, \
        default=[False], help='Determine whether save combined npz file.')
    
    parse.add_argument('--keepDataOption', action='store', nargs=1, type=int, \
        default=[0], help='Determine kept data point. 0 for z = 0; 1 for z = +-; 2 for all.')
    
    parse.add_argument('--writeTimeStepCombine', action='store', nargs=1, type=utils.str2bool, \
        default=[False], help='Determine whether write combined file for each step.')
    
    args = parse.parse_args()
    print(args)
    
    # Print summary of current input.
    print('*******************************')
    path_root = os.path.abspath(args.case[0])
    print("Merge timeseries data in {0} from time {1} to time {2} for iso {3} = {4}.".format(
        path_root, args.timeStart, args.timeEnd, args.isoSdSurface, args.isoSdSurfaceValue))
    tecplotMerge(path_root, args.timeStart[0], args.timeEnd[0], args.rfPosMaximal[0], \
        args.isoSdSurface[0], args.isoSdSurfaceValue, args.isoZBilgerValue[0], \
        args.savenpz[0], args.writeTimeStepCombine[0], args.averageWidth[0], args.keepDataOption[0])

    
