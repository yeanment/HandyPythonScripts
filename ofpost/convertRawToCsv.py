#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This is a converter from 1D line timeseries data to Tecplot 2D data format and
get flame position (defined with maximum Qdot).
"""

import numpy as np
import scipy.signal, scipy.interpolate
import sys, os
import shutil
import argparse
import utils

def convertCsvToRaw(pathRoot, timeStart, timeEnd):
    timeSeries = utils.obtain_sorted_timeseries(pathRoot)
    timeSeriesN = list(set(timeSeries))
    timeSeriesN.sort(key=lambda item:eval(item))

    # find time start and end index
    iTimeSeriesNS = utils.first(timeSeriesN, condition=lambda x:eval(x)>=timeStart)
    if iTimeSeriesNS == None:
        print("No data for time after {0}.".format(timeStart))
        return
    if (timeEnd == -1):
        iTimeSeriesNE = len(timeSeriesN)
    else:
        iTimeSeriesNE = utils.first(timeSeriesN, condition=lambda x:eval(x)>=timeEnd) + 1
    if (iTimeSeriesNE < iTimeSeriesNS):
        print("End time is smaller than start time.")
        return
    
    itemSuffixDict = utils.obtain_sorted_itemseries(os.path.join(pathRoot, timeSeriesN[iTimeSeriesNS]))
    # Only iterate the items contained in both directory
    itemSeries = list(itemSuffixDict.keys())

    for iTimeSeriesN in range(iTimeSeriesNS, iTimeSeriesNE):
        timeNow = eval(timeSeriesN[iTimeSeriesN])
        # Output for progress
        if((iTimeSeriesN - iTimeSeriesNS)%50 == 0):
            print("Working on time {0}.".format(timeSeriesN[iTimeSeriesN]))
        for iItemSeries in range(len(itemSeries)):
            item = itemSeries[iItemSeries]
            # Output for item
            itemSuffixSeries = list(itemSuffixDict[item])
            try:
                for iItemSuffixSeries in range(len(itemSuffixSeries)):
                    fileNow = os.path.join(pathRoot, timeSeriesN[iTimeSeriesN], item + "_" + itemSuffixSeries[iItemSuffixSeries])
                    fileToNow = os.path.join(pathRoot, timeSeriesN[iTimeSeriesN], item + "_" + itemSuffixSeries[iItemSuffixSeries][0:-2] + "csv")
                    # print(fileNow)
                    with open(fileNow, "r") as f:
                        # dataItemNameFileNow = f.readline()[0:-1].split(',')
                        dataItemFileNow = f.read()
                        dataItemFileNow = dataItemFileNow.replace('\t', ',')
                        # dataItemFileNow = np.loadtxt(f, delimiter='\t', skiprows=0)
                        with open(fileToNow, "w") as fTo:
                            if(itemSuffixSeries[iItemSuffixSeries] == "K1_sdH2O_K2_sdH2O_K3_sdH2O_K4_sdH2O_cur1_sdH2O_cur2_sdH2O_igradU_sdH2O_nngradU_sdH2O_sd_conv_sdH2O_sd_corr_sdH2O_sd_diff_sdH2O_sd_rr_sdH2O_sd_sdH2O_sd_unsteady_sdH2O.xy"):
                                fTo.writelines(["x,y,z,K1_sdH2O,K2_sdH2O,K3_sdH2O,K4_sdH2O,cur1_sdH2O,cur2_sdH2O,igradU_sdH2O,nngradU_sdH2O,sd_conv_sdH2O,sd_corr_sdH2O,sd_diff_sdH2O,sd_rr_sdH2O,sd_sdH2O,sd_unsteady_sdH2O\r\n"])
                            else:
                                fTo.writelines(["x,y,z,W_sdH2O_0,W_sdH2O_1,W_sdH2O_2,grad_sdH2O_0,grad_sdH2O_1,grad_sdH2O_2,n_sdH2O_0,n_sdH2O_1,n_sdH2O_2\r\n"])
                            fTo.write(dataItemFileNow)
                            # np.savetxt(f, dataItemFileNow, delimiter=',')
                            
            except OSError or IOError:
                print("Failed to load data at time {0}.".format(timeNow))
                continue


if __name__ == "__main__":
    
    parse = argparse.ArgumentParser( \
        description='Process some input parameters for post-processing.')
    parse.add_argument('--timeStart', action='store', nargs=1, type=float, \
        default = [0.0], help='Starting time.')
    parse.add_argument('--timeEnd', action='store', nargs=1, type=float, \
        default = [-1.0], help='Ending time.')
    parse.add_argument('--case', action='store', nargs=1, type=utils.check_path, \
        default = ['.'], help='Root path.')
    
    args = parse.parse_args()
    print(args)
    pathRoot = os.path.abspath(args.case[0])

    # Print summary of current input.
    print('*******************************')
    print("Merge timeseries data in {0} from time {1} to time {2}.".format(
        pathRoot, args.timeStart, args.timeEnd))
    convertCsvToRaw(pathRoot, args.timeStart[0], args.timeEnd[0])
    # , args.averageWidth[0])

    
