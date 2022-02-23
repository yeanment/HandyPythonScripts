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

def writeCombineCsv(pathRoot, timeStart, timeEnd, isoSdSurface):

    lineSamplePath = os.path.join(pathRoot, 'lineSample')
    lineSampleIsoPath = os.path.join(pathRoot, 'lineSampleSd' + isoSdSurface)
    timeSeries = utils.obtain_sorted_timeseries(lineSamplePath)
    timeIsoSeries = utils.obtain_sorted_timeseries(lineSampleIsoPath)
    timeSeriesN = list(set(timeSeries) & set(timeIsoSeries))
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
    
    itemSuffixDict = utils.obtain_sorted_itemseries(os.path.join(lineSamplePath, timeSeriesN[iTimeSeriesNS]))
    itemIsoSuffixDict = utils.obtain_sorted_itemseries(os.path.join(lineSampleIsoPath, timeSeriesN[iTimeSeriesNS]))
    # Only iterate the items contained in both directory
    itemSeries = list(itemSuffixDict.keys() & itemIsoSuffixDict.keys())

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
                    fileNow = os.path.join(lineSamplePath, timeSeriesN[iTimeSeriesN], item + "_" + itemSuffixSeries[iItemSuffixSeries])
                    with open(fileNow) as f:
                        dataItemNameFileNow = f.readline()[0:-1].split(',')
                        dataItemFileNow = np.loadtxt(f, delimiter=',', skiprows=0)
                    if iItemSuffixSeries == 0:
                        dataItemNow = dataItemFileNow
                        dataItemNameNow = list(dataItemNameFileNow)
                    else:
                        dataItemNow = np.concatenate((dataItemNow, dataItemFileNow[:, 3:]), axis=1)
                        dataItemNameNow.extend(dataItemNameFileNow[3:])
                itemIsoSuffixSeries = list(itemIsoSuffixDict[item])
                for iItemIsoSuffixSeries in range(len(itemIsoSuffixSeries)):
                    fileNow = os.path.join(lineSampleIsoPath, timeSeriesN[iTimeSeriesN], item + "_" + itemIsoSuffixSeries[iItemIsoSuffixSeries])
                    with open(fileNow) as f:
                        dataItemNameFileNow = f.readline()[0:-1].split(',')
                        dataItemFileNow = np.loadtxt(f, delimiter=',', skiprows=0)
                    dataItemNow = np.concatenate((dataItemNow, dataItemFileNow[:, 3:]), axis=1)
                    dataItemNameNow.extend(dataItemNameFileNow[3:])
            except OSError or IOError:
                print("Failed to load data at time {0}.".format(timeNow))
                continue
            
            # Writing of combined data for each time
            fileNow = os.path.join(lineSamplePath, timeSeriesN[iTimeSeriesN], item + "iso" + isoSdSurface +".csv")
            with open(fileNow, 'w') as f:
                f.write("{0}\r\n".format(', '.join(dataItemNameNow)))
                np.savetxt(f, dataItemNow, delimiter=',', newline='\n')

if __name__ == "__main__":
    
    parse = argparse.ArgumentParser( \
        description='Process some input parameters for post-processing.')
    parse.add_argument('--timeStart', action='store', nargs=1, type=float, \
        default = [0.0], help='Starting time.')
    parse.add_argument('--timeEnd', action='store', nargs=1, type=float, \
        default = [-1.0], help='Ending time.')
    parse.add_argument('--case', action='store', nargs=1, type=utils.check_path, \
        default = ['.'], help='Root path.')
    parse.add_argument('--isoSdSurface', action='store', nargs=1, type=str, \
        default = [], help='Iso-surfaces for Sd evaluation.')

    args = parse.parse_args()
    print(args)
    pathRoot = os.path.abspath(args.case[0])

    # Print summary of current input.
    print('*******************************')
    print("Merge timeseries data in {0} from time {1} to time {2} for iso {3}.".format(
        pathRoot, args.timeStart, args.timeEnd, args.isoSdSurface))
    writeCombineCsv(pathRoot, args.timeStart[0], args.timeEnd[0], args.isoSdSurface[0])
    # , args.averageWidth[0])

    
