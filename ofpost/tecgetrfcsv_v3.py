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

def tecplotMerge(pathRoot, timeStart, timeEnd, isoSdSurface, isoSdSurfaceValue,\
    writeCombine=False, leftMaxPos=False,
    averageWdith=64, enableSavgolFilter=False, savgolFilterWidth=25):

    lineSamplePath = os.path.join(pathRoot, 'lineSample')
    lineSampleIsoPath = os.path.join(pathRoot, 'lineSampleSd' + isoSdSurface)

    timeSeries = utils.obtain_sorted_timeseries(lineSamplePath)
    timeIsoSeries = utils.obtain_sorted_timeseries(lineSampleIsoPath)
    timeSeriesN = list(set(timeSeries) & set(timeIsoSeries))
    timeSeriesN.sort(key=lambda item: eval(item))

    # find time start and end index
    iTimeSeriesNS = utils.first(timeSeriesN,
                                condition=lambda x: eval(x) >= timeStart)
    if iTimeSeriesNS == None:
        print("No data for time after {0}.".format(timeStart))
        return
    if (timeEnd == -1):
        iTimeSeriesNE = len(timeSeriesN)
    else:
        iTimeSeriesNE = utils.first(timeSeriesN,
                                    condition=lambda x: eval(x) >= timeEnd) + 1
    if (iTimeSeriesNE < iTimeSeriesNS):
        print("End time is smaller than start time.")
        return

    itemSuffixDict = utils.obtain_sorted_itemseries(
        os.path.join(lineSamplePath, timeSeriesN[iTimeSeriesNS]))
    itemIsoSuffixDict = utils.obtain_sorted_itemseries(
        os.path.join(lineSampleIsoPath, timeSeriesN[iTimeSeriesNS]))
    # Only iterate the items contained in both directory
    itemSeries = list(itemSuffixDict.keys() & itemIsoSuffixDict.keys())

    # Proceeding for different items
    for iItemSeries in range(len(itemSeries)):
        item = itemSeries[iItemSeries]
        itemSuffixSeries = list(itemSuffixDict[item])
        itemIsoSuffixSeries = list(itemIsoSuffixDict[item])

        if (writeCombine):
            if (not os.path.exists(
                    os.path.join(pathRoot, item + "-iso" + isoSdSurface))):
                os.mkdir(os.path.join(pathRoot, item + "-iso" + isoSdSurface))

        cntItemTimeSeriesLoad = np.zeros(len(isoSdSurfaceValue) + 1,
                                         dtype=np.int64)
        dataItemRFSeriesName = ["time"]
        for iTimeSeriesN in range(iTimeSeriesNS, iTimeSeriesNE):
            timeNow = eval(timeSeriesN[iTimeSeriesN])
            # Output for progress
            if ((iTimeSeriesN - iTimeSeriesNS) % 50 == 0):
                print("{0}: Working on time {1}.".format(item, timeNow))

            # Output for item
            try:
                for iItemSuffixSeries in range(len(itemSuffixSeries)):
                    fileNow = os.path.join(
                        lineSamplePath, timeSeriesN[iTimeSeriesN],
                        item + "_" + itemSuffixSeries[iItemSuffixSeries])
                    with open(fileNow) as f:
                        dataItemNameFileNow = f.readline()[0:-1].split(',')
                        dataItemFileNow = np.loadtxt(f,
                                                     delimiter=',',
                                                     skiprows=0)
                    if iItemSuffixSeries == 0:
                        dataItemNow = dataItemFileNow
                        dataItemNameNow = list(dataItemNameFileNow)
                    else:
                        dataItemNow = np.concatenate(
                            (dataItemNow, dataItemFileNow[:, 3:]), axis=1)
                        dataItemNameNow.extend(dataItemNameFileNow[3:])

                for iItemIsoSuffixSeries in range(len(itemIsoSuffixSeries)):
                    fileNow = os.path.join(
                        lineSampleIsoPath, timeSeriesN[iTimeSeriesN],
                        item + "_" + itemIsoSuffixSeries[iItemIsoSuffixSeries])
                    with open(fileNow) as f:
                        dataItemNameFileNow = f.readline()[0:-1].split(',')
                        dataItemFileNow = np.loadtxt(f,
                                                     delimiter=',',
                                                     skiprows=0)
                    dataItemNow = np.concatenate(
                        (dataItemNow, dataItemFileNow[:, 3:]), axis=1)
                    dataItemNameNow.extend(dataItemNameFileNow[3:])
            except OSError or IOError:
                print("{0}: Failed to load data at time {1}.".format(
                    item, timeNow))
                continue

            if (writeCombine):
                # Writing of combined data for each time
                fileNow = os.path.join(pathRoot, item + "-iso" + isoSdSurface,
                                       timeSeriesN[iTimeSeriesN] + ".dat")
                with open(fileNow, 'w') as f:
                    otFileHeader = 'TITLE="Item {0} at {1} s."\r\n'.format(
                        item, timeSeriesN[iTimeSeriesN])
                    otFileVariables = 'VARIABLES="{0}"\r\n'.format(
                        '" "'.join(dataItemNameNow))
                    otZoneHeader = 'ZONE T="{2} s" I={0} J={1}\r\n'.format(
                        dataItemNow.shape[0], 1, timeSeriesN[iTimeSeriesN]
                    ) + 'ZONETYPE = Ordered, DATAPACKING = BLOCK\r\n'
                    f.writelines([otFileHeader, otFileVariables, otZoneHeader])
                    for i in range(dataItemNow.shape[1]):
                        np.savetxt(f,
                                   dataItemNow[:, i],
                                   delimiter=',',
                                   newline='\n')

            # postprocess
            iIsoSdSurface = utils.first(dataItemNameNow,
                                        condition=lambda x: x == isoSdSurface)
            iQdot = utils.first(dataItemNameNow,
                                condition=lambda x: x == 'Qdot')
            QdotItemNow = dataItemNow[:, iQdot]
            indexQdotItemNow = np.where(QdotItemNow < 10)[0]

            QdotItemNow[indexQdotItemNow] = 0
            iFlameDataItemNow = scipy.signal.argrelextrema(
                QdotItemNow, np.greater)[0]
            if (len(iFlameDataItemNow) == 0):
                print("    No intense reaction zones.")
                continue
            if (leftMaxPos):
                iFlameDataItemNow = iFlameDataItemNow[0]
            else:
                iFlameDataItemNow = iFlameDataItemNow[-1]

            if cntItemTimeSeriesLoad[0] == 0:
                dataItemRFSeriesName.extend(dataItemNameNow)
                dataItemRFSeries = np.zeros(
                    (cntItemTimeSeriesLoad.shape[0],
                     iTimeSeriesNE - iTimeSeriesN, len(dataItemRFSeriesName)))

            dataItemRFSeries[0, cntItemTimeSeriesLoad[0], 0] = timeNow
            dataItemRFSeries[0, cntItemTimeSeriesLoad[0],
                             1:] = dataItemNow[iFlameDataItemNow, :]
            cntItemTimeSeriesLoad[0] += 1
            for iIsoSdSurfaceValue in range(len(isoSdSurfaceValue)):
                # Conduct interpolate for flame value
                dsItemNow = np.sqrt(
                    np.power(np.diff(dataItemNow[:, 0]), 2.0) +
                    np.power(np.diff(dataItemNow[:, 1]), 2.0))
                sItemNow = np.append(0, np.cumsum(dsItemNow))
                isoSdSurfaceLine = np.zeros_like(
                    dataItemNow[:, iIsoSdSurface]
                ) + isoSdSurfaceValue[iIsoSdSurfaceValue]
                try:
                    isoX, isoY = utils.interpolated_intercepts(
                        sItemNow, dataItemNow[:, iIsoSdSurface],
                        isoSdSurfaceLine)
                except:
                    continue
                if isoY.size == 0:
                    continue
                if (leftMaxPos):
                    sFlameItemNow = isoX[0]
                else:
                    sFlameItemNow = isoX[-1]

                if (len(indexQdotItemNow) < 2 or (not enableSavgolFilter)):
                    # Consider all interpolate
                    for isoIndex in range(1, len(dataItemRFSeriesName)):
                        dataItemRFSeries[
                            iIsoSdSurfaceValue + 1,
                            cntItemTimeSeriesLoad[iIsoSdSurfaceValue + 1],
                            isoIndex] = np.interp(sFlameItemNow, sItemNow,
                                                  dataItemNow[:, isoIndex - 1])
                else:
                    # Consider first filter and then interpolate
                    for isoIndex in range(1, len(dataItemRFSeriesName)):
                        windowLength = min(
                            savgolFilterWidth,
                            int(
                                int((indexQdotItemNow[-1] -
                                     indexQdotItemNow[0]) / 2) * 2 + 1))
                        order = 0 if windowLength < 4 else 3
                        filterData = scipy.signal.savgol_filter(
                            dataItemNow[
                                indexQdotItemNow[0]:indexQdotItemNow[-1],
                                isoIndex - 1], windowLength, order)
                        dataItemRFSeries[
                            iIsoSdSurfaceValue + 1,
                            cntItemTimeSeriesLoad[iIsoSdSurfaceValue + 1],
                            isoIndex] = np.interp(
                                sFlameItemNow, sItemNow[
                                    indexQdotItemNow[0]:indexQdotItemNow[-1]],
                                filterData)

                dataItemRFSeries[iIsoSdSurfaceValue + 1,
                                 cntItemTimeSeriesLoad[iIsoSdSurfaceValue + 1],
                                 0] = timeNow
                cntItemTimeSeriesLoad[iIsoSdSurfaceValue + 1] += 1

    # Proceeding for output
        if (dataItemRFSeriesName == ["time"]):
            print("No file is read for item {0}.".format(item))
            continue
        dataItemRFSeriesName.extend([
            "x1", "y1", "z1", "dx1/dt", "dy1/dt", "dz1/dt", \
            "x2", "y2", "z2", "dx2/dt", "dy2/dt", "dz2/dt"
        ])

        otFileItem = os.path.join(
            pathRoot, "{0}_S{1}_E{2}_iso{3}_rf{4}{5}.dat".format(
                item, timeSeriesN[iTimeSeriesNS],
                timeSeriesN[iTimeSeriesNE - 1], isoSdSurface,
                "L" if leftMaxPos else "R", "_sgF{0}".format(savgolFilterWidth)
                if enableSavgolFilter else ""))
        with open(otFileItem, "w") as f:
            otFileHeader = 'TITLE="Rf{3} for {0} from {1} to {2}."\r\n'.format(
                item, timeSeriesN[iTimeSeriesNS],
                timeSeriesN[iTimeSeriesNE - 1], "L" if leftMaxPos else "R")
            otFileVariables = 'VARIABLES="{0}"\r\n'.format(
                '" "'.join(dataItemRFSeriesName))
            f.writelines([otFileHeader, otFileVariables])

        for iIsoSdSurfaceValue in range(cntItemTimeSeriesLoad.shape[0]):
            dataItemIsoRFSeries = dataItemRFSeries[
                iIsoSdSurfaceValue,
                0:cntItemTimeSeriesLoad[iIsoSdSurfaceValue], :]
            if (dataItemIsoRFSeries.shape[0] == 0):
                continue
            dataItemIsoRFSeriesAppend1 = np.zeros(
                (dataItemIsoRFSeries.shape[0], 6))
            dataItemIsoRFSeriesAppend2 = np.zeros(
                (dataItemIsoRFSeries.shape[0], 6))
            try:
                dataItemIsoRFSeriesTXYZ = dataItemIsoRFSeries[:, 0:4]
                # x1, y1, z1 are from interpolate of time data to
                index = np.logical_and(
                    np.abs(np.diff(dataItemIsoRFSeriesTXYZ[:, 1])) < 1E-10,
                    np.abs(np.diff(dataItemIsoRFSeriesTXYZ[:, 2])) < 1E-10)
                splits = np.where(~index)[0] + 1
                # COnduct average for same XY data
                if splits.size != 0:
                    dataItemIsoRFSeriesTXYZ = np.stack(list(
                        map(lambda item: np.mean(item, axis=0),
                            np.split(dataItemIsoRFSeriesTXYZ, splits,
                                     axis=0))),
                                                       axis=0)
                dataItemIsoRFSeriesUVW = utils.moving_derivative(
                    dataItemIsoRFSeriesTXYZ[:, 1:4], 10, 3,
                    dataItemIsoRFSeriesTXYZ[:, 0])
                dataItemIsoRFSeriesAppend10 = np.hstack(
                    (dataItemIsoRFSeriesTXYZ, dataItemIsoRFSeriesUVW))

                for iAppend1 in range(1, dataItemIsoRFSeriesAppend10.shape[1]):
                    pchipAppend1 = scipy.interpolate.PchipInterpolator(
                        dataItemIsoRFSeriesAppend10[:, 0],
                        dataItemIsoRFSeriesAppend10[:, iAppend1])
                    dataItemIsoRFSeriesAppend1[:, iAppend1 - 1] = pchipAppend1(
                        dataItemIsoRFSeries[:, 0])

                # x2, y2, z2 are from interpolate of time data to
                windowLength = 15
                for i in range(3):
                    dataItemIsoRFSeriesAppend2[:,
                                               i] = scipy.signal.savgol_filter(
                                                   dataItemIsoRFSeries[:,
                                                                       i + 1],
                                                   windowLength, 3)
                dataItemIsoRFSeriesAppend2[:, 3:6] = utils.moving_derivative(
                    dataItemIsoRFSeriesAppend2[:, 0:3], windowLength, 3,
                    dataItemIsoRFSeries[:, 0])
            except:
                pass

            dataItemIsoRFSeries = np.concatenate(
                (dataItemIsoRFSeries, dataItemIsoRFSeriesAppend1), axis=1)
            dataItemIsoRFSeries = np.concatenate(
                (dataItemIsoRFSeries, dataItemIsoRFSeriesAppend2), axis=1)

            with open(otFileItem, "a") as f:
                if iIsoSdSurfaceValue == 0:
                    otZoneHeader = 'ZONE T="Qdot max" I={0} J={1}\r\n'.format(
                        dataItemIsoRFSeries.shape[0],
                        1) + 'ZONETYPE = Ordered, DATAPACKING = BLOCK\r\n'
                else:
                    otZoneHeader = 'ZONE T="{2} = {3}" I={0} J={1}\r\n'.format(
                        dataItemIsoRFSeries.shape[0], 1, isoSdSurface,
                        isoSdSurfaceValue[iIsoSdSurfaceValue - 1]
                    ) + 'ZONETYPE = Ordered, DATAPACKING = BLOCK\r\n'

                f.writelines([otZoneHeader])
                for i in range(dataItemIsoRFSeries.shape[1]):
                    np.savetxt(f,
                               dataItemIsoRFSeries[:, i],
                               delimiter=',',
                               newline='\n')


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
    parse.add_argument('--isoSdSurfaceValue', action='store', nargs='*', type=float, \
        help='Iso-surfaces values in Sd evaluation.')
    parse.add_argument('--enableSavgolFilter', action='store', nargs=1, type=utils.str2bool, \
         default=[False], help='Determine whether enable savgol filter.')
    parse.add_argument('--savgolFilterWidth', action='store', nargs=1, type=int, \
         default=[25], help='Determine the savgol filter width.')
    parse.add_argument('--leftMaxPos', action='store', nargs=1, type=utils.str2bool, \
         default=[False], help='Determine whether obtain left max position.')
    parse.add_argument('--writeTimeStepCombine', action='store', nargs=1, type=utils.str2bool, \
        default=[False], help='Determine whether write combined file for each step.')

    args = parse.parse_args()
    print(args)
    pathRoot = os.path.abspath(args.case[0])

    # Print summary of current input.
    print('******************************************************************')
    print(
        "Merge timeseries data in {0} from time {1} to time {2} for iso{3} = {4}."
        .format(pathRoot, args.timeStart, args.timeEnd, args.isoSdSurface,
                args.isoSdSurfaceValue))
    tecplotMerge(pathRoot,
                 args.timeStart[0],
                 args.timeEnd[0],
                 args.isoSdSurface[0],
                 args.isoSdSurfaceValue,
                 args.writeTimeStepCombine[0],
                 leftMaxPos=args.leftMaxPos[0],
                 enableSavgolFilter=args.enableSavgolFilter[0],
                 savgolFilterWidth=args.savgolFilterWidth[0])
    # , args.averageWidth[0])
