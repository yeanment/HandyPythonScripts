#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This is a converter from isoSurface Data to obtain the distribution of certain 
variables at specific time.
"""

import numpy as np
import scipy.interpolate, scipy.signal, scipy.optimize, scipy.stats
import seaborn as sns
import matplotlib.pyplot as plt
import sys, os, shutil, re, argparse
import utils

plt.style.use('/mnt/d/Work/pythonscripts/matplotlib/styles/science.mplstyle')
sns.set_style("whitegrid", rc=plt.rcParams)


def obtainFlux(pathRoot, itemName=""):
    itemListDict = utils.obtain_sorted_itemseries(pathRoot)
    itemName = list(itemListDict.keys())[0]
    timeListSuffix = list(itemListDict[itemName])
    timeListSuffix.sort(key=lambda item: eval(os.path.splitext(item)[0]))
    # Obtain the 2D data
    for iTimeSeries in range(len(timeListSuffix)):
        timeNow = timeListSuffix[iTimeSeries]
        try:
            fileNow = os.path.join(
                pathRoot, itemName + '_' + timeListSuffix[iTimeSeries])
            with open(fileNow) as f:
                dataNameFileNow = re.split(r"[\t, ]+",
                                           f.readline()[1:-1].strip())[1:]
                dataFileNow = np.loadtxt(f, delimiter='\t', skiprows=0)
                if iTimeSeries == 0:
                    dataNow = dataFileNow
                    dataNameNow = list(dataNameFileNow)
                    dataRGridLen = dataNow.shape[0]
                else:
                    dataNow = np.concatenate((dataNow, dataFileNow), axis=0)
        except OSError or IOError:
            print(
                "Failed to load data of fluxTEXT at time {0}.".format(timeNow))
            continue
    dataNameNow = [item + "Flux" for item in dataNameNow]
    dataTGrid = dataNow[:, 0].reshape((-1, dataRGridLen))[:, 0]
    dataRGrid = dataNow[:, 1].reshape((-1, dataRGridLen))[0, :]
    # interpDataNow = scipy.interpolate.interp2d(dataRGrid, dataTGrid, dataNow[:,3])
    return dataRGrid, dataTGrid, dataNow, dataNameNow


def funCToFit(x, a, b):
    return a * x + b / x


def tecplotIsoSurfaceToDistribution(pathRoot,
                                    timeStart,
                                    timeEnd,
                                    isoSdSurface,
                                    isoSdSurfaceValue,
                                    enableSavgolFilter=False,
                                    averageWidth=32,
                                    savgolFilterWidth=25,
                                    keepDataOption=2,
                                    writeCombine=False,
                                    savenpz=False):

    timeSeries = utils.obtain_sorted_timeseries(
        os.path.join(pathRoot, "sampleDictIso" + isoSdSurface))
    iTimeSeriesS = utils.first(timeSeries,
                               condition=lambda x: eval(x) >= timeStart)
    if iTimeSeriesS == None:
        print("No data for time after {0}.".format(timeStart))
        return
    if (timeEnd == -1):
        iTimeSeriesE = len(timeSeries)
    else:
        iTimeSeriesE = utils.first(timeSeries,
                                   condition=lambda x: eval(x) >= timeEnd) + 1
    if (iTimeSeriesE < iTimeSeriesS):
        print("End time is smaller than start time.")
        return

    itemPrefixDict = utils.obtain_sorted_surfaceitemseries(
        os.path.join(pathRoot, "sampleDictIso" + isoSdSurface,
                     timeSeries[iTimeSeriesS]))
    itemSeries = list(itemPrefixDict.keys())

    for iItemSeries in range(len(itemSeries)):
        itemNow = itemSeries[iItemSeries]
        # Output for item
        itemPrefixSeries = list(itemPrefixDict[itemNow])
        itemPrefixSeries.sort()

        cntItemTimeSeriesLoad = np.zeros(len(isoSdSurfaceValue) + 1,
                                         dtype=np.int32)

        if (writeCombine):
            if (not os.path.exists(os.path.join(pathRoot, itemNow))):
                os.mkdir(os.path.join(pathRoot, itemNow))
        if (not os.path.exists(os.path.join(pathRoot, "PDF_" + itemNow))):
            os.mkdir(os.path.join(pathRoot, "PDF_" + itemNow))
        if (not os.path.exists(os.path.join(pathRoot, "PDFPLOT_" + itemNow))):
            os.mkdir(os.path.join(pathRoot, "PDFPLOT_" + itemNow))

        for iTimeSeriesN in range(iTimeSeriesS, iTimeSeriesE):
            timeNow = eval(timeSeries[iTimeSeriesN])
            # Output for progress
            if ((iTimeSeriesN - iTimeSeriesS) % 50 == 0):
                print("Working on time {0}.".format(timeSeries[iTimeSeriesN]))

            # trying reading the datafile
            try:
                for iItemPrefixSeries in range(len(itemPrefixSeries)):
                    fileNow = os.path.join(
                        pathRoot, "sampleDictIso" + isoSdSurface,
                        timeSeries[iTimeSeriesN],
                        itemPrefixSeries[iItemPrefixSeries] + '_' + itemNow +
                        '.raw')
                    with open(fileNow) as f:
                        dataItemNameFileNow = f.readline()[0:-1].split(',')
                        dataItemNameFileNow = re.split(
                            r"[ ]+",
                            f.readline()[1:-1].strip())
                        dataItemFileNow = np.loadtxt(f,
                                                     delimiter=' ',
                                                     skiprows=0)
                    if iItemPrefixSeries == 0:
                        dataItemNow = dataItemFileNow
                        dataItemNameNow = list(dataItemNameFileNow)
                    else:
                        dataItemNow = np.concatenate(
                            (dataItemNow, dataItemFileNow[:, 3:]), axis=1)
                        dataItemNameNow.extend(dataItemNameFileNow[3:])
            except OSError or IOError:
                print("Failed to load data at time {0}.".format(timeNow))
                continue

            # Select to keep which parts of data
            # x, y the same, average z
            if (keepDataOption == 1 or keepDataOption == 2):
                index = np.logical_and(
                    np.abs(np.diff(dataItemNow[:, 0])) < 1E-20,
                    np.abs(np.diff(dataItemNow[:, 1])) < 1E-20)
                splits = np.where(~index)[0] + 1
                if splits.size != 0:
                    dataItemNow = np.stack(list(
                        map(lambda item:np.mean(item, axis=0), \
                            filter(lambda item:item.shape[0] == keepDataOption, \
                                   np.split(dataItemNow, splits, axis=0)))), axis=0)
            # sortIndex = np.argsort(dataItemNow[:, 1])
            # dataItemNow = dataItemNow[sortIndex, :]
            # determine the local tangential vector
            tck, u = scipy.interpolate.splprep(dataItemNow[:, 0:3].transpose(),
                                               s=0)
            dx, dy, dz = scipy.interpolate.splev(u, tck, der=1)
            dz *= 0
            ds = np.sqrt(dx * dx + dy * dy + dz * dz)
            dataItemNow = np.concatenate(
                (dataItemNow, (dx / ds)[:, np.newaxis]), axis=1)
            dataItemNow = np.concatenate(
                (dataItemNow, (dy / ds)[:, np.newaxis]), axis=1)
            dataItemNow = np.concatenate(
                (dataItemNow, (dz / ds)[:, np.newaxis]), axis=1)
            dataItemNameNow.extend(['t_x', 't_y', 't_z'])

            if (writeCombine):
                # Writing of combined data for each time
                fileNow = os.path.join(pathRoot, itemNow,
                                       timeSeries[iTimeSeriesN] + ".dat")
                with open(fileNow, 'w') as f:
                    otFileHeader = 'TITLE="{0} = {1} at {2} s."\r\n'.format(
                        itemNow[0:len("iso" + isoSdSurface)],
                        itemNow[len("iso" + isoSdSurface):],
                        timeSeries[iTimeSeriesN])
                    otFileVariables = 'VARIABLES="{0}"\r\n'.format(
                        '" "'.join(dataItemNameNow))
                    otZoneHeader = 'ZONE T="{2} s" I={0} J={1}\r\n'.format(
                        dataItemNow.shape[0], 1, timeSeries[iTimeSeriesN]
                    ) + 'ZONETYPE = Ordered, DATAPACKING = BLOCK\r\n'
                    f.writelines([otFileHeader, otFileVariables, otZoneHeader])
                    for i in range(dataItemNow.shape[1]):
                        np.savetxt(f,
                                   dataItemNow[:, i],
                                   delimiter=',',
                                   newline='\n')

            # Postprocess for the frequence
            if (not os.path.exists(
                    os.path.join(pathRoot, "PDFPLOT_" + itemNow,
                                 timeSeries[iTimeSeriesN]))):
                os.mkdir(
                    os.path.join(pathRoot, "PDFPLOT_" + itemNow,
                                 timeSeries[iTimeSeriesN]))

            # Writing of combined data for each time
            fileNow = os.path.join(pathRoot, "PDF_" + itemNow,
                                   timeSeries[iTimeSeriesN] + ".dat")
            with open(fileNow, 'w') as f:
                otFileHeader = 'TITLE="PDF of {0} = {1} at {2} s."\r\n'.format(
                    itemNow[0:len("iso" + isoSdSurface)],
                    itemNow[len("iso" + isoSdSurface):],
                    timeSeries[iTimeSeriesN])
                otFileVariables = 'VARIABLES="Var" "Density"\r\n'.format(
                    '" "'.join(dataItemNameNow))
                f.writelines([otFileHeader, otFileVariables])
                for i in range(dataItemNow.shape[1]):
                    if (re.match(".*_sd" + isoSdSurface + "$",
                                 dataItemNameNow[i]) == None):
                        continue
                    otZoneHeader = 'ZONE T="{2}" I={0} J={1}\r\n'.format(
                        500, 1, dataItemNameNow[i]
                    ) + 'ZONETYPE = Ordered, DATAPACKING = BLOCK\r\n'
                    f.writelines([otZoneHeader])
                    density = scipy.stats.gaussian_kde(dataItemNow[:, i])
                    dist = np.max(dataItemNow[:, i]) - np.min(dataItemNow[:,
                                                                          i])
                    var = np.linspace(
                        np.min(dataItemNow[:, i]) - dist * 5,
                        np.max(dataItemNow[:, i]) + dist * 5, 500)
                    np.savetxt(f, var, delimiter=',', newline='\n')
                    np.savetxt(f, density(var), delimiter=',', newline='\n')

                    ## Plot using the seaborn
                    pdfPlotFileNow = os.path.join(pathRoot,
                                                  "PDFPLOT_" + itemNow,
                                                  timeSeries[iTimeSeriesN],
                                                  dataItemNameNow[i] + ".png")
                    g = plt.figure()
                    ax = sns.kdeplot(dataItemNow[:, i])
                    ax = sns.rugplot(dataItemNow[:, i], height=.2)
                    ax.set_xlabel(dataItemNameNow[i])
                    ax.set_ylabel("Density")
                    ax.set_title("{0} s, mean({1}) = {2}".format(
                        timeSeries[iTimeSeriesN], dataItemNameNow[i],
                        np.mean(dataItemNow[:, i])))
                    # g = sns.displot(dataItemNow[:, i], kind="kde", rug=True, height=4, aspect=1.25)
                    # g.despine=False
                    # g.set_axis_labels(dataItemNameNow[i], "Density")
                    # g.set_titles(col_template=timeSeries[iTimeSeriesN]+" s")
                    g.tight_layout()
                    plt.show()
                    g.savefig(pdfPlotFileNow)
                    plt.close(g)


if __name__ == "__main__":
    parse = argparse.ArgumentParser( \
        description='Process some input parameters for post-processing.')
    parse.add_argument('--case', action='store', nargs=1, type=utils.check_path, \
        default = '.', help='Root path.')
    parse.add_argument('--timeStart', action='store', nargs=1, type=float, \
        default = [0.0], help='Starting time.')
    parse.add_argument('--timeEnd', action='store', nargs=1, type=float, \
        default = [-1.0], help='Ending time.')

    parse.add_argument('--isoSdSurface', action='store', nargs=1, type=str, \
        help='Iso-surfaces for Sd evaluation.')
    parse.add_argument('--isoSdSurfaceValue', action='store', nargs='*', type=float, \
        help='Iso-surfaces values in Sd evaluation.')

    parse.add_argument('--enableSavgolFilter', action='store', nargs=1, type=utils.str2bool, \
         default=[False], help='Determine whether enable savgol filter.')
    parse.add_argument('--averageWidth', action='store', nargs=1, type=int, \
        default=[32], help='Determine whether write combined file for each step.')
    parse.add_argument('--savgolFilterWidth', action='store', nargs=1, type=int, \
        default=[25], help='Determine whether write combined file for each step.')

    parse.add_argument('--keepDataOption', action='store', nargs=1, type=int, \
        default=[2], help='Determine kept data point. 1 for z = 0; 2 for z = +-; 0 for all.')

    parse.add_argument('--writeTimeStepCombine', action='store', nargs=1, type=utils.str2bool, \
        default=[False], help='Determine whether write combined file for each step.')
    parse.add_argument('--savenpz', action='store', nargs=1, type=utils.str2bool, \
        default=[False], help='Determine whether save combined npz file.')

    args = parse.parse_args()
    print(args)

    # Print summary of current input.
    print('*******************************')
    pathRoot = os.path.abspath(args.case[0])

    print(
        "Transform timeseries data in {0} from time {1} to time {2} for iso {3} = {4}."
        .format(pathRoot, args.timeStart, args.timeEnd, args.isoSdSurface,
                args.isoSdSurfaceValue))
    tecplotIsoSurfaceToDistribution(
        pathRoot, args.timeStart[0], args.timeEnd[0], args.isoSdSurface[0],
        args.isoSdSurfaceValue, args.enableSavgolFilter[0],
        args.averageWidth[0], args.savgolFilterWidth[0],
        args.keepDataOption[0], args.writeTimeStepCombine[0], args.savenpz[0])
