#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This is a converter from 1D line timeseries data to Tecplot 2D data format and
get flame position (defined with maximum Qdot).
"""

import numpy as np
import scipy.interpolate, scipy.signal, scipy.optimize
import sys, os, shutil, re, argparse
import utils


def obtainFlux(pathRoot, itemName = ""):
    itemListDict = utils.obtain_sorted_itemseries(pathRoot)
    itemName = list(itemListDict.keys())[0]
    timeListSuffix = list(itemListDict[itemName])
    timeListSuffix.sort(key=lambda item:eval(os.path.splitext(item)[0]))
    # Obtain the 2D data
    for iTimeSeries in range(len(timeListSuffix)):
        timeNow = timeListSuffix[iTimeSeries]
        try:
            fileNow = os.path.join(pathRoot, itemName + '_' + timeListSuffix[iTimeSeries])
            with open(fileNow) as f:
                dataNameFileNow = re.split(r"[\t, ]+", f.readline()[1:-1].strip())[1:]
                dataFileNow = np.loadtxt(f, delimiter='\t', skiprows=0)
                if iTimeSeries == 0:
                    dataNow = dataFileNow
                    dataNameNow = list(dataNameFileNow)
                    dataRGridLen = dataNow.shape[0]
                else:
                    dataNow = np.concatenate((dataNow, dataFileNow), axis=0)
        except OSError or IOError:
            print("Failed to load data of fluxTEXT at time {0}.".format(timeNow))
            continue
    dataNameNow = [item + "Flux" for item in dataNameNow]
    dataTGrid = dataNow[:, 0].reshape((-1, dataRGridLen))[:, 0]
    dataRGrid = dataNow[:, 1].reshape((-1, dataRGridLen))[0,:]
    # interpDataNow = scipy.interpolate.interp2d(dataRGrid, dataTGrid, dataNow[:,3])
    return dataRGrid, dataTGrid, dataNow, dataNameNow


def funCToFit(x, a, b):
    return a*x + b/x


def tecplotMerge(pathRoot, timeStart, timeEnd, isoSdSurface, isoSdSurfaceValue, \
    savenpz=False, writeCombine=False, averageWdith=32, keepDataOption=0, savgolFilterWidth=25):
    
    timeSeries = utils.obtain_sorted_timeseries(pathRoot)
    iTimeSeriesS = utils.first(timeSeries, condition=lambda x:eval(x)>=timeStart)
    if iTimeSeriesS == None:
        print("No data for time after {0}.".format(timeStart))
        return
    if (timeEnd == -1):
        iTimeSeriesE = len(timeSeries)
    else:
        iTimeSeriesE = utils.first(timeSeries, condition=lambda x:eval(x)>=timeEnd) + 1
    if (iTimeSeriesE < iTimeSeriesS):
        print("End time is smaller than start time.")
        return
    itemPrefixDict = utils.obtain_sorted_surfaceitemseries(os.path.join(pathRoot, timeSeries[iTimeSeriesS]))
    itemSeries = list(itemPrefixDict.keys())

    isoSdSurfaceValue.sort()
    lenIsoSdSurfaceValue = len(isoSdSurfaceValue)
    cntTimeSeriesLoad = np.zeros((len(itemSeries), lenIsoSdSurfaceValue + 1), dtype=np.int32)

    dataRFSeries = dict()
    dataRFSeriesName = dict()

    for iTimeSeriesN in range(iTimeSeriesS, iTimeSeriesE):
        timeNow = eval(timeSeries[iTimeSeriesN])
        # Output for progress
        if((iTimeSeriesN - iTimeSeriesS)%50 == 0):
            print("Working on time {0}.".format(timeSeries[iTimeSeriesN]))
        for iItemSeries in range(len(itemSeries)):
            itemNow = itemSeries[iItemSeries]
            # Output for item
            itemPrefixSeries = list(itemPrefixDict[itemNow])
            itemPrefixSeries.sort()
            # trying reading the datafile
            try:
                for iItemPrefixSeries in range(len(itemPrefixSeries)):
                    fileNow = os.path.join(pathRoot, timeSeries[iTimeSeriesN], \
                        itemPrefixSeries[iItemPrefixSeries] + '_' + itemNow + '.raw')
                    with open(fileNow) as f:
                        dataItemNameFileNow = f.readline()[0:-1].split(',')
                        dataItemNameFileNow = re.split(r"[ ]+", f.readline()[1:-1].strip())
                        dataItemFileNow = np.loadtxt(f, delimiter=' ', skiprows=0)
                    if iItemPrefixSeries == 0:
                        dataItemNow = dataItemFileNow
                        dataItemNameNow = list(dataItemNameFileNow)
                    else:
                        dataItemNow = np.concatenate((dataItemNow, dataItemFileNow[:, 3:]), axis=1)
                        dataItemNameNow.extend(dataItemNameFileNow[3:])
            except OSError or IOError:
                print("Failed to load data at time {0}.".format(timeNow))
                continue
            
            # postProcess for imported data
            # x, y the same, average z
            if(keepDataOption):
                index = np.logical_and(np.abs(np.diff(dataItemNow[:, 0])) < 1E-20, np.abs(np.diff(dataItemNow[:, 1])) < 1E-20)
                splits = np.where(~index)[0] + 1
                if splits.size != 0:    
                    dataItemNow = np.stack(list(
                        map(lambda item:np.mean(item, axis=0), \
                            filter(lambda item:item.shape[0] == keepDataOption, \
                                   np.split(dataItemNow, splits, axis=0)))), axis=0)
            sortIndex = np.argsort(dataItemNow[:, 1])
            dataItemNow = dataItemNow[sortIndex, :]
            # determine the local tangential vector
            tck, u = scipy.interpolate.splprep(dataItemNow[:, 0:3].transpose(), s=0)
            dx, dy, dz = scipy.interpolate.splev(u, tck,der=1)
            dz *= 0
            ds = np.sqrt(dx*dx + dy*dy + dz*dz)
            dataItemNow = np.concatenate((dataItemNow, (dx/ds)[:,np.newaxis]), axis=1)
            dataItemNow = np.concatenate((dataItemNow, (dy/ds)[:,np.newaxis]), axis=1)
            dataItemNow = np.concatenate((dataItemNow, (dz/ds)[:,np.newaxis]), axis=1)
            dataItemNameNow.extend(['n_x', 'n_y', 'n_z'])


            # Obtain the u_m from interpolated value
            # Use data from 0.03 to 0.04
            yNow = dataItemNow[:, 1]
            UyNow = dataItemNow[:, utils.first(dataItemNameNow, condition=lambda x:x=='U_y')]
            iToInterop = np.where( (yNow > 0.03) & (yNow < 0.035))[0]
            popt, pcov = scipy.optimize.curve_fit(funCToFit, yNow[iToInterop], UyNow[iToInterop], method='trf')
            A = np.vstack([yNow[iToInterop], np.ones(len(yNow[iToInterop]))]).T
            m, c = np.linalg.lstsq(A, UyNow[iToInterop], rcond=None)[0]

            if(writeCombine):
                # write to file 
                tecDataPathN = os.path.join(pathRoot, timeSeries[iTimeSeriesN], \
                                            'tecComb-Time_' + timeSeries[iTimeSeriesN] + \
                                            '-' + itemNow + '.dat')
                with open(tecDataPathN, 'w') as f:
                    file_header = 'TITLE="{0}-{1}"\r\n'.format(timeNow, itemNow)
                    variables = 'VARIABLES="{0}" "U_m (kr + b)" "U_m (Ar + B/r)" \r\n'.format('" "'.join(dataItemNameNow))
                    zone_header = 'ZONE T="{0}-{1}" I={2} J={3}\r\n \
                                ZONETYPE = Ordered, DATAPACKING = BLOCK\r\n'.format( \
                                timeNow, itemNow, dataItemNow.shape[0], 1)
                    f.writelines([file_header, variables, zone_header])
                with open(tecDataPathN, 'ab') as f:
                    for i in range(dataItemNow.shape[1]):
                        np.savetxt(f, dataItemNow[:, i], delimiter=',', newline='\n')
                    np.savetxt(f, m*dataItemNow[:, 1] + c, delimiter=',', newline='\n')
                    np.savetxt(f, funCToFit(dataItemNow[:, 1] + 1E-30, *popt), delimiter=',', newline='\n')
            
            # postProcess for ZBilger,
            # generate the data item rf file
            if cntTimeSeriesLoad[iItemSeries, 0] == 0:
                dataItemRFSeriesName = ["time"]
                dataItemRFSeriesName.extend(dataItemNameNow)
                # Obtain the interpolated local velocity
                dataItemRFSeriesName.extend(["u_m"])
                dataItemRFSeries = np.zeros((cntTimeSeriesLoad.shape[1], iTimeSeriesE - iTimeSeriesS, len(dataItemRFSeriesName)))
                dataRFSeries.update({itemNow:dataItemRFSeries})
                dataRFSeriesName.update({itemNow:dataItemRFSeriesName})
            dataItemRFSeries = dataRFSeries[itemNow]
            dataItemRFSeriesName = dataRFSeriesName[itemNow]

            # First consider iso surface? or maximal curvature/Qdot
            iIsoSdSurface = utils.first(dataItemNameNow, condition=lambda x:x==isoSdSurface)
            iQdot = utils.first(dataItemNameNow, condition=lambda x:x=='Qdot')

            # matplotlib to see the difference


            for iIsoSdSurfaceValue in range(lenIsoSdSurfaceValue):
                # Conduct interpolate for flame value
                dsItemNow = np.sqrt(np.power(np.diff(dataItemNow[:, 0]), 2.0) + np.power(np.diff(dataItemNow[:, 1]), 2.0))
                sItemNow = np.append(0, np.cumsum(dsItemNow))
                isoSdSurfaceLine = np.zeros_like(dataItemNow[:, iIsoSdSurface]) + isoSdSurfaceValue[iIsoSdSurfaceValue]
                try:
                    isoX, isoY = utils.interpolated_intercepts(sItemNow, dataItemNow[:, iIsoSdSurface], isoSdSurfaceLine)
                except:
                    continue
                if isoY.size == 0:
                    continue
                sFlameItemNow = isoX[-1]
                QdotItemNow = dataItemNow[:, iQdot]
                indexQdotItemNow = np.where(QdotItemNow > 0.01)[0]
                # if(len(indexQdotItemNow) < 2):
                #     # Consider all interpolate
                for iIsoItemDataIndex in range(1, len(dataRFSeriesName[itemNow]) - 1):
                    dataItemRFSeries[iIsoSdSurfaceValue, cntTimeSeriesLoad[iItemSeries, iIsoSdSurfaceValue], iIsoItemDataIndex] \
                                          = np.interp(sFlameItemNow, sItemNow, dataItemNow[:, iIsoItemDataIndex - 1])
                # else:
                #     # Consider first filter and then interpolate
                #     for iIsoItemDataIndex in range(1, len(dataRFSeriesName[itemNow])):
                #         if (indexQdotItemNow[-1] - indexQdotItemNow[0])%2 == 1:
                #             windowLength = min(savgolFilterWidth, indexQdotItemNow[-1] - indexQdotItemNow[0])
                #         else:
                #             windowLength = min(savgolFilterWidth, indexQdotItemNow[-1] - indexQdotItemNow[0] - 1)
                #         if windowLength < 4:
                #             order = 0
                #         else:
                #             order = 3
                #         filterData = scipy.signal.savgol_filter(dataItemNow[indexQdotItemNow[0]:indexQdotItemNow[-1], iIsoItemDataIndex - 1], windowLength, order)
                #         dataItemRFSeries[iIsoSdSurfaceValue, cntTimeSeriesLoad[iItemSeries, iIsoSdSurfaceValue], iIsoItemDataIndex] \
                #                          = np.interp(sFlameItemNow, sItemNow[indexQdotItemNow[0]:indexQdotItemNow[-1]], filterData)                    
                
                
                # dataItemRFSeries[iIsoSdSurfaceValue, cntTimeSeriesLoad[iItemSeries, iIsoSdSurfaceValue], len(dataRFSeriesName[itemNow]) - 1] \
                #         = m*dataItemRFSeries[iIsoSdSurfaceValue, cntTimeSeriesLoad[iItemSeries, iIsoSdSurfaceValue], 2] + c
                dataItemRFSeries[iIsoSdSurfaceValue, cntTimeSeriesLoad[iItemSeries, iIsoSdSurfaceValue], len(dataRFSeriesName[itemNow]) - 1] \
                        = funCToFit(dataItemRFSeries[iIsoSdSurfaceValue, cntTimeSeriesLoad[iItemSeries, iIsoSdSurfaceValue], 2] + 1E-30, *popt)
                        

                dataItemRFSeries[iIsoSdSurfaceValue, cntTimeSeriesLoad[iItemSeries, iIsoSdSurfaceValue], 0] = timeNow
                cntTimeSeriesLoad[iItemSeries, iIsoSdSurfaceValue] += 1
            # Here for Qdot or CurlMax
            QdotNm = np.max(dataItemNow[:, iQdot]) # np.argmax
            iCurl = utils.first(dataItemNameNow, condition=lambda x:x==('cur_sd' + isoSdSurface))
            index_QdotmList = np.where(\
                np.logical_and(np.diff(np.diff(dataItemNow[:, iCurl])) < 0, dataItemNow[1:-1, iQdot] > 0.1*QdotNm))[0]
            if index_QdotmList.size == 0:
                continue
            iFlameDataItemNow = index_QdotmList[np.argmax(dataItemNow[index_QdotmList, iCurl])]
            dataItemRFSeries[lenIsoSdSurfaceValue, cntTimeSeriesLoad[iItemSeries, lenIsoSdSurfaceValue], 0] = timeNow
            dataItemRFSeries[lenIsoSdSurfaceValue, cntTimeSeriesLoad[iItemSeries, lenIsoSdSurfaceValue], 1:len(dataRFSeriesName[itemNow]) - 1] = dataItemNow[iFlameDataItemNow, :]

            # dataItemRFSeries[lenIsoSdSurfaceValue, cntTimeSeriesLoad[iItemSeries, lenIsoSdSurfaceValue], len(dataRFSeriesName[itemNow]) - 1] \
            #         = m*dataItemRFSeries[lenIsoSdSurfaceValue, cntTimeSeriesLoad[iItemSeries, lenIsoSdSurfaceValue], 2] + c
            dataItemRFSeries[lenIsoSdSurfaceValue, cntTimeSeriesLoad[iItemSeries, lenIsoSdSurfaceValue], len(dataRFSeriesName[itemNow]) - 1] \
                    = funCToFit(dataItemRFSeries[lenIsoSdSurfaceValue, cntTimeSeriesLoad[iItemSeries, lenIsoSdSurfaceValue], 2] + 1E-30, *popt)
            print("fit at time {0}: {1}".format(timeNow, popt))

            cntTimeSeriesLoad[iItemSeries, lenIsoSdSurfaceValue] += 1
            
            # # determine the interpolated local velocity
            # indexColdS = utils.first(dataN[:, 1], condition=lambda x:x>0.03)
            # indexColdE = utils.first(dataN[:, 1], condition=lambda x:x>0.04)
            # A = np.vstack([dataN[indexColdS:indexColdE, 1], np.ones(indexColdE - indexColdS)]).T
            # m, c = np.linalg.lstsq(A, \
            #     dataN[indexColdS:indexColdE, utils.first(item_data, condition=lambda x:x=='U_y')])[0]
            # dataf[0, iter_e[0], len(item_data)+1] = m*dataf[0, iter_e[0], 2] + c
    
    dataRGrid, dataTGrid, dataRTFlux, dataRTFluxName = obtainFlux(os.path.join(pathRoot, "../fluxTest", timeSeries[-1]))
    areaInterp = scipy.interpolate.interp2d(dataRGrid, dataTGrid, dataRTFlux[:,2].reshape((len(dataTGrid), -1)), kind='cubic')
    rhoInterp = scipy.interpolate.interp2d(dataRGrid, dataTGrid, dataRTFlux[:,3].reshape((len(dataTGrid), -1)), kind='cubic')
    rhoUInterp = scipy.interpolate.interp2d(dataRGrid, dataTGrid, dataRTFlux[:, 4].reshape((len(dataTGrid), -1)), kind='cubic')


    # Proceding for output
    for iItemSeries in range(len(itemSeries)):
        itemNow = itemSeries[iItemSeries]
        dataItemRFSeries = dataRFSeries[itemNow]
        dataItemRFSeriesName = dataRFSeriesName[itemNow]

        dataItemRFSeriesName.extend(["x1", "y1", "z1", "dx1/dt", "dy1/dt", "dz1/dt", \
                                     "x2", "y2", "z2", "dx2/dt", "dy2/dt", "dz2/dt"])
        dataItemRFSeriesName.extend(dataRTFluxName[2:])
        
        otFileItem = os.path.join(pathRoot, "{0}_S{1}_E{2}_iso{3}_rf_sgF{4}_avWid{5}.dat".format(\
            itemNow, timeSeries[iTimeSeriesS], timeSeries[iTimeSeriesE-1], isoSdSurface, savgolFilterWidth, averageWdith))
        with open(otFileItem, "w") as f:
            otFileHeader = 'TITLE="Rf for {0} from {1} to {2}."\r\n'.format(itemNow, timeSeries[iTimeSeriesS], timeSeries[iTimeSeriesE-1])
            otFileVariables = 'VARIABLES="{0}"\r\n'.format('" "'.join(dataItemRFSeriesName))
            f.writelines([otFileHeader, otFileVariables])
        
        for iIsoSdSurfaceValue in range(cntTimeSeriesLoad.shape[1]):
            dataItemIsoRFSeries = dataItemRFSeries[iIsoSdSurfaceValue, 0:cntTimeSeriesLoad[iItemSeries, iIsoSdSurfaceValue], :]
            if(dataItemIsoRFSeries.shape[0] == 0):
                continue
            dataItemIsoRFSeriesAppend1 = np.zeros((dataItemIsoRFSeries.shape[0], 6))
            dataItemIsoRFSeriesAppend2 = np.zeros((dataItemIsoRFSeries.shape[0], 6))
            try:
                dataItemIsoRFSeriesTXYZ = dataItemIsoRFSeries[:, 0:4]
                # x1, y1, z1 are from interpolate of time data to 
                index = np.logical_and(np.abs(np.diff(dataItemIsoRFSeriesTXYZ[:, 1])) < 1E-20, \
                                       np.abs(np.diff(dataItemIsoRFSeriesTXYZ[:, 2])) < 1E-20)
                splits = np.where(~index)[0] + 1
                # Conduct average for same XY data
                if splits.size != 0:
                    dataItemIsoRFSeriesTXYZ = np.stack(list(\
                        map(lambda item:np.mean(item, axis=0), np.split(dataItemIsoRFSeriesTXYZ, splits, axis=0))), axis=0)
                dataItemIsoRFSeriesUVW = utils.moving_derivative(dataItemIsoRFSeriesTXYZ[:,1:4], averageWdith, 2, dataItemIsoRFSeriesTXYZ[:, 0])
                dataItemIsoRFSeriesAppend10 = np.hstack((dataItemIsoRFSeriesTXYZ, dataItemIsoRFSeriesUVW))
                
                for iAppend1 in range(1, dataItemIsoRFSeriesAppend10.shape[1]):
                    pchipAppend1 = scipy.interpolate.PchipInterpolator(dataItemIsoRFSeriesAppend10[:, 0], dataItemIsoRFSeriesAppend10[:, iAppend1])
                    dataItemIsoRFSeriesAppend1[:, iAppend1 - 1] = pchipAppend1(dataItemIsoRFSeries[:, 0])            
            except:
                pass
            
            try:
                # x2, y2, z2 are from interpolate of time data to 
                for i in range(3):
                    dataItemIsoRFSeriesAppend2[:, i] = scipy.signal.savgol_filter(dataItemIsoRFSeries[:, i+1], savgolFilterWidth, 3)
                    # dataItemIsoRFSeriesAppend2[:, i+3] = scipy.signal.savgol_filter(dataItemIsoRFSeries[:, i+1], savgolFilterWidth, 3, deriv=1)
                dataItemIsoRFSeriesAppend2[:, 3:6] = utils.moving_derivative(dataItemIsoRFSeries[:, 1:4], averageWdith, 2, dataItemIsoRFSeries[:, 0])
            except:
                pass

            dataItemIsoRFSeries = np.concatenate((dataItemIsoRFSeries, dataItemIsoRFSeriesAppend1), axis = 1)
            dataItemIsoRFSeries = np.concatenate((dataItemIsoRFSeries, dataItemIsoRFSeriesAppend2), axis = 1)

            ## Interpolate for RTFlux
            dataItemIsoRFSeriesAppendRT = np.zeros((dataItemIsoRFSeries.shape[0], 3))
            dataItemIsoRFSeriesAppendRT[:, 0] = \
                np.squeeze([areaInterp(*p) for p in zip(dataItemIsoRFSeries[:, 2], dataItemIsoRFSeries[:, 0])])
            dataItemIsoRFSeriesAppendRT[:, 1] = \
                np.squeeze([rhoInterp(*p) for p in zip(dataItemIsoRFSeries[:, 2], dataItemIsoRFSeries[:, 0])])
            dataItemIsoRFSeriesAppendRT[:, 2] = \
                np.squeeze([rhoUInterp(*p) for p in zip(dataItemIsoRFSeries[:, 2], dataItemIsoRFSeries[:, 0])])
            
            dataItemIsoRFSeries = np.concatenate((dataItemIsoRFSeries, dataItemIsoRFSeriesAppendRT), axis = 1)

            with open(otFileItem, "a") as f:
                if iIsoSdSurfaceValue == lenIsoSdSurfaceValue:
                    otZoneHeader = 'ZONE T="Curl max" I={0} J={1}\r\nZONETYPE = Ordered, DATAPACKING = BLOCK\r\n'.format( \
                                    dataItemIsoRFSeries.shape[0], 1)
                else:
                    otZoneHeader = 'ZONE T="{2} = {3}" I={0} J={1}\r\n  \
                                    ZONETYPE = Ordered, DATAPACKING = BLOCK\r\n'.format( \
                                    dataItemIsoRFSeries.shape[0], 1, isoSdSurface, isoSdSurfaceValue[iIsoSdSurfaceValue])
                f.writelines([otZoneHeader])
                for i in range(dataItemIsoRFSeries.shape[1]):
                    np.savetxt(f, dataItemIsoRFSeries[:, i], delimiter=',', newline='\n')

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
    
    parse.add_argument('--averageWidth', action='store', nargs=1, type=int, \
        default=[32], help='Determine whether write combined file for each step.')

    parse.add_argument('--savgolFilterWidth', action='store', nargs=1, type=int, \
        default=[25], help='Determine whether write combined file for each step.')

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
    pathRoot = os.path.abspath(args.case[0])
    print("Merge timeseries data in {0} from time {1} to time {2} for iso {3} = {4}.".format(
        pathRoot, args.timeStart, args.timeEnd, args.isoSdSurface, args.isoSdSurfaceValue))
    tecplotMerge(pathRoot, args.timeStart[0], args.timeEnd[0], \
        args.isoSdSurface[0], args.isoSdSurfaceValue, \
        args.savenpz[0], args.writeTimeStepCombine[0], \
        args.averageWidth[0], args.keepDataOption[0], args.savgolFilterWidth[0])

    
