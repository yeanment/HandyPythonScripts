#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This is a post program for 1D flame position data (defined with maximum Qdot).
"""

import numpy as np
import sys, os
import shutil
# import argparse

def item_file_remove(item_path):
    for item in item_path:
        try:
            os.remove(item)
        except FileNotFoundError as e:
            pass
        else:
            pass

if __name__ == "__main__":
    # parse = argparse.ArgumentParser()
    # parse.add_argument('--timeStart', action=store, nargs=1, type=float)
    # args = parse.parse_args()
    
    if '--file' in sys.argv:
        file_index = sys.argv.index('--file') + 1
        dataFilePath = sys.argv[file_index]
        if (not os.access(dataFilePath, os.R_OK)):
            print('Error read dataFile.\n\r')
            exit()
        else:
            dataFileDir = os.path.dirname(dataFilePath)
            dataFileBase = os.path.basename(dataFilePath)
            dataFileName = os.path.splitext(dataFileBase)[0]
            dataFileExt = os.path.splitext(dataFileBase)[1]
    else:
        print("Please input rf file position.")
        exit()

    # Conduct data postprocessing.

    fDataFile = open(dataFilePath, "r+")
    rfTitle = fDataFile.readline()
    rfVariable = fDataFile.readline()
    varSIndex = rfVariable.find('=') 
    if (varSIndex == -1):
        exit()
    else:
        # print(type(varSIndex))
        varSList = rfVariable[varSIndex+1:-1].strip().split(" ")
        varItemList = [(item).strip()[1:-1] for item in varSList]
        itemNum = len(varItemList)
        # print(varItemList)
    rfZone = fDataFile.readline()
    zoneHeaderList = [item.strip() for item in rfZone[0:-1].split(" ")]
    try:
        iHeaderIndex = next(index for index, value in enumerate(zoneHeaderList) if (value[0].upper() == 'I'))
        jHeaderIndex = next(index for index, value in enumerate(zoneHeaderList) if (value[0].upper() == 'J'))
    except StopIteration:
        print("I/J not found.\n")
        exit()
    else:
        try:
            iNum = int(zoneHeaderList[iHeaderIndex].split("=")[1])
            jNum = int(zoneHeaderList[jHeaderIndex].split("=")[1])
        except:
            print("Not identify i/j num correctly.\r\n")
            exit()
    # read datafile
    fDataFile.close()
    data = np.loadtxt(dataFilePath, skiprows=4)
    data = np.reshape(data, [itemNum, iNum])
    # remove duplicate rf file
    yIndex = varItemList.index('Y')
    # print(yIndex)
        
    ydiffIndex = ~(np.diff(data[yIndex, :]*1E+6) < 10)
    yuIndex = np.append(True, ydiffIndex)
    ydIndex = np.append(ydiffIndex, True)
    datau = data[:, yuIndex]
    datad = data[:, ydIndex]
    dataa = (datau + datad)/2
    # average 
    print(data.shape)
    dataSplit =[np.average(itemS, axis=1) for itemS in np.split(data, np.where(ydiffIndex)[0] + 1, axis=1)]
    dataS = np.transpose(np.vstack(dataSplit))

    # save file
    tecDatarfPostPath = os.path.join(dataFileDir, "{0}_postM{1}".format(dataFileName, dataFileExt))
    with open(tecDatarfPostPath, 'w') as f:
        file_header = 'TITLE="{0}_post"\r\n'.format(dataFileName)
        variables = 'VARIABLES="{0}"\r\n'.format('" "'.join(varItemList))
        zone_header = 'ZONE T="{0}_post" I={1} J={2}\r\n  \
                    ZONETYPE = Ordered, DATAPACKING = BLOCK\r\n'.format(dataFileName, \
                    dataa.shape[1], jNum)
        f.writelines([file_header, variables, zone_header])
        for i in range(dataa.shape[0]):
            np.savetxt(f, dataa[i, :], delimiter=',', newline='\n')
        tecDatarfPostPath = os.path.join(dataFileDir, "{0}_post{1}".format(dataFileName, dataFileExt))
    
    tecDatarfPostPath = os.path.join(dataFileDir, "{0}_postA{1}".format(dataFileName, dataFileExt))
    with open(tecDatarfPostPath, 'w') as f:
        file_header = 'TITLE="{0}_post"\r\n'.format(dataFileName)
        variables = 'VARIABLES="{0}"\r\n'.format('" "'.join(varItemList))
        zone_header = 'ZONE T="{0}_post" I={1} J={2}\r\n  \
                    ZONETYPE = Ordered, DATAPACKING = BLOCK\r\n'.format(dataFileName, \
                    dataS.shape[1], jNum)
        f.writelines([file_header, variables, zone_header])
        for i in range(dataS.shape[0]):
            np.savetxt(f, dataS[i, :], delimiter=',', newline='\n')


    

