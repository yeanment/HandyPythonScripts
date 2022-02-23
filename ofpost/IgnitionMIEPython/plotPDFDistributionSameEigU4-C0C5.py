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
sys.path.append(os.getcwd() + "/..")
sys.path.append(os.getcwd() + "/.")

import utils


def obtainCombinedData(pathRoot,
                       timeNow,
                       isoSdSurface,
                       isoSdSurfaceValue,
                       keepDataOption=2):
    itemNow = "iso" + isoSdSurface + isoSdSurfaceValue
    pathRootTimeNow = os.path.join(pathRoot, "sampleDictIso" + isoSdSurface,
                                   timeNow)
    itemPrefixDict = utils.obtain_sorted_surfaceitemseries(pathRootTimeNow)
    itemPrefixSeries = list(itemPrefixDict[itemNow])
    itemPrefixSeries.sort()

    for iItemPrefixSeries in range(len(itemPrefixSeries)):
        fileNow = os.path.join(
            pathRootTimeNow,
            itemPrefixSeries[iItemPrefixSeries] + '_' + itemNow + '.raw')
        with open(fileNow) as f:
            dataItemNameFileNow = f.readline()[0:-1].split(',')
            dataItemNameFileNow = re.split(r"[ ]+", f.readline()[1:-1].strip())
            dataItemFileNow = np.loadtxt(f, delimiter=' ', skiprows=0)
        if iItemPrefixSeries == 0:
            dataItemNow = dataItemFileNow
            dataItemNameNow = list(dataItemNameFileNow)
        else:
            dataItemNow = np.concatenate((dataItemNow, dataItemFileNow[:, 3:]),
                                         axis=1)
            dataItemNameNow.extend(dataItemNameFileNow[3:])

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
    sortIndex = np.argsort(dataItemNow[:, 0])
    dataItemNow = dataItemNow[sortIndex, :]
    # determine the local tangential vector
    tck, u = scipy.interpolate.splprep(dataItemNow[:, 0:3].transpose(), s=0)
    dx, dy, dz = scipy.interpolate.splev(u, tck, der=1)
    dz *= 0
    ds = np.sqrt(dx * dx + dy * dy + dz * dz)
    dataItemNow = np.concatenate((dataItemNow, (dx / ds)[:, np.newaxis]),
                                 axis=1)
    dataItemNow = np.concatenate((dataItemNow, (dy / ds)[:, np.newaxis]),
                                 axis=1)
    dataItemNow = np.concatenate((dataItemNow, (dz / ds)[:, np.newaxis]),
                                 axis=1)

    dataItemNameNow.extend(['t_x', 't_y', 't_z'])

    iCur = utils.first(dataItemNameNow,
                       condition=lambda x: x == "cur_sd" + isoSdSurface)
    iIgradU = utils.first(dataItemNameNow,
                          condition=lambda x: x ==
                          ("igradU_sd" + isoSdSurface))
    iNNgradU = utils.first(dataItemNameNow,
                           condition=lambda x: x ==
                           ("nngradU_sd" + isoSdSurface))
    curnormalized1 = dataItemNow[:, iCur] * 1.1344992
    curnormalized2 = dataItemNow[:, iCur] * 0.5431198 / 1000
    dataItemNow = np.concatenate((dataItemNow, curnormalized1[:, np.newaxis]),
                                 axis=1)
    dataItemNow = np.concatenate((dataItemNow, curnormalized2[:, np.newaxis]),
                                 axis=1)

    Kas = (dataItemNow[:, iIgradU] -
           dataItemNow[:, iNNgradU]) * 0.5431198 / 1000 / 1.1344992
    iN0 = utils.first(dataItemNameNow,
                      condition=lambda x: x == ("n_sd" + isoSdSurface + "_x"))
    iN1 = utils.first(dataItemNameNow,
                      condition=lambda x: x == ("n_sd" + isoSdSurface + "_y"))
    iN2 = utils.first(dataItemNameNow,
                      condition=lambda x: x == ("n_sd" + isoSdSurface + "_z"))
    iU0 = utils.first(dataItemNameNow, condition=lambda x: x == "U_x")
    iU1 = utils.first(dataItemNameNow, condition=lambda x: x == "U_y")
    iU2 = utils.first(dataItemNameNow, condition=lambda x: x == "U_z")
    Kasn = (dataItemNow[:, iN0] * dataItemNow[:, iU0] +
            dataItemNow[:, iN1] * dataItemNow[:, iU1] + dataItemNow[:, iN2] *
            dataItemNow[:, iU2]) * 0.5431198 / 1000 / 1.1344992
    Kast = Kas - Kasn
    Ka = (dataItemNow[:,
                      utils.first(dataItemNameNow,
                                  condition=lambda x: x ==
                                  ("Ks_sd" + isoSdSurface))]
          ) * 0.5431198 / 1000 / 1.1344992
    dataItemNow = np.concatenate((dataItemNow, Kas[:, np.newaxis]), axis=1)
    dataItemNow = np.concatenate((dataItemNow, Kasn[:, np.newaxis]), axis=1)
    dataItemNow = np.concatenate((dataItemNow, Kast[:, np.newaxis]), axis=1)
    dataItemNow = np.concatenate((dataItemNow, Ka[:, np.newaxis]), axis=1)

    SdSL = (dataItemNow[:,
                        utils.first(dataItemNameNow,
                                    condition=lambda x: x ==
                                    ("sd_sd" + isoSdSurface))] *
            dataItemNow[:,
                        utils.
                        first(dataItemNameNow, condition=lambda x: x == "rho")]
            ) / 0.43165866 / 1.1344992
    dataItemNow = np.concatenate((dataItemNow, SdSL[:, np.newaxis]), axis=1)

    dataItemNameNow.extend(
        ["S_L-cur", "d_f-cur", "Ka_s", "Ka_sn", "Ka_st", "Ka", "S_d-S_L"])

    return dataItemNow, dataItemNameNow


if __name__ == "__main__":
    # Print summary of current input.
    print('*******************************')
    plt.style.use(
        '/mnt/d/Work/pythonscripts/matplotlib/styles/science.mplstyle')
    sns.set_style("whitegrid", rc=plt.rcParams)

    # pathRootU0 = os.path.abspath(
    #     '/mnt/m/OpenFOAM-1712-Diff/H2Ig-phi_05.10-dilution_00.00/' +
    #     'IgnitionMIE_Upwind/H2Ig_V0.0_L40.0T2_D4.0E-4_T2.0E-4_C0.0E-03_PLOT/' +
    #     'S9.0E10/postProcessing')
    pathRootU4C0 = os.path.abspath(
        '/mnt/m/OpenFOAM-1712-Diff/H2Ig-phi_05.10-dilution_00.00/' +
        'IgnitionMIE_Upwind/H2Ig_V4.0_L20.0T2_D4.0E-4_T2.0E-4_C0.0E-03_PLOT/' +
        'S9.0E10/postProcessing')
    pathRootU4C5 = os.path.abspath(
        '/mnt/m/OpenFOAM-1712-Diff/H2Ig-phi_05.10-dilution_00.00/' +
        'IgnitionMIE_Upwind/H2Ig_V4.0_L20.0T2_D4.0E-4_T2.0E-4_C5.0E-03_PLOT/' +
        'S9.0E10/postProcessing')

    pathToSave = os.path.abspath(
        '/mnt/m/OpenFOAM-1712-Diff/H2Ig-phi_05.10-dilution_00.00/' +
        'tecplot/IgnitionMIEUpwind-DifferentPos/' +
        'H2IgMIEUpwind_L20.0T2_D4.0E-4_T2.0E-4_U4_C0C5/' + 'pythonPDF')

    isoSdSurface = 'H2O'
    isoSdSurfaceValue = '0.082377'
    timeNow = '0.0015'
    keepDataOption = 2
    itemNow = "iso" + isoSdSurface + isoSdSurfaceValue

    # dataItemNowU0, dataItemNameNowU0 = obtainCombinedData(
    #     pathRootU0, timeNow, isoSdSurface, isoSdSurfaceValue, keepDataOption)
    dataItemNowU4C0, dataItemNameNowU4C0 = obtainCombinedData(
        pathRootU4C0, timeNow, isoSdSurface, isoSdSurfaceValue, keepDataOption)
    dataItemNowU4C5, dataItemNameNowU4C5 = obtainCombinedData(
        pathRootU4C5, timeNow, isoSdSurface, isoSdSurfaceValue, keepDataOption)

    if (not os.path.exists(os.path.join(pathToSave, "PDF_" + itemNow))):
        os.makedirs(os.path.join(pathToSave, "PDF_" + itemNow))
    if (not os.path.exists(os.path.join(pathToSave, "PDFPLOT_" + itemNow))):
        os.makedirs(os.path.join(pathToSave, "PDFPLOT_" + itemNow))

    toPlotLabel = [
        "S_{L} {\\nabla} \cdot n", "\delta_{f} {\\nabla} \cdot n", "Ka_{s}",
        "Ka_{s, n}", "Ka_{s, t}", "Ka", "S_{d}^{*}/S_{L}"
    ]

    for i in range(dataItemNowU4C0.shape[1] - 7, dataItemNowU4C0.shape[1]):
        varName = dataItemNameNowU4C0[i]
        # Writing of combined data for each time
        fileNow = os.path.join(pathToSave, "PDF_" + itemNow, timeNow + ".dat")
        with open(fileNow, 'w') as f:
            otFileHeader = 'TITLE="PDF of {0} = {1} at {2} s."\r\n'.format(
                itemNow[0:len("iso" + isoSdSurface)],
                itemNow[len("iso" + isoSdSurface):], timeNow)
            otFileVariables = 'VARIABLES="Var" "Density"\r\n'
            f.writelines([otFileHeader, otFileVariables])

            otZoneHeader = 'ZONE T="{2}" I={0} J={1}\r\n'.format(
                500, 1, varName +
                " U4C0") + 'ZONETYPE = Ordered, DATAPACKING = BLOCK\r\n'
            f.writelines([otZoneHeader])
            density = scipy.stats.gaussian_kde(dataItemNowU4C0[:, i])
            dist = np.max(dataItemNowU4C0[:, i]) - np.min(dataItemNowU4C0[:, i])
            var = np.linspace(
                np.min(dataItemNowU4C0[:, i]) - dist * 2,
                np.max(dataItemNowU4C0[:, i]) + dist * 2, 500)
            np.savetxt(f, var, delimiter=',', newline='\n')
            np.savetxt(f, density(var), delimiter=',', newline='\n')

            otZoneHeader = 'ZONE T="{2}" I={0} J={1}\r\n'.format(
                500, 1, varName +
                " U4C5") + 'ZONETYPE = Ordered, DATAPACKING = BLOCK\r\n'
            f.writelines([otZoneHeader])
            density = scipy.stats.gaussian_kde(dataItemNowU4C5[:, i])
            dist = np.max(dataItemNowU4C5[:, i]) - np.min(dataItemNowU4C5[:, i])
            var = np.linspace(
                np.min(dataItemNowU4C5[:, i]) - dist * 2,
                np.max(dataItemNowU4C5[:, i]) + dist * 2, 500)
            np.savetxt(f, var, delimiter=',', newline='\n')
            np.savetxt(f, density(var), delimiter=',', newline='\n')

        ## Plot using the seaborn
        if (not os.path.exists(
                os.path.join(pathToSave, "PDFPLOT_" + itemNow, timeNow))):
            os.makedirs(os.path.join(pathToSave, "PDFPLOT_" + itemNow,
                                     timeNow))

        pdfPlotFileNow = os.path.join(pathToSave, "PDFPLOT_" + itemNow,
                                      timeNow, dataItemNameNowU4C0[i] + ".svg")
        # g = plt.figure()
        # ax = sns.kdeplot(dataItemNowU4C0[:, i],
        #                  label="$a_g = 400 s^{{-1}}$, mean = {0:.4f}".format(
        #                      np.mean(dataItemNowU4C0[:, i])))
        # ax = sns.kdeplot(dataItemNowU4C5[:, i],
        #                  label="$a_g = 800 s^{{-1}}$, mean = {0:.4f}".format(
        #                      np.mean(dataItemNowU4C5[:, i])))
        # ybound = ax.get_ylim()
        # plt.close(g)
        g = plt.figure()
        # ax = sns.kdeplot(
        #     dataItemNowU0[:, i],
        #     color='red',
        #     label="$a_g = 0 s^{{-1}}$"
        #     #  + "mean = {0:.4f}".format(np.mean(dataItemNowU0[:, i]))
        # )
        ax = sns.kdeplot(
            dataItemNowU4C0[:, i],
            color='black',
            label="$\\bf{x_{{ign}} =  0\ mm}$"
            # + ", mean = {0:.4f}".format(np.mean(dataItemNowU4C0[:, i]))
        )
        ax = sns.kdeplot(
            dataItemNowU4C5[:, i],
            color='red',
            label="$\\bf{x_{{ign}} =  5\ mm}$"
            # + ", mean = {0:.4f}".format(np.mean(dataItemNowU4C5[:, i]))
        )
        # ax.set_ylim(ybound)
        legend = ax.legend(loc='best')
        # ax = sns.rugplot(dataItemNowU4C0[:, i], height=.2, color='red')
        # ax = sns.rugplot(dataItemNowU4C0[:, i], height=.15, color='black')
        # ax = sns.rugplot(dataItemNowU4C5[:, i], height=.1, color='red')
        ax.set_xlabel("$\\bf{" + toPlotLabel[i - dataItemNowU4C0.shape[1] + 7] + "}$")
        ax.set_ylabel("\\bf{PDF}")
        # g = sns.displot(dataItemNow[:, i], kind="kde", rug=True, height=4, aspect=1.25)
        # g.despine=False
        # g.set_axis_labels(dataItemNameNow[i], "Density")
        # g.set_titles(col_template=timeSeries[iTimeSeriesN]+" s")
        g.tight_layout()
        plt.show()
        g.savefig(pdfPlotFileNow)
        plt.close(g)
