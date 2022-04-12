#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Compute ignition delay time for differendt initial temperature.
"""

import cantera as ct
import numpy as np
from scipy.signal import argrelextrema
import sys, os
import csv
import matplotlib.pylab as plt
import argparse

plt.switch_backend('agg')
plt.style.use('/mnt/d/Work/pythonscripts/matplotlib/styles/science.mplstyle')

sys.path.append(os.getcwd() + "/..")
sys.path.append(os.getcwd() + "/.")
import utils

if __name__ == "__main__":

    parse = argparse.ArgumentParser( \
        description='Process some input parameters.')
    parse.add_argument('--Mech', action='store', nargs=1, type=str, \
        default = ['./chem/Li_CXY/H2Li.xml'], help='Cantera mechanism file.')
    parse.add_argument('--Fuel', action='store', nargs=1, type=str, \
        default = ['H2'], help='Fuel for the flame.')
    # parse.add_argument('--Oxidizer', action='store', nargs=1, type=str, \
    #     default = ['Air'], help='Oxidizer for the flame.')
    parse.add_argument('--EquivalenceRatio', action='store', nargs=1, type=float, \
        default = [5.1], help='Fuel:Oxid')
    parse.add_argument('--Pressure', action='store', nargs=1, type=float, \
        default = [1.0], help='Default 1 atm.')

    args = parse.parse_args()
    print(args)

    # rxnmech = 'chem/DME_LuChen2015/DME_LuChen2015.cti'
    rxnmech = args.Mech[0]  # reaction mechanism file
    fuel = args.Fuel[0]
    phi = args.EquivalenceRatio[0]
    gas = ct.Solution(rxnmech)

    gas.set_equivalence_ratio(phi, fuel, 'O2:1.0, N2:3.76', basis='mole')
    comp = gas.mole_fraction_dict()

    Tmin, Tmax, nlgT = 0.5, 2.00, 51  # 1000/T min, max
    lgTi = 1000 / np.linspace(Tmin, Tmax, nlgT)
    tauig, Tequil = np.zeros_like(lgTi), np.zeros_like(lgTi)
    tauig1, tauig2 = np.zeros_like(lgTi), np.zeros_like(lgTi)
    tauig3, tauig4 = np.zeros_like(lgTi), np.zeros_like(lgTi)

    Pin = args.Pressure[0] * ct.one_atm

    pathRootSave = os.path.join(
        './data/IgnitionDelay',
        '{0}-Air-phi{1:05.2f}-P{2}atm'.format(fuel, phi, args.Pressure[0]))
    if not os.path.exists(pathRootSave):
        os.makedirs(pathRootSave)

    ndt = int(1e+6)
    timeArray, tempArray = np.zeros(ndt), np.zeros(ndt)

    for i in range(len(lgTi)):
        Tin = lgTi[i]
        gas.TPX = Tin, Pin, comp
        gas.equilibrate('HP')
        Teq = gas.T
        gas.TPX = Tin, Pin, comp

        r = ct.IdealGasConstPressureReactor(gas)
        sim = ct.ReactorNet([r])

        time = 0.0
        timeEstimate = 100
        n = ndt
        dt = 1.e-7
        for idt in range(ndt):
            time += dt
            sim.advance(time)
            timeArray[idt] = time
            tempArray[idt] = r.T

            if (time > 5e-1):
                dt = 2e-5
            elif (time > 2e-1):
                dt = 1e-5
            elif (time > 1e-1):
                dt = 5e-6
            elif (time > 1e-2):
                dt = 2e-6
            elif (time > 1e-3):
                dt = 1e-6
            else:
                dt = 1e-7

            if (abs(r.T - Teq) < 1e-1 or (Tin < 800 and Teq - r.T < 5)
                    or (Tin < 700 and Teq - r.T < 10)
                    or (Tin < 600 and Teq - r.T < 20)):
                timeArray[idt] = time
                tempArray[idt] = r.T
                n = idt
                # print(time, dt)
                break

        plt.figure()
        plt.plot(timeArray[0:n], tempArray[0:n], '^', color='orange')
        plt.xlabel(r'Time [s]')
        plt.ylabel(r'Temp [K]')
        plt.savefig(os.path.join(pathRootSave, 'Tu{0}.png'.format(Tin)),
                    bbox_inches='tight')
        plt.close()

        dTdt = np.gradient(tempArray[0:n], timeArray[0:n])
        try:
            tau1 = np.where(tempArray[0:n] > Tin + 20)[0][0]
        except IndexError as e:
            print("Not enough time to calculate tau 1.")
            tau1 = -1

        try:
            tau2 = np.where(tempArray[0:n] > Tin + 500)[0][0]
        except IndexError as e:
            print("Not enough time to calculate tau 2.")
            tau2 = -1

        try:
            tau3 = np.where(tempArray[0:n] > Tin + 800)[0][0]
        except IndexError as e:
            print("Not enough time to calculate tau 2.")
            tau3 = -1

        try:
            tau4 = np.where(tempArray[0:n] > Tin + 1000)[0][0]
        except IndexError as e:
            print("Not enough time to calculate tau 2.")
            tau4 = -1

        # print(tau1, tempArray[tau1])
        tau = argrelextrema(dTdt, np.greater, axis=0)[0]
        print('Tu = {0} K: tau =  {1} s ({2} K), tau1 = {3} s ({4} K)'.format(
            Tin, timeArray[np.argmax(dTdt)], tempArray[np.argmax(dTdt)],
            timeArray[tau1], tempArray[tau1]))
        # Save data to
        np.savetxt(os.path.join(pathRootSave, 'Tu{0}.csv'.format(Tin)),
                   np.stack((timeArray[0:n], tempArray[0:n]), axis=-1),
                   delimiter=',')

        tauig[i] = timeArray[np.argmax(dTdt)]
        tauig1[i] = timeArray[tau1]
        tauig2[i] = timeArray[tau2]
        tauig3[i] = timeArray[tau3]
        tauig4[i] = timeArray[tau4]
        Tequil[i] = tempArray[n - 1]

    np.savetxt(os.path.join(pathRootSave, 'Tu-tauig.csv'),
               np.stack((lgTi, tauig, tauig1, tauig2, Tequil), axis=-1),
               delimiter=',')

    #################################################################
    # Plot results
    #################################################################
    # create plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.semilogy(1000 / lgTi, tauig, '^', color='orange')
    ax.semilogy(1000 / lgTi, tauig1, '+', color='blue')
    ax.semilogy(1000 / lgTi, tauig2, '+', color='red')
    ax.semilogy(1000 / lgTi, tauig3, 'o', color='black')
    ax.semilogy(1000 / lgTi, tauig4, 'o', color='purple')
    ax.set_xlabel(r'1000/Temp [1000/K]')
    ax.set_ylabel(r'Autoignition delay [ms]')
    ax.axis([Tmin, Tmax, 1E-8, 100.0])

    ax2 = ax.twiny()
    ticks = ax.get_xticks()
    ax2.set_xticks(ticks)
    ax2.set_xticklabels((1000 / ticks).round(1))
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xlabel(r'Temp [K]')
    plt.grid()
    plt.savefig(os.path.join(pathRootSave, 'Tu-tauig.png'),
                bbox_inches='tight')
