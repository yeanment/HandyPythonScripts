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

    Pin = 1.0 * ct.one_atm

    dt = 2.e-7
    ndt = 10000000

    timeArray, tempArray = np.zeros(ndt), np.zeros(ndt)

    pathRootSave = os.path.join('./data/IgnitionDelay',
                                '{0}-Air-phi{1:05.2f}'.format(fuel, phi))
    if not os.path.exists(pathRootSave):
        os.makedirs(pathRootSave)

    for i in range(len(lgTi)):
        Tin = lgTi[i]
        gas.TPX = Tin, Pin, comp
        r = ct.IdealGasConstPressureReactor(gas)
        sim = ct.ReactorNet([r])
        time = 0.0
        for n in range(ndt):
            time += dt
            sim.advance(time)
            timeArray[n] = time
            tempArray[n] = r.T
        plt.figure()
        plt.plot(timeArray, tempArray, '^', color='orange')
        plt.xlabel(r'Time [s]')
        plt.ylabel(r'Temp [K]')
        plt.savefig(os.path.join(pathRootSave, 'Tu{0}.png'.format(Tin)),
                    bbox_inches='tight')
        plt.close()

        dTdt = np.gradient(tempArray, timeArray)
        tau = argrelextrema(dTdt, np.greater, axis=0)[0]
        print(
            'Ignition delay time at {0} K is {1} s with temperature of {2} K'.
            format(Tin, timeArray[tau], tempArray[tau]))
        # Save data to
        np.savetxt(os.path.join(pathRootSave, 'Tu{0}.csv'.format(Tin)),
                   np.stack((timeArray, tempArray), axis=-1),
                   delimiter=',')
        
        tauig[i] = timeArray[np.argmax(dTdt)]
        Tequil[i] = tempArray[-1]

    np.savetxt(os.path.join(pathRootSave, 'Tu-tauig.csv'),
               np.stack((lgTi, tauig, Tequil), axis=-1),
               delimiter=',')

    #################################################################
    # Plot results
    #################################################################
    # create plot
    plt.semilogy(1000 / lgTi, tauig, '^', color='orange')
    plt.xlabel(r'Temp [1000/K]')
    plt.ylabel(r'Autoignition delay [ms]')
    plt.axis([Tmin, Tmax, 1E-8, 100.0])
    plt.grid()
    plt.savefig(os.path.join(pathRootSave, 'Tu-tauig.png'),
                bbox_inches='tight')
