#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
An opposed-flow premixed strained flame in evaluating ESR.
"""

import cantera as ct
import numpy as np
import sys, os
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
    parse.add_argument('--TransportModel', action='store', nargs=1, type=str, \
        default = ['Mix'], help='Transportmodel, e.g. Mix, Multi, UnityLeiws.')

    args = parse.parse_args()
    print(args)

    Pin, Tin = 1.0 * ct.one_atm, 600.0  # 1 atm & 298 K
    rxnmech = args.Mech[0]  # reaction mechanism file
    fuel = args.Fuel[0]
    phi = args.EquivalenceRatio[0]
    gas = ct.Solution(rxnmech)
    gas.TP = Tin, Pin
    gas.set_equivalence_ratio(phi, fuel, 'O2:1.0, N2:3.76', basis='mole')
    comp = gas.mass_fraction_dict()

    pathRootSave = './data/CFPremixedTwinFlame/{0}-Air-phi_{1:05.2f}'.format(
        fuel, phi)
    # Set up the problem
    gas = ct.Solution(rxnmech)
    gas.TPX = Tin, Pin, comp

    # Compute the extinction strain rates
    # Comute the Extinction strain rates
    exp_d_a, exp_u_a = -1. / 2., 1. / 2.
    # Set normalized initial strain rate
    alpha = [1.]
    # Initial relative strain rate increase, refinement factor, limit factor
    delta_alpha, delta_alpha_factor, delta_alpha_min = 5., 2., .001
    # Limit of the Temperature decrease
    delta_T_min = 1  # K

    # Indicator of iteration and the latest flame still burning
    iter, iter_last_burning = 0, 0
    # Init iter
    iter_Uin = 1.0  # m/s
    halfWidth = 10.00e-3  # m
    # Umiform increase
    delta_iter_Uin, delta_iter_Uin_factor = 1.0, 2.0
    delta_iter_Uin_min = 0.001
    # Umiform increase

    # Initial for restartFile
    iter_Tmax, iter_posf, iter_ag, iter_Sc = utils.computeCFPremixedTwinFlame(
        gas,
        Pin,
        Tin,
        iter_Uin,
        halfWidth,
        restartFlag=False,
        loglevel=0,
        pathRootSave=pathRootSave)

    # List of peak temperatures， global strain rate, inlet velocity, flamepos
    T_maxArray, ag_Array, Sc_Array = [iter_Tmax], [iter_ag], [iter_Sc]
    Uin_Array, xf_Array = [iter_Uin], [iter_posf]

    # extinction temperature
    temperature_limit_extinction = 500
    restartFlag = True
    while True:
        iter += 1
        gas.TPX = Tin, Pin, comp
        # Update relative strain rates
        alpha.append(alpha[iter_last_burning] + delta_alpha)
        strain_factor = alpha[-1] / alpha[iter_last_burning]
        # Update axial_velocity
        # iter_Uin *= strain_factor**exp_u_a  # m/s
        iter_Uin += delta_iter_Uin
        # iter_ag = 4 * iter_Uin / (2 * halfWidth)

        try:
            iter_Tmax, iter_posf, iter_ag, iter_Sc = utils.computeCFPremixedTwinFlame(
                gas,
                Pin,
                Tin,
                iter_Uin,
                halfWidth,
                restartFlag=restartFlag,
                loglevel=1,
                pathRootSave=pathRootSave)

            # Runtime progress output
            print(
                'Flame established at iter = {0} with iter_Uin = {1}, T_f = {2}'
                .format(iter, iter_Uin, iter_Tmax))
            # Flame still burning, so go to next strain rate
            iter_last_burning = iter
            T_maxArray.append(iter_Tmax)
            ag_Array.append(iter_ag)
            Uin_Array.append(iter_Uin)
            xf_Array.append(iter_posf)
            Sc_Array.append(iter_Sc)

        except Exception as e:
            # Throw Exception if solution fails
            print('Error: Did not converge at iter =', iter, e)
            print(
                'Flame extinguished at iter = {0} with iter_Uin = {1}.'.format(
                    iter, iter_Uin))
            # Reduce relative strain rate increase
            iter_Uin -= delta_iter_Uin
            delta_iter_Uin /= delta_iter_Uin_factor
            # iter_UO = UO_Array[-1]
            delta_alpha = delta_alpha / delta_alpha_factor
        # if ((delta_alpha < delta_alpha_min)):
        if ((delta_iter_Uin < delta_iter_Uin_min)):
            print(
                'Flame extinguished at iter = {0} (delta_iter_Uin = {1}).'.
                format(iter, delta_iter_Uin), 'Abortion criterion satisfied.')
            break
        # if ((abs(T_maxArray[-1] - T_maxArray[-2]) < delta_T_min)):
        #     print('Delta T factor satisfied')
        #     break

    # Print some parameters at the extinction point
    print('----------------------------------------------------------------')
    print('Parameters at the extinction point:')
    print('Peak temperature T={0:4.0f} K'.format(T_maxArray[-1]))
    print('Global strain rate a_g={0:.2e} 1/s'.format(ag_Array[-1]))
    print('Mixture inlet velocity U_in={0:.4f} m/s'.format(Uin_Array[-1]))

    # Plot the maximum temperature over the maximum axial velocity gradient
    plt.figure()
    plt.semilogx(ag_Array, T_maxArray)
    plt.xlabel(r'$a_{g}$ [1/s]')
    plt.ylabel(r'$T_{max}$ [K]')
    plt.tight_layout()
    plt.savefig(pathRootSave + '/T_max_ag.png')
    np.savetxt(pathRootSave + '/T_max_ag.csv',
               np.stack((ag_Array, Uin_Array, T_maxArray, xf_Array, Sc_Array),
                        axis=-1),
               delimiter=',')
