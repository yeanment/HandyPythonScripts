#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
An opposed-flow premixed strained flame in evaluating ESR.
"""

import cantera as ct
import numpy as np
import sys, os, re
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
    parse.add_argument('--Tin', action='store', nargs=1, type=float, \
        default = [300.], help='Fuel:Oxid')
    parse.add_argument('--TransportModel', action='store', nargs=1, type=str, \
        default = ['Mix'], help='Transportmodel, e.g. Mix, Multi, UnityLeiws.')
    parse.add_argument('--InitialFile', action='store', nargs=1, type=str, \
        default = [''], help='InitialFile for ESR estimation.')

    args = parse.parse_args()
    print(args)

    Pin, Tin = 1.0 * ct.one_atm, args.Tin[0]  # 1 atm & 298 K
    rxnmech = args.Mech[0]  # reaction mechanism file
    fuel = args.Fuel[0]
    phi = args.EquivalenceRatio[0]
    gas = ct.Solution(rxnmech)
    gas.TP = Tin, Pin
    gas.set_equivalence_ratio(phi, fuel, 'O2:1.0, N2:3.76', basis='mole')
    comp = gas.mole_fraction_dict()

    pathRootSave = os.path.join(
        './data/CFPremixedTwinFlame',
        '{0}-Air-phi{1:05.2f}-{2}-Tu{3:.1f}'.format(fuel, phi,
                                                    args.TransportModel[0],
                                                    Tin))
    if not os.path.exists(pathRootSave):
        os.makedirs(pathRootSave)

    # Set up the problem
    gas = ct.Solution(rxnmech)
    gas.TPX = Tin, Pin, comp
    print('// {0}/Air: phi = {1}'.format(fuel, phi))
    print('// Mass fraction dict: {0}'.format(gas.mass_fraction_dict()))
    print('// Mole fraction dict: {0}'.format(gas.mole_fraction_dict()))

    gas.transport_model = args.TransportModel[0]
    LeX = gas.thermal_conductivity / (gas.density_mass * gas.cp_mass *
                                      gas.mix_diff_coeffs)
    rhou = gas.density
    gas.equilibrate('HP')
    rhob = gas.density
    print('// LeF = {0:.3f}, LeO = {1:.3f}'.format(
        LeX[gas.species_index(fuel)], LeX[gas.species_index('O2')]))
    print('// phou = {0:.3f}, rhob = {1:.3f}'.format(rhou, rhob))
    print('// alpha = {0:.3f}, Tf = {1:.3f}.'.format(rhou / rhob, gas.T))

    gas.TPX = Tin, Pin, comp
    initial_grid = 5 * np.array(
        [0.0, 0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005, 0.01], 'd')  # m
    width = 20.e-3  # 10mm wide
    tol_ss = [1.0e-6, 1.0e-10]  # [rtol atol] for steady-state problem
    tol_ts = [1.0e-6, 1.0e-10]  # [rtol atol] for time stepping
    loglevel = 1  # amount of diagnostic output (0 to 5)
    refine_grid = True  # True to enable refinement, False to disable
    sim = ct.FreeFlame(gas, initial_grid)
    sim.flame.set_steady_tolerances(default=tol_ss)
    sim.flame.set_transient_tolerances(default=tol_ts)
    sim.soret_enabled = False
    sim.max_time_step_count, sim.max_grid_points = 1E+4, 1E+5
    gas.transport_model = args.TransportModel[0]
    sim.transport_model = args.TransportModel[0]
    if (sim.transport_model == 'Multi'):
        sim.soret_enabled = True
    sim.energy_enabled = True
    sim.inlet.X, sim.inlet.T = comp, Tin
    sim.set_max_jac_age(50, 50)  # Max number of times the Jacobian
    # sim.set_time_step(0.1e-06, [2, 5, 10, 20, 80])  # Time steps (s)
    # sim.set_refine_criteria(ratio=3.0, slope=0.05, curve=0.05, prune=0.02)
    sim.set_refine_criteria(ratio=3.0, slope=0.2, curve=0.2, prune=0.001)
    sim.solve(loglevel, refine_grid)
    sim.write_csv(os.path.join(pathRootSave,
                               '{0}-Air-phi_{1:05.2f}.csv'.format(fuel, phi)),
                  species='Y')
    Sc = utils.computeConsumptionSpeed(sim)
    deltaF = utils.computeFlameThickness(sim)
    rho = np.zeros_like(sim.grid)
    for n in range(sim.flame.n_points):
        sim.set_gas_state(n)
        rho[n] = sim.gas.density_mass
    print('// SL = {0:.4f} m/s, Sc = {1:.4f} m/s, deltaF = {2:.4f} mm'.format(
        sim.velocity[0], Sc, deltaF * 1000))
    print('// rhou = {0:.4f}, rhob = {1:.4f}'.format(rho[0], rho[-1]))
    utils.plotFlamePNG(sim, gas, pathRootSave)

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
    delta_iter_Uin, delta_iter_Uin_factor = 0.1, 2.0
    delta_iter_Uin_min = 0.001
    # Umiform increase

    # Initial for restartFile
    gas.TPX = Tin, Pin, comp
    if (args.InitialFile[0] == ''):
        iter_Tmax, iter_posf, iter_ag, iter_Sc = utils.computeCFPremixedTwinFlame(
            gas,
            Pin,
            Tin,
            iter_Uin,
            halfWidth,
            restartFlag=False,
            loglevel=0,
            pathRootSave=pathRootSave)
    else:
        # Provide initial file for
        with open(args.InitialFile[0]) as f:
            dataName = re.split(r"[\t, ]+", f.readline().strip())
            dataName = [item[1:-1] for item in dataName]
            dataInput = np.loadtxt(f, delimiter=",", skiprows=0)
            dataInput = dataInput[::-1, :]
        iGrid = utils.first(dataName, condition=lambda x: x == 'Points:0')
        # Generate restart file
        sim = ct.CounterflowTwinPremixedFlame(gas,
                                              grid=dataInput[0, iGrid] -
                                              dataInput[:, iGrid])
        locRelative = sim.grid / sim.grid[-1]
        sim.set_profile('T', locRelative, dataInput[:, dataName.index('T')])
        # sim.set_profile('pressure', locRelative, np.ones_like(sim.P) * ct.one_atm)
        sim.set_profile('velocity', locRelative,
                        -dataInput[:, dataName.index('U:0')])
        for item in gas.species_names:
            sim.set_profile(item, locRelative, dataInput[:,
                                                         dataName.index(item)])

        gas.TPX = Tin, Pin, comp
        sim.reactants.mdot = gas.density * iter_Uin
        # sim.set_interrupt(interrupt_extinction)
        sim.max_grid_points = 1e5
        sim.soret_enabled = False
        sim.energy_enabled = True
        sim.transport_model = args.TransportModel[0]
        if (sim.transport_model == 'Multi'):
            sim.soret_enabled = True
        sim.save(args.InitialFile[0][:-3] + 'xml', name='initial restart')
        sim.save(os.path.join(pathRootSave, 'restart.xml'),
                 name='restart')
        iter_Tmax, iter_posf, iter_ag, iter_Sc = utils.computeCFPremixedTwinFlame(
            gas,
            Pin,
            Tin,
            iter_Uin,
            halfWidth,
            restartFlag=True,
            loglevel=1,
            pathRootSave=pathRootSave)

    # List of peak temperatures， global strain rate, inlet velocity, flamepos
    T_maxArray, ag_Array, Sc_Array = [iter_Tmax], [iter_ag], [iter_Sc]
    Uin_Array, xf_Array = [iter_Uin], [iter_posf]

    # extinction temperature
    # temperature_limit_extinction = Tin + 100
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
                loglevel=0,
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
    plt.savefig(os.path.join(pathRootSave, 'T_max_ag.png'))
    np.savetxt(os.path.join(pathRootSave, 'T_max_ag.csv'),
               np.stack((ag_Array, Uin_Array, T_maxArray, xf_Array, Sc_Array),
                        axis=-1),
               delimiter=',')
