#! /usr/bin/env python
# -*- coding: utf-8 -*-
import cantera as ct
import numpy as np
import sys, os
import matplotlib.pylab as plt

plt.switch_backend('agg')
plt.style.use('/mnt/d/Work/pythonscripts/matplotlib/styles/science.mplstyle')

sys.path.append(os.getcwd() + "/..")
sys.path.append(os.getcwd() + "/.")
import utils

if __name__ == "__main__":
    rxnmech = r'chem/H2He/H2He.xml'  # r'chem/H2He/H2He.xml'
    gas = ct.Solution(rxnmech)

    # Initial temperature and pressure
    Tin, Pin = 298.0, 1.0 * ct.one_atm

    # Initial fuel and oxidier profile
    fuel, oxidizer, diluent = 'H2', 'O2', 'Ar'
    iFuel, iOxid = gas.species_index(fuel), gas.species_index(oxidizer)
    iDilu, iN2 = gas.species_index(diluent), gas.species_index('N2')

    MWF, MWO = gas.molecular_weights[iFuel], gas.molecular_weights[iOxid]
    MWD, MWN2 = gas.molecular_weights[iDilu], gas.molecular_weights[iN2]

    nuOF = gas.stoich_air_fuel_ratio(fuel, oxidizer,
                                     basis='mole')  # Mass based
    stoichOF = nuOF * MWF / MWO  # stoichOF = 0.5  # Mole based
    airN2O2 = 3.76  # Mole based

    # Given X_H2 in fuel stream, obtain corresponding properties
    XFF = 0.2
    XFD = 1 - XFF
    YFF = XFF * MWF / (XFF * MWF + XFD * MWD)
    YFD = 1 - YFF

    XOO = 0.21
    XON2 = 1 - XOO
    YOO = XOO * MWO / (XOO * MWO + XON2 * MWN2)
    YON2 = 1 - YOO

    # Obtain the initial mixture strength and stoichiometric mixture fraction
    phiS = nuOF * YFF / YOO
    Zst = 1 / (1 + phiS)
    qN2 = stoichOF * MWO * YON2 / (MWN2 * YOO)
    qDilu = MWF * YFD / (MWD * YFF)
    # qDF = MWF * YFD / (MWD * YFF) + stoichOF * MWO * YOD / (MWD * YOO)
    # qDF = MWF / MWD * (1 / (Zst * YFF) - 1 - nuOF)  # Equivalently

    data_directory = 'data/EdgeFlameFuelDiluXAir/' + fuel + '-' + diluent + '-Air/'
    if not os.path.exists(data_directory):
        os.makedirs(data_directory)

    # Compute the laminar flame speed of stoichiometric mixture
    X = np.zeros(gas.n_species)
    X[[iFuel, iOxid, iDilu, iN2]] = 1, stoichOF, qDilu, qN2
    gas.TPX = Tin, Pin, X
    gas.transport_model = 'Mix'
    LeX = gas.thermal_conductivity / (gas.density_mass * gas.cp_mass *
                                      gas.mix_diff_coeffs)
    rhou = gas.density
    gas.equilibrate('HP')
    rhob = gas.density
    print('// {0} : {1} : {2} : {3} = {4:.4f} : {5:.4f} : {6:.4f} : {7:.4f}'.
          format(fuel, diluent, oxidizer, "N2", 1, qDilu, stoichOF, qN2))
    print('// LeF = {0:.3f}, LeO = {1:.3f}, phou = {2:.3f}, rhob = {3:.3f}.'.
          format(LeX[iFuel], LeX[iOxid], rhou, rhob))
    print('// alpha = {0:.3f}, Tf = {1:.3f}.'.format(rhou / rhob, gas.T))

    gas.TPX = Tin, Pin, X
    initial_grid = 5 * np.array(
        [0.0, 0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005, 0.01], 'd')  # m
    width = 20.e-3  # 10mm wide
    tol_ss = [1.0e-6, 1.0e-10]  # [rtol atol] for steady-state problem
    tol_ts = [1.0e-6, 1.0e-10]  # [rtol atol] for time stepping
    loglevel = 0  # amount of diagnostic output (0 to 5)
    refine_grid = True  # True to enable refinement, False to disable
    sim = ct.FreeFlame(gas, initial_grid)
    sim.flame.set_steady_tolerances(default=tol_ss)
    sim.flame.set_transient_tolerances(default=tol_ts)
    sim.soret_enabled = False
    gas.transport_model = 'Mix'
    sim.transport_model = 'Mix'
    sim.energy_enabled = True
    sim.max_grid_points = 1E+5
    sim.inlet.X, sim.inlet.T = X, Tin
    sim.set_max_jac_age(50, 50)  # Max number of times the Jacobian
    sim.set_time_step(0.1e-06, [2, 5, 10, 20, 80])  # Time steps (s)
    sim.set_refine_criteria(ratio=3.0, slope=0.005, curve=0.005, prune=0.005)
    sim.solve(loglevel, refine_grid)
    sim.write_csv(
        data_directory +
        'phiS_{0:.2f}-{1}-XH2_{2:.2f}-Premixed.csv'.format(phiS, diluent, XFF),
        species='Y')
    Sc = utils.computeConsumptionSpeed(sim)
    deltaF = utils.computeFlameThickness(sim)
    rho = np.zeros_like(sim.grid)
    for n in range(sim.flame.n_points):
        sim.set_gas_state(n)
        rho[n] = sim.gas.density_mass
    print(
        '// SL = {0:.4f}, Sc = {1:.4f}, deltaF = {2:.4f} mm, rhou = {3:.4f}, rhob = {4:.4f}'
        .format(sim.velocity[0], Sc, deltaF * 1000, rho[0], rho[-1]))
    # sim.show_stats
    dataSaveDir = data_directory + 'phiS_{0:.2f}-{1}-XH2_{2:.2f}-Premixed'.format(
        phiS, diluent, XFF)
    utils.plotFlamePNG(sim, gas, dataSaveDir)

    # For fuel and oxider stream
    gasF = ct.Solution(rxnmech)
    Y1 = np.zeros(gas.n_species)
    Y1[[iFuel, iDilu]] = ('%.4f' % YFF), ('%.4f' % YFD)
    gasF.TPY = Tin, Pin, Y1
    gasF.transport_model = 'Mix'
    rhoF = gasF.density
    LeF = gasF.thermal_conductivity / (gasF.density_mass * gasF.cp_mass *
                                       gasF.mix_diff_coeffs[iFuel])

    gasO = ct.Solution(rxnmech)
    Y2 = np.zeros(gas.n_species)
    Y2[[iOxid, iN2]] = ('%.4f' % YOO), ('%.4f' % YON2)
    gasO.TPY = Tin, Pin, Y2
    gasO.transport_model = 'Mix'
    rhoO = gasO.density
    LeO = gasO.thermal_conductivity / (gasO.density_mass * gasO.cp_mass *
                                       gasO.mix_diff_coeffs[iOxid])

    print('// phiS = {0}, Z_st = {1}'.format(phiS, Zst))
    print(
        '// Fuel: YFF = {0:.4f}, YFD = {1:.4f}, rhoF = {2:.4f}, LeF = {3:.4f}'.
        format(YFF, YFD, rhoF, LeF))
    print(
        '// Oxid: YOO = {0:.4f}, YON2 = {1:.4f}, rhoO = {2:.4f}, LeO = {3:.4f}'
        .format(YOO, YON2, rhoO, LeO))

    # Generate Output
    width = 15.e-03
    agN = 100
    # Considering equal momentum, rho U ^2 = const
    UO = agN * width / 4
    UF = UO * np.sqrt(rhoO / rhoF)
    # UF = agN / (2 * (1 + np.sqrt(rhoO / rhoF)) / width), UO = UF
    print('// rho UF UF = rho UO UO, a_g = {0}, width = {1}'.format(
        agN, width))
    print('// UF = {0:.6f}, UO = {1:.6f}'.format(UF, UO))
    print('// rhoF = {0:.6f}, rhoO = {1:.6f}'.format(rhoF, rhoO))

    # Compute the extinction strain rates

    # For Counterflow configuration
    # Comute the Extinction strain rates
    exp_d_a, exp_u_a = -1. / 2., 1. / 2.

    # Set normalized initial strain rate
    alpha = [1.]
    # Initial relative strain rate increase, refinement factor, limit factor
    delta_alpha, delta_alpha_factor, delta_alpha_min = 2., 2., .001
    # Limit of the Temperature decrease
    delta_T_min = 1  # K

    # Indicator of iteration and the latest flame still burning
    iter, iter_last_burning = 0, 0
    # Init iter
    iter_UO = 0.05  # m/s
    delta_iter_UO = 0.05
    delta_iter_UO_factor = 2.
    delta_iter_UO_min = 0.001

    restartFlag = True
    gas = ct.Solution(rxnmech)
    gas.TP = Tin, Pin
    gas.transport_model = 'Mix'
    iter_Tmax, iter_pos_f, iter_ag = utils.computeCFDiffusionFlame(
        gas,
        Pin,
        Tin,
        Tin,
        Y1,
        Y2,
        iter_UO,
        width=width,
        loglevel=0,
        restartFlag=False,
        pathRootSave=dataSaveDir,
        velocityBalanceOption='momentum')
    T_maxArray, ag_Array = [iter_Tmax], [iter_ag]
    UO_Array, xf_Array = [iter_UO], [iter_pos_f]

    while True:
        iter += 1
        alpha.append(alpha[iter_last_burning] + delta_alpha)
        strain_factor = alpha[-1] / alpha[iter_last_burning]
        # # Update axial_velocity
        # iter_UO *= strain_factor**exp_u_a  # m/s
        iter_UO += delta_iter_UO
        gas = ct.Solution(rxnmech)
        gas.TP = Tin, Pin
        gas.transport_model = 'Mix'
        try:
            iter_Tmax, iter_pos_f, iter_ag = utils.computeCFDiffusionFlame(
                gas,
                Pin,
                Tin,
                Tin,
                Y1,
                Y2,
                iter_UO,
                width=width,
                loglevel=0,
                restartFlag=restartFlag,
                pathRootSave=dataSaveDir,
                velocityBalanceOption='momentum')
            # Runtime progress output
            print(
                'Flame established at iter = {0} with iter_UO = {1}, flame temperature = {2}'
                .format(iter, iter_UO, iter_Tmax))
            # Flame still burning, so go to next strain rate
            iter_last_burning = iter
            T_maxArray.append(iter_Tmax)
            ag_Array.append(iter_ag)
            UO_Array.append(iter_UO)
            xf_Array.append(iter_pos_f)
        except Exception as e:
            # Corresponding to failed ignition
            # Procedure if flame extinguished but abortion criterion is not satisfied
            print(e)
            print(
                'Flame extinguished at iter = {0} with iter_UO = {1}.'.format(
                    iter, iter_UO))
            # Reduce relative strain rate increase
            iter_UO -= delta_iter_UO
            delta_iter_UO /= delta_iter_UO_factor
            # iter_UO = UO_Array[-1]
            delta_alpha = delta_alpha / delta_alpha_factor
        # If the temperature difference is too small and the minimum relative
        # strain rate increase is reached, abort
        # if ((delta_alpha < delta_alpha_min)):
        if ((delta_iter_UO_factor < delta_iter_UO_min)):
            print(
                'Flame extinguished at iter = {0} (delta_iter_UO_factor = {1}).'
                .format(iter,
                        delta_iter_UO_factor), 'Abortion criterion satisfied.')
            break
        if ((abs(T_maxArray[-1] - T_maxArray[-2]) < delta_T_min)):
            print('Delta T factor satisfied')
            break
    # Print some parameters at the extinction point
    print('------------------------------------------------------------------')
    print('Parameters at the extinction point:')
    print('Peak temperature T={0:4.0f} K'.format(T_maxArray[-1]))
    print('Global strain rate a_g={0:.4e} 1/s'.format(ag_Array[-1]))
    print('Oxidizer inlet velocity U_O={0:.4f} m/s'.format(UO_Array[-1]))

    # Plot the maximum temperature over the maximum axial velocity gradient
    import matplotlib.pyplot as plt
    plt.figure()
    plt.semilogx(ag_Array, T_maxArray)
    plt.xlabel(r'$a_{g}$ [1/s]')
    plt.ylabel(r'$T_{max}$ [K]')
    plt.tight_layout()
    plt.savefig(dataSaveDir + '/figure_T_max_a_g.png')
    np.savetxt(dataSaveDir + '/T_max_a_g.csv',
               np.stack((ag_Array, UO_Array, T_maxArray, xf_Array), axis=-1),
               delimiter=',')
