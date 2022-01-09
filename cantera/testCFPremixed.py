# -*- coding: utf-8 -*-
"""
An opposed-flow premixed strained flame
"""

import cantera as ct
import numpy as np
import sys, os

# Differentiation function for data that has variable grid spacing Used here to
# compute normal strain-rate
def computeConsumptionSpeed(sim):
    Tb = max(sim.T)
    Tu = min(sim.T)
    rho_u = max(sim.density)
    integrand = sim.heat_release_rate/sim.cp
    I = np.trapz(integrand, sim.grid)
    Sc = I/(Tb - Tu)/rho_u
    return Sc


def computeStrainRates(sim):
    # Compute the derivative of axial velocity to obtain normal strain rate
    strainRates = np.gradient(sim.u, sim.grid)
    # Obtain the location of the max. strain rate upstream of the pre-heat zone.
    # This is the characteristic strain rate
    maxStrLocation = abs(strainRates).argmax()
    diffIndex = np.diff(sim.u) < 0
    minVelocityPoint = np.logical_and(np.append(False, diffIndex), np.append(np.logical_not(diffIndex), False))
    # Characteristic Strain Rate = K
    # print(minVelocityPoint)
    if len(strainRates[minVelocityPoint]) == 0:
        return None, None, 0
    print(np.where(minVelocityPoint == True)[0][0])
    strainRatePoint = abs(strainRates[:np.where(minVelocityPoint == True)[0][0]]).argmax()
    K = abs(strainRates[strainRatePoint])
    print("computeStrainRate: {0}-{1}-{2}.".format(sim.u[minVelocityPoint], strainRates[minVelocityPoint], K))
    return strainRates, strainRatePoint, K


def plotFlamePNG(sim, pathSave):
    import matplotlib.pyplot as plt
    plt.switch_backend('agg')
    if not os.path.exists(pathSave):
        os.makedirs(pathSave)
    else:
        pass
    
    fontTickLabel = {'family':'Tahoma', 'weight':'normal', 'size': 10, }
    fontAxisLabel = {'family':'Tahoma', 'weight':'normal', 'size': 12, }
    fontLegendLabel = {'family':'Tahoma', 'weight':'normal', 'size': 12, }
    
    # Plot temperature and velocity profiles
    pathSaveTu = os.path.join(pathSave, "Tu.png")
    fig_Tu = plt.figure(figsize=(5,4), facecolor='white', dpi=300)
    fig_Tu.subplots_adjust(right=0.75, top=0.95)
    ax_T = fig_Tu.subplots()
    ax_u = ax_T.twinx()
    ax_Q = ax_T.twinx()

    ax_u.spines["right"].set_position(("axes", 1.2))
    ax_u.set_frame_on(True)
    ax_u.patch.set_visible(False)
    for sp in ax_u.spines.values():
        sp.set_visible(False)
    ax_u.spines["right"].set_visible(True)

    ax_T.plot(sim.grid, sim.T, 'b-', label='Temperature')
    ax_u.plot(sim.grid, sim.u, 'r-', label='Axial Velocity')
    ax_Q.plot(sim.grid, sim.heat_release_rate, 'g-', label='Heat Release Rate')
    ax_T.set_xlim(sim.grid[0], sim.grid[-1])
    ax_T.set_xlabel('Distance [m]')
    ax_T.set_ylabel('Temperature [K]')
    ax_u.set_ylabel('Axial Velocity [m/s]')
    ax_Q.set_ylabel('Heat Release Rate [W/$\mathregular{m^3}$]')
    fig_Tu.legend(loc='lower left', bbox_to_anchor=(0.2, 0.2))
    fig_Tu.tight_layout()
    # fig_Tu.draw()
    fig_Tu.savefig(pathSaveTu)
    plt.close(fig_Tu)

    # Plot species mole fractions
    pathSaveX = os.path.join(pathSave, "X_i.png")
    fig_X = plt.figure(figsize=(5,4), facecolor='white', dpi=300)
    ax_X = fig_X.subplots()
    ax_X.plot(sim.grid, sim.X[gas.species_index('H2')], label='$\mathregular{H_2}$')
    ax_X.plot(sim.grid, sim.X[gas.species_index('O2')], label='$\mathregular{O_2}$')
    ax_X.plot(sim.grid, sim.X[gas.species_index('H2O')], label='$\mathregular{H_2O}$')
    ax_X.set_xlim(sim.grid[0], sim.grid[-1])
    ax_X.set_xlabel('Distance (m)')
    ax_X.set_ylabel('$\mathregular{X_i}$')
    fig_X.legend(loc='lower left', bbox_to_anchor=(0.2, 0.2))
    fig_X.tight_layout()
    # fig_X.draw()
    fig_X.savefig(pathSaveX)
    plt.close(fig_X)
    
def saveFlameCSV(sim, pathSave):
    if not os.path.exists(pathSave):
        os.makedirs(pathSave)
    else:
        pass    
    pathSaveCSV = os.path.join(pathSave, "data.csv")
    sim.write_csv(pathSaveCSV)


def computePremixedTwinFlame(premixedGas, axial_velocity, width = 7.50e-3, loglevel=1, restartFlag=False, pathRootSave='./data'):
    temperature_limit_extinction = 500
    def interrupt_extinction(t):
        # temperature_limit_extinction=600
        if np.max(sim.T) < temperature_limit_extinction:
            raise Exception('Flame extinguished')
        return 0.

    if not os.path.exists(pathRootSave):
        os.makedirs(pathRootSave)
    else:
        pass

    # Create the flame simulation object
    sim = ct.CounterflowTwinPremixedFlame(gas=premixedGas, width=width)

    fileRestart = os.path.join(pathRootSave, 'restart.xml')
    if restartFlag:
        sim.restore(filename=fileRestart, name='restart', loglevel=1)
    else:
        pass

    mdot = premixedGas.density * axial_velocity  # kg/m^2/s
    # Set grid refinement parameters
    sim.set_refine_criteria(ratio=3, slope=0.1, curve=0.2, prune=0.02)
    # set the boundary flow rates
    sim.reactants.mdot = mdot
    # set axisymetric
    
    # set the extinction criteria
    # sim.set_interrupt(interrupt_extinction)
    sim.solve(loglevel, auto=True)

    T_max = np.max(sim.T)
    pos_flame = sim.grid[np.argmax(sim.heat_release_rate)]
    if T_max < temperature_limit_extinction:
        return T_max, pos_flame

    sim.save(fileRestart, name='restart', 
            description='Cantera version '+ct.__version__)
    
    # sim.show_solution()
    # sim.show_stats()

    Sc = computeConsumptionSpeed(sim)
    _, _, strain = computeStrainRates(sim)
    # print("Peak temperature: {0:.1f} K".format(T))
    print("Consumption Speed: {0:.2f} cm/s. Strain Rate: {1:.3f} s^-1.".format(Sc*100, strain))

    strPathu = "vu_%06.3f" % (axial_velocity)

    pathSave = os.path.join(pathRootSave, strPathu)

    # Generate plots to see results, if user desires
    if '--plot' in sys.argv:
        plotFlamePNG(sim, pathSave)
    else:
        print('************')
        print('Plotting option not enabled. Re-run script with --plot to see key plots.')
        print('************')
    
    if '--save' in sys.argv:
        saveFlameCSV(sim, pathSave)
    else:
        print('************')
        print('Saving option not enabled. Re-run script with --save to see key plots.')
        print('************')
    
    return T_max, pos_flame, Sc, strain



if __name__ == "__main__":
    p = 1.0 * ct.one_atm  # pressure
    T_in = 298.0  # inlet temperature
    rxnmech = './chem/Li_CXY/H2Li.xml'  # reaction mechanism file
    phi = 1.0
    dilution_ratio = 0.0
    stoich_O2_H2 = 0.5
    air_N2_O2_ratio = 0.79/0.21
    X_H2 = (1-dilution_ratio)*phi/(phi + stoich_O2_H2 + stoich_O2_H2*air_N2_O2_ratio)
    X_O2 = (1-dilution_ratio)*stoich_O2_H2/(phi + stoich_O2_H2 + stoich_O2_H2*air_N2_O2_ratio)
    X_N2 = (1-dilution_ratio)*stoich_O2_H2*air_N2_O2_ratio/(phi + stoich_O2_H2 + stoich_O2_H2*air_N2_O2_ratio)+dilution_ratio
    comp = 'H2:%.4f, O2:%.4f, N2:%.4f' % (X_H2, X_O2, X_N2) # premixed gas composition
    width = 10.00e-3* 2 # m
    pathRootSave = './data/phi_%05.2f-dilution_%05.2f' % (phi, dilution_ratio)
    # Set up the problem
    gas = ct.Solution(rxnmech)
    gas.TPX = T_in, p, comp

    exp_d_a = -1. / 3.
    exp_u_a = 1. / 3.
    
    # Set normalized initial strain rate
    alpha = [1.]
    # Initial relative strain rate increase
    delta_alpha = 5
    # Factor of refinement of the strain rate increase
    delta_alpha_factor = 5.
    # Limit of the refinement: Minimum normalized strain rate increase
    delta_alpha_min = .001
    # Limit of the Temperature decrease
    delta_T_min = 1  # K

    # Iteration indicator
    iter = 0
    # Indicator of the latest flame still burning
    iter_last_burning = 0
    # Init iter
    iter_axialVelocity = 1.0 # m/s
    iter_globalStrainRate = 4*iter_axialVelocity/(2*width)
    # Initial for restartFile
    iter_Tmax, iter_posf, iter_Sc, iter_strain = computePremixedTwinFlame(gas, iter_axialVelocity, width, restartFlag=False, pathRootSave=pathRootSave)
    
    # List of peak temperatures， global strain rate, inlet velocity, flamepos
    T_max = [iter_Tmax]
    a_g = [iter_globalStrainRate]
    v_in = [iter_axialVelocity]
    pos_f = [iter_posf]
    # List of maximum axial velocity gradients
    a_g = [iter_globalStrainRate]
    Sc= [iter_Sc]
    strain = [iter_strain]
    # extinction temperature
    temperature_limit_extinction = 500
    restartFlag = True

    while True:
        iter += 1
        gas.TPX = T_in, p, comp
        # Update relative strain rates
        alpha.append(alpha[iter_last_burning] + delta_alpha)
        strain_factor = alpha[-1] / alpha[iter_last_burning]
        # Update axial_velocity
        iter_axialVelocity += min(0.1, strain_factor ** exp_u_a)# m/s
        iter_globalStrainRate = 4*iter_axialVelocity/(2*width)
        
        try:
            iter_Tmax, iter_posf, iter_Sc, iter_strain  = computePremixedTwinFlame(gas, iter_axialVelocity, width, restartFlag=restartFlag, pathRootSave=pathRootSave)
        except Exception as e:
            # Throw Exception if solution fails
            print('Error: Did not converge at iter =', iter, e)
            iter_Tmax = 0
            iter_posf = 0

        if iter_Tmax > temperature_limit_extinction:
            # Runtime progress output
            print('Flame established at iter = {0} with axial_velocity = {1}, flame temperature = {2}, alpha = {3}, delta_alpha = {4}'.format(
                iter, iter_axialVelocity, iter_Tmax, alpha[iter], delta_alpha))
            # Flame still burning, so go to next strain rate
            iter_last_burning = iter
            T_max.append(iter_Tmax)
            a_g.append(iter_globalStrainRate)
            v_in.append(iter_axialVelocity)
            pos_f.append(iter_posf)
            Sc.append(iter_Sc)
            strain.append(iter_strain)
            # If the temperature difference is too small and the minimum relative
            # strain rate increase is reached, abort
            if ((delta_alpha < delta_alpha_min)):
                print('Flame extinguished at iter = {0}.'.format(iter),
                    'Abortion criterion satisfied.')
                break
            if (iter_axialVelocity < 1.0):
                break
        else:
            # Procedure if flame extinguished but abortion criterion is not satisfied
            print('Flame extinguished at iter = {0} with axial_velocity = {4:.6f}. Restoring iter = {1} with alpha = {2:.5f} with delta_alpha = {3:.5f} '.format(
                iter, iter_last_burning, alpha[iter_last_burning], delta_alpha, iter_axialVelocity))
            # Reduce relative strain rate increase
            iter_globalStrainRate = a_g[-1]
            iter_axialVelocity = v_in[-1]
            delta_alpha = delta_alpha / delta_alpha_factor
            if ((delta_alpha < delta_alpha_min)):
                print('Flame extinguished at iter = {0}.'.format(iter),
                    'Abortion criterion satisfied.')
                break

    # Print some parameters at the extinction point
    print('----------------------------------------------------------------------')
    print('Parameters at the extinction point:')
    print('Peak temperature T={0:4.0f} K'.format(T_max[-1]))
    print('Global strain rate a_g={0:.2e} 1/s'.format(a_g[-1]))
    
    np.savetxt(pathRootSave + 'T_max_a_g.csv', np.stack((a_g, v_in, T_max, pos_f, Sc, strain), axis = -1), delimiter=',')
    # Plot the maximum temperature over the maximum axial velocity gradient
    import matplotlib.pyplot as plt
    plt.switch_backend('egg')
    plt.figure()
    plt.semilogx(a_g, T_max)
    plt.xlabel(r'$a_{g}$ [1/s]')
    plt.ylabel(r'$T_{max}$ [K]')
    plt.tight_layout()
    plt.savefig(pathRootSave + 'figure_T_max_a_g.png')
    # np.savetxt(pathRootSave + 'T_max_a_g.csv', np.stack((a_g, v_in, T_max, pos_f, Sc, strain), axis = -1), delimiter=',')


        


    
    
    