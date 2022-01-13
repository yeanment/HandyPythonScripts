#! /usr/bin/env python
# -*- coding: utf-8 -*-
import cantera as ct
import numpy as np
import scipy, scipy.special
import sys, os
import matplotlib.pylab as plt

plt.switch_backend('agg')


def computeConsumptionSpeed(sim):
    Tb = max(sim.T)
    Tu = min(sim.T)
    rho_u = max(sim.density)
    integrand = sim.heat_release_rate / sim.cp
    I = np.trapz(integrand, sim.grid)
    Sc = I / (Tb - Tu) / rho_u
    return Sc


def computeFlameThickness(sim):
    Tb = max(sim.T)
    Tu = min(sim.T)
    dT = np.gradient(sim.T, sim.grid)
    deltaf = (Tb - Tu) / np.max(np.abs(dT))
    return deltaf


def computeStrainRates(sim):
    # Compute the derivative of axial velocity to obtain normal strain rate
    strainRates = np.gradient(sim.velocity, sim.grid)
    # Obtain the location of the max. strain rate upstream of the pre-heat zone.
    # This is the characteristic strain rate
    maxStrLocation = abs(strainRates).argmax()
    diffIndex = np.diff(sim.velocity) < 0
    minVelocityPoint = np.logical_and(
        np.append(False, diffIndex), np.append(np.logical_not(diffIndex),
                                               False))
    # Characteristic Strain Rate = K
    # print(minVelocityPoint)
    if len(strainRates[minVelocityPoint]) == 0:
        return None, None, 0
    print(np.where(minVelocityPoint == True)[0][0])
    strainRatePoint = abs(
        strainRates[:np.where(minVelocityPoint == True)[0][0]]).argmax()
    K = abs(strainRates[strainRatePoint])
    print("computeStrainRate: {0}-{1}-{2}.".format(
        sim.velocity[minVelocityPoint], strainRates[minVelocityPoint], K))
    return strainRates, strainRatePoint, K


def plotFlamePNG(sim, gas, pathSave):
    if not os.path.exists(pathSave):
        os.makedirs(pathSave)

    # fontTickLabel = {'family':'sans-serif', 'weight':'normal', 'size': 10}
    # fontAxisLabel = {'family':'sans-serif', 'weight':'normal', 'size': 12}
    # fontLegendLabel = {'family':'sans-serif', 'weight':'normal', 'size': 10}

    # Plot temperature and velocity profiles
    pathSaveTu = os.path.join(pathSave, "Tu-Position.png")
    fig_Tu = plt.figure(facecolor='white', figsize=[12, 8])
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

    ax_T.plot(sim.grid, sim.T, 'b-', label='Temperature')  #sim.flame.grid
    ax_u.plot(sim.grid, sim.velocity, 'r-', label='Axial Velocity')
    ax_Q.plot(sim.grid, sim.heat_release_rate, 'g-', label='Heat Release Rate')
    ax_T.set_xlim(sim.grid[0], sim.grid[-1])
    ax_T.set_xlabel('Distance [$m$]')
    # this controls the number of tick marks on the axis
    ax_T.yaxis.set_major_locator(plt.MaxNLocator(4))
    ax_T.set_ylabel('Temperature [$K$]')
    ax_u.set_ylabel('Axial Velocity [$m/s$]')
    ax_Q.set_ylabel('Heat Release Rate [$W/m^{{3}}$]')
    # fig_Tu.legend(loc='lower right', bbox_to_anchor=(0.5, 0., 0.5, 0.5))
    line1, label1 = ax_T.get_legend_handles_labels()
    line2, label2 = ax_u.get_legend_handles_labels()
    line3, label3 = ax_Q.get_legend_handles_labels()
    ax_u.legend(line1 + line2 + line3, label1 + label2 + label3, loc='best')
    # fig_Tu.legend(prop=fontLegendLabel)

    fig_Tu.tight_layout()
    # fig_Tu.draw()
    fig_Tu.savefig(pathSaveTu)
    plt.close(fig_Tu)

    # Plot species mole fractions
    pathSaveX = os.path.join(pathSave, "X_i.png")
    fig_X = plt.figure(facecolor='white')
    ax_X = fig_X.subplots()
    ax_X.plot(sim.grid, sim.X[gas.species_index('H2')], label='$H_{{2}}$')
    ax_X.plot(sim.grid, sim.X[gas.species_index('O2')], label='$O_{{2}}$')
    ax_X.plot(sim.grid, sim.X[gas.species_index('H2O')], label='$H_{{2}}O$')
    ax_X.set_xlim(sim.grid[0], sim.grid[-1])
    ax_X.set_xlabel('Distance [$m$]')
    ax_X.set_ylabel('$X_{{i}}$')
    ax_X.legend(loc='best')
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


def computeCFPremixedTwinFlame(gas,
                               Pin,
                               Tin,
                               Uin,
                               halfWidth=10.0e-3,
                               loglevel=1,
                               restartFlag=False,
                               pathRootSave='./data',
                               transportModel='Mix'):
    temperature_limit_extinction = Tin + 10

    class FlameExtinguished(Exception):
        pass

    def interrupt_extinction(t):
        # temperature_limit_extinction=600
        if np.max(sim.T) < temperature_limit_extinction:
            raise FlameExtinguished('Flame extinguished')
        return 0.

    if not os.path.exists(pathRootSave):
        os.makedirs(pathRootSave)

    gas.transport_model = transportModel
    gas.TP = Tin, Pin

    sim = ct.CounterflowTwinPremixedFlame(gas, width=halfWidth)
    fileRestart = os.path.join(pathRootSave, 'restart.xml')
    if restartFlag:
        sim.restore(filename=fileRestart, name='restart', loglevel=1)
    ag = 2 * Uin / halfWidth
    sim.reactants.mdot = gas.density * Uin

    # sim.set_interrupt(interrupt_extinction)
    sim.max_grid_points = 1e5
    sim.soret_enabled = False
    sim.energy_enabled = True
    sim.transport_model = transportModel
    if (sim.transport_model == 'Multi'):
        sim.soret_enabled = True

    try:
        sim.set_refine_criteria(ratio=5.0, slope=0.2, curve=0.2, prune=0.2)
        sim.solve(loglevel=loglevel, refine_grid=True, auto=True)
        print("Data size of run 0:", len(sim.grid))
        sim.set_refine_criteria(ratio=3.0, slope=0.05, curve=0.05, prune=0.02)
        sim.solve(loglevel=loglevel, auto=True)
        print("Data size of run 1:", len(sim.grid))
        T_max = np.max(sim.T)
        if T_max < temperature_limit_extinction:
            raise FlameExtinguished('Flame extinguished')

        pos_flame = sim.grid[np.argmax(sim.heat_release_rate)]
        sim.save(fileRestart, name='restart', description='ag = ' + str(ag))
        # sim.show_solution()
        # sim.show_stats()

        Sc = computeConsumptionSpeed(sim)
        _, _, strain = computeStrainRates(sim)

        strPathu = "Uin_%06.3f" % (Uin)
        pathSave = os.path.join(pathRootSave, strPathu)

        plotFlamePNG(sim, gas, pathSave)
        saveFlameCSV(sim, pathSave)
        sim.save(os.path.join(pathSave, 'data.xml'),
                 name='strPathu',
                 description='ag = ' + str(ag))

        print('// Tf = {0:.4f} K, a_g = {1:.4e} s^-1'.format(
            np.max(sim.T), ag))
        print("// Sc = {0:.2f} cm/s, a = {1:.3f} s^-1.".format(
            Sc * 100, strain))

        return T_max, pos_flame, ag, Sc
    except FlameExtinguished:
        print('Flame extinguished')
        return None
    except ct.CanteraError as e:
        print('Error occurred while solving:', e)
        return None
    return None


def computeCFDiffusionFlame(gas,
                            Pin,
                            Tfuel,
                            Toxid,
                            Yfuel,
                            Yoxid,
                            UO,
                            width=10.0e-3,
                            loglevel=1,
                            restartFlag=False,
                            pathRootSave='./data',
                            velocityBalanceOption='velocity',
                            transportModel='Mix'):
    temperature_limit_extinction = 400

    class FlameExtinguished(Exception):
        pass

    def interrupt_extinction(t):
        # temperature_limit_extinction=600
        if np.max(sim.T) < temperature_limit_extinction:
            raise FlameExtinguished('Flame extinguished')
        return 0.

    if not os.path.exists(pathRootSave):
        os.makedirs(pathRootSave)

    gas.transport_model = transportModel
    T, P, Y = gas.TPY

    gas.TPY = Tfuel, Pin, Yfuel
    rhoF = gas.density

    gas.TPY = Toxid, Pin, Yoxid
    rhoO = gas.density

    gas.TPY = T, P, Y

    if (velocityBalanceOption == 'velocity'):
        UF = UO
    elif (velocityBalanceOption == 'mass'):
        UF = UO * rhoO / rhoF
    elif (velocityBalanceOption == 'momentum'):
        UF = UO * np.sqrt(rhoO / rhoF)
    else:
        exit("The option is not avialable for velocityBalanceOption.")

    ag = 2 * UO * (1 + UF / UO * np.sqrt(rhoF / rhoO)) / width

    sim = ct.CounterflowDiffusionFlame(gas, width=width)
    fileRestart = os.path.join(pathRootSave, 'restart.xml')
    if restartFlag:
        sim.restore(
            filename=fileRestart,
            name='restart',
            loglevel=1,
        )

    sim.P = Pin  # 1 bar
    sim.fuel_inlet.mdot = rhoF * UF  # kg/m^2/s
    sim.fuel_inlet.T, sim.fuel_inlet.Y = Tfuel, Yfuel

    sim.oxidizer_inlet.mdot = rhoO * UO  # kg/m^2/s
    sim.oxidizer_inlet.T, sim.oxidizer_inlet.Y = Toxid, Yoxid

    # sim.set_interrupt(interrupt_extinction)
    sim.max_grid_points = 1e5
    sim.max_time_step_count = 1E+6
    sim.soret_enabled = False
    sim.energy_enabled = True
    sim.transport_model = transportModel
    if (sim.transport_model == 'Multi'):
        sim.soret_enabled = True

    try:
        sim.set_refine_criteria(ratio=5.0, slope=0.2, curve=0.2, prune=0.2)
        sim.solve(loglevel=loglevel, refine_grid=True, auto=True)
        print("Data size of run 0:", len(sim.grid))
        sim.set_refine_criteria(ratio=3.0, slope=0.05, curve=0.05, prune=0.01)
        sim.solve(loglevel=loglevel, auto=True)
        print("Data size of run 1:", len(sim.grid))
        # sim.set_refine_criteria(ratio=3.0, slope=0.02, curve=0.02, prune=0.01)
        # sim.solve(loglevel=loglevel, auto=True)
        # print("Data size of run 2:", len(sim.grid))
        T_max = np.max(sim.T)
        if T_max < temperature_limit_extinction:
            raise FlameExtinguished('Flame extinguished')
        pos_flame = sim.grid[np.argmax(sim.heat_release_rate)]
        sim.save(fileRestart, name='restart', description='ag = ' + str(ag))
        # sim.show_solution()
        # sim.show_stats()

        strainrate = sim.strain_rate('stoichiometric',
                                     fuel='H2',
                                     oxidizer='O2',
                                     stoich=0.5)
        strPathu = velocityBalanceOption.capitalize() + "-UO_%06.3f" % (UO)
        pathSave = os.path.join(pathRootSave, strPathu)
        plotFlamePNG(sim, gas, pathSave)
        saveFlameCSV(sim, pathSave)
        sim.save(os.path.join(pathSave, 'data.xml'),
                 name='strPathu',
                 description='ag = ' + str(ag))
        print('// velocityBalanceOption: ' + velocityBalanceOption)
        print('// UF = {0:.6f}, UO = {1:.6f}'.format(UF, UO))
        print('// Tf = {0:.4f}, a_g = {1:.4e}'.format(np.max(sim.T), ag))
        print(
            '// a_mean = {0:.4e}, a_max = {1:.4e}, a_fuel = {2:.4e}, a_oxid = {3:.4e}, a_stochi = {4:.4e}'
            .format(sim.strain_rate('mean'), sim.strain_rate('max'),
                    sim.strain_rate('potential_flow_fuel'),
                    sim.strain_rate('potential_flow_oxidizer'), strainrate))
        return T_max, pos_flame, ag
    except FlameExtinguished:
        print('Flame extinguished')
        return None
    except ct.CanteraError as e:
        print('Error occurred while solving:', e)
        return None
    return None
