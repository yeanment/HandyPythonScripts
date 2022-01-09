#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This example computes the planar laminar flame speed of premixtures.
"""

import cantera as ct
import numpy as np
import matplotlib.pylab as plt 

def computeConsumptionSpeed(sim):
    Tb = max(sim.T)
    Tu = min(sim.T)
    rho_u = max(sim.density)
    integrand = sim.heat_release_rate/sim.cp
    I = np.trapz(integrand, sim.grid)
    Sc = I/(Tb - Tu)/rho_u
    return Sc

def computeFlameThickness(sim):
    Tb = max(sim.T)
    Tu = min(sim.T)
    dT = np.gradient(sim.T, sim.grid)
    deltaf = (Tb - Tu)/np.max(np.abs(dT))
    return deltaf

if __name__ == "__main__":
    p = 1.0 * ct.one_atm  # pressure
    T_in = 298.0  # inlet temperature
    rxnmech = './chem/Li_CXY/H2Li.xml'  # reaction mechanism file
    phi = 5.1
    dilution_ratio = 0.0
    stoich_O2_H2 = 0.5
    air_N2_O2_ratio = 0.79/0.21
    X_H2 = (1-dilution_ratio)*phi/(phi + stoich_O2_H2 + stoich_O2_H2*air_N2_O2_ratio)
    X_O2 = (1-dilution_ratio)*stoich_O2_H2/(phi + stoich_O2_H2 + stoich_O2_H2*air_N2_O2_ratio)
    X_N2 = (1-dilution_ratio)*stoich_O2_H2*air_N2_O2_ratio/(phi + stoich_O2_H2 + stoich_O2_H2*air_N2_O2_ratio)+dilution_ratio
    comp = 'H2:%.4f, O2:%.4f, N2:%.4f' % (X_H2, X_O2, X_N2) # premixed gas composition
    # width = 7.50e-3 # m
    pathRootSave = './data_/phi_%05.2f-dilution_%05.2f' % (phi, dilution_ratio)

    gas = ct.Solution(rxnmech)
    gas.TPX = T_in, p, comp

    initial_grid = 5*np.array([0.0, 0.001, 0.01, 0.02, 0.029, 0.03],'d')/3 # m
    width = 18.e-3 # 18mm wide
    print("End of init.")

    #Set tolerance properties
    tol_ss    = [1.0e-5, 1.0e-10]    # [rtol atol] for steady-state problem
    tol_ts    = [1.0e-5, 1.0e-10]    # [rtol atol] for time stepping
    loglevel  = 1                    # amount of diagnostic output (0 to 5)	
    refine_grid = True               # True to enable refinement, False to disable 

    sim = ct.FreeFlame(gas, initial_grid)
    sim.flame.set_steady_tolerances(default=tol_ss)
    sim.flame.set_transient_tolerances(default=tol_ts)

    sim.inlet.X = comp
    sim.inlet.T = T_in
    sim.soret_enabled = False
    sim.energy_enabled = True
    sim.transport_model = 'Mix'
    sim.set_max_jac_age(50, 50)      # Max number of times the Jacobian before re-evaluated
    sim.set_time_step(0.1e-06, [2, 5,10, 20, 80]) # Time steps (s)
    sim.set_refine_criteria(ratio = 3.0, slope = 0.02, curve = 0.02, prune = 0.01)
    sim.solve(loglevel, refine_grid)
    
    rnet = sim.net_production_rates   # Unités Kmol/m3/s
    hnet = sim.heat_release_rate      # Unités W/m3
    print ('Laminar flame speed = {0:.4f} m/s'.format(sim.u[0]))
    Sc = computeConsumptionSpeed(sim)
    print("Consumption Speed: {0:.4f} m/s".format(Sc))
    deltaF = computeFlameThickness(sim)
    print("Flame Thickness: {0:.4f} mm".format(deltaF*1000))
    # get density field
    rho = np.zeros_like(sim.grid)
    for n in range(sim.flame.n_points):
        sim.set_gas_state(n)
        rho[n] = sim.gas.density_mass
    print("rho_u: {0:.4f} kg/m^3; rho_b: {1:.4f} kg/m^3".format(rho[0], rho[-1]))
    
    sim.show_stats

    # Plot data
    fig=plt.figure(1)
    # create first subplot - adiabatic flame temperature
    a=fig.add_subplot(221)
    a.plot(sim.flame.grid, sim.T)
    plt.title(r'T vs. Position')
    plt.xlabel(r'Position [m]', fontsize=12)
    plt.ylabel(r'Temperature [K]',fontsize=12)
    a.xaxis.set_major_locator(plt.MaxNLocator(4)) # this controls the number of tick marks on the axis

    # create second subplot - velocity
    b=fig.add_subplot(222)
    b.plot(sim.flame.grid, sim.u, label='$V$')
    plt.title(r'Velocity vs. Position')
    plt.xlabel(r'Position [m]', fontsize=12)
    plt.ylabel(r'Velocity [m/s]',fontsize=12)
    b.xaxis.set_major_locator(plt.MaxNLocator(4)) 

    # create third subplot - rho
    c=fig.add_subplot(223)
    p = np.zeros(sim.flame.n_points,'d')
    for n in range(sim.flame.n_points):
        sim.set_gas_state(n)
        p[n]= gas.density_mass
    c.plot(sim.flame.grid, p)
    plt.title(r'Rho vs. Position')
    plt.xlabel(r'Position [m]', fontsize=12)
    plt.ylabel(r'Rho [kg/m^3]',fontsize=12)
    c.xaxis.set_major_locator(plt.MaxNLocator(4)) 

    # create fourth subplot - specie CH4
    d=fig.add_subplot(224)
    fuel = np.zeros(sim.flame.n_points,'d')
    for n in range(sim.flame.n_points):
        sim.set_gas_state(n)
        fuel[n]= gas.Y[gas.species_index('H2')]
    d.plot(sim.flame.grid, fuel)
    plt.title(r'FUEL vs. Position')
    plt.xlabel(r'Position [m]', fontsize=12)
    plt.ylabel(r'FUEL Mole Fraction',fontsize=12)
    d.xaxis.set_major_locator(plt.MaxNLocator(4))

    #plt.subplots_adjust(left=0, right=1, wspace=0.5,top=0.9)
    plt.tight_layout()
    plt.savefig("test.png")


