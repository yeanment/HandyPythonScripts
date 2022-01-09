# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
import cantera as ct
import numpy as np
import os
import utils

# reaction_mechanism = r'chem/H2He/H2He.xml'
reaction_mechanism = r'gri30.xml'
gas = ct.Solution(reaction_mechanism)

Tin = 298
Pin= 1*101325

phi=1
stoich_O2_H2 = 2
iH2 = gas.species_index('CH4')
iO2 = gas.species_index('O2')
iN2 = gas.species_index('N2')
X = np.zeros(gas.n_species)
X[iH2] = phi
X[iO2] = stoich_O2_H2
X[iN2] = 3.76 * stoich_O2_H2

gas.TPX = Tin, Pin, X
initial_grid = 5*np.array([0.0, 0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005, 0.01],'d') # m
width = 20.e-3 # 10mm wide
tol_ss    = [1.0e-6, 1.0e-10]    # [rtol atol] for steady-state problem
tol_ts    = [1.0e-6, 1.0e-10]    # [rtol atol] for time stepping
loglevel  = 0                    # amount of diagnostic output (0 to 5)	
refine_grid = True               # True to enable refinement, False to disable 
sim = ct.FreeFlame(gas, initial_grid)
sim.flame.set_steady_tolerances(default=tol_ss)
sim.flame.set_transient_tolerances(default=tol_ts)
sim.soret_enabled = False
sim.inlet.X = X
sim.inlet.T = Tin
sim.energy_enabled = True
sim.max_grid_points = 1E+5
sim.set_max_jac_age(50, 50)      # Max number of times the Jacobian before re-evaluated
sim.set_time_step(0.1e-06, [2, 5,10, 20, 80]) # Time steps (s)
sim.set_refine_criteria(ratio = 3.0, slope = 0.005, curve = 0.005, prune = 0.005)
sim.solve(loglevel, refine_grid)
Sc = utils.computeConsumptionSpeed(sim)
deltaF = utils.computeFlameThickness(sim)
rho = np.zeros_like(sim.grid)
for n in range(sim.flame.n_points):
    sim.set_gas_state(n)
    rho[n] = sim.gas.density_mass
print('// SL = {0:.4f}, Sc = {1:.4f}, deltaF = {2:.4f} mm, rhou = {3:.4f}, rhob = {4:.4f}'.format(
          sim.u[0], Sc, deltaF*1000, rho[0], rho[-1]))
# sim.show_stats
dataSaveDir = data_directory + 'phiS_{0:.2f}-qDF_{1:.2f}-Premixed'.format(phiS, qDF)
utils.plotFlamePNG(sim, gas, dataSaveDir)

diluent_N2_H2_ratio = [x/2.0 for x in range(50)]

for q in diluent_N2_H2_ratio:
    X[iN2] = q
    gas.TPX = Tin, Pin, X
    gas.transport_model = 'Multi'
    DXi = gas.mix_diff_coeffs
    DT = gas.thermal_conductivity/(gas.density_mass*gas.cp_mass)
    LeX = DT/DXi
    rhou = gas.density
    gas.equilibrate('HP')
    rhob = gas.density
    print('q = {6:.1f}: LeF = {0:.3f}, LeO = {1:.3f}, phou = {2:.3f}, rhob = {3:.3f}, alpha = {4:.3f}, Tf = {5:.3f}.'.format(LeX[iH2], LeX[iO2], rhou, rhob, rhou/rhob, gas.T, q))

