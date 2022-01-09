#! /usr/bin/env python
# -*- coding: utf-8 -*-
from genericpath import isfile
import cantera as ct
import numpy as np
import sys, os
import matplotlib.pylab as plt 

sys.path.append(os.getcwd() + "/..")
sys.path.append(os.getcwd() + "/.")
import utils

if __name__ == "__main__":
    plt.switch_backend('agg')
    rxnmech = r'chem/Li_CXY/H2Li.xml'# r'chem/H2He/H2He.xml'
    gas = ct.Solution(rxnmech)

    # Initial temperature and pressure
    Tin = 298.0
    Pin= 1.0 * ct.one_atm 

    # Initial fuel and oxidier profile
    fuel = 'H2'
    oxidizer = 'O2'
    diluent = 'N2'
    
    stoichOF = 0.5
    airN2O2 = 3.76
    
    iFuel = gas.species_index(fuel)
    iOxid = gas.species_index(oxidizer)
    iDilu = gas.species_index(diluent)

    # Obtain Mass fraction for fuel and oxidizer
    phiS = 0.5 # Initial mixture strength
    Zst = 1/(1 + phiS)
    # Zst = 0.2 
    # phiS = 1/Zst - 1 # Initial mixture strength
    
    MW = gas.molecular_weights
    MWF = MW[iFuel]
    MWO = MW[iOxid]
    MWD = MW[iDilu]
    nuOF = MWO*stoichOF/MWF  # nu = MWO * stochiOF / MWF    
    YOO_YFF = nuOF/phiS # Y_OO / Y_FF

    YOO = MWO/(MWO + airN2O2*MWD)
    YFF = YOO/YOO_YFF
    YFD = 1 - YFF
    YOD = 1 - YOO

    qDF = MWF*YFD/(MWD*YFF) + stoichOF*MWO*YOD/(MWD*YOO)
    qDF = MWF/MWD*(1/(Zst*YFF) - 1 - nuOF)
    MWFD = 1/(YFF/MWF + YFD/MWD) # Fuel stream
    MWOD = 1/(YOO/MWO + YOD/MWD) # Oxidizer stream
    MWFD_MWOD = (YOO/MWO + YOD/MWD)/(YFF/MWF + YFD/MWD) 

    data_directory = 'EdgeFlameDNS/' + fuel + '-Air-' + diluent + '/'
    if not os.path.exists(data_directory):
        os.makedirs(data_directory)
    
    # Compute the laminar flame speed of stochiometrci mixture
    X = np.zeros(gas.n_species)
    X[iFuel] = 1
    X[iOxid] = stoichOF
    X[iDilu] = qDF
    gas.TPX = Tin, Pin, X
    gas.transport_model = 'Multi'
    LeX = gas.thermal_conductivity/(gas.density_mass*gas.cp_mass*gas.mix_diff_coeffs)
    rhou = gas.density
    gas.equilibrate('HP')
    rhob = gas.density
    print('// Fuel : Oxidizer : Diluent = {0} : {1} : {2} = {3:.4f} : {4:.4f} : {5:.4f}'.format( 
            fuel, oxidizer, diluent, 1, stoichOF, qDF))
    print('// LeF = {0:.3f}, LeO = {1:.3f}, phou = {2:.3f}, rhob = {3:.3f}, alpha = {4:.3f}, Tf = {5:.3f}.'.format(
            LeX[iFuel], LeX[iOxid], rhou, rhob, rhou/rhob, gas.T,))

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
    sim.write_csv(data_directory + 'phiS_{0:.2f}-qDF_{1:.2f}-Premixed.csv'.format(phiS, qDF), species='Y')
    Sc = utils.computeConsumptionSpeed(sim)
    deltaF = utils.computeFlameThickness(sim)
    rho = np.zeros_like(sim.grid)
    for n in range(sim.flame.n_points):
        sim.set_gas_state(n)
        rho[n] = sim.gas.density_mass
    print('// SL = {0:.4f}, Sc = {1:.4f}, deltaF = {2:.4f} mm, rhou = {3:.4f}, rhob = {4:.4f}'.format(
              sim.velocity[0], Sc, deltaF*1000, rho[0], rho[-1]))
    # sim.show_stats
    dataSaveDir = data_directory + 'phiS_{0:.2f}-qDF_{1:.2f}-Premixed'.format(phiS, qDF)
    utils.plotFlamePNG(sim, gas, dataSaveDir)
    
    # For fuel and oxider stream
    gasF = ct.Solution(rxnmech)
    Y1 = np.zeros(gas.n_species)
    Y1[iFuel] = ('%.4f' % YFF) # YFF
    Y1[iDilu] = ('%.4f' % YFD) # YFD
    gasF.TPY = Tin, Pin, Y1
    gasF.transport_model = 'Multi'
    rhoF = gasF.density

    gasO = ct.Solution(rxnmech)
    Y2 = np.zeros(gas.n_species)
    Y2[iOxid] = ('%.4f' % YOO) # YOO
    Y2[iDilu] = ('%.4f' % YOD) # YOD
    gasO.TPY = Tin, Pin, Y2
    gasO.transport_model = 'Multi'
    rhoO = gasO.density
    
    print('// phiS = {0}, Z_st = {1}'.format(phiS, Zst))
    print('// Fuel: YFF = {0:.4f}, YFD = {1:.4f}, rhoF = {2:.4f}'.format(YFF, YFD, rhoF))
    print('// Oxid: YOO = {0:.4f}, YOD = {1:.4f}, rhoO = {2:.4f}'.format(YOO, YOD, rhoO))

    # Generate Output
    agN = 100
    UF = agN/(2*(1 + np.sqrt(rhoO/rhoF))/width)
    UO = UF
    print('// UF = UO, a_g = {0}, width = {1}'.format(agN, width))
    print('// UO = {0:.6f}, UO = {1:.6f}'.format(UF, UO))

    # For Counterflow configuration
    # Comute the Extinction strain rates
    exp_d_a = -1. / 2.
    exp_u_a = 1. / 2.
    
    # Set normalized initial strain rate
    alpha = [1.]
    # Initial relative strain rate increase
    delta_alpha = 2.
    # Factor of refinement of the strain rate increase
    delta_alpha_factor = 2.
    # Limit of the refinement: Minimum normalized strain rate increase
    delta_alpha_min = .001
    # Limit of the Temperature decrease
    delta_T_min = 1  # K

    # Iteration indicator
    iter = 0
    # Indicator of the latest flame still burning
    iter_last_burning = 0
    # Init iter
    iter_UF = 0.05 # m/s
    width = 10.e-3 # 10mm wide
    dataSaveDir = data_directory + 'phiS_{0:.2f}-qDF_{1:.2f}-Premixed'.format(
        phiS, qDF) + '-Width{0:.0f}mm'.format(width*1000)
    
    T_max = []
    ag_max = []
    UF = []
    pos_f = [] 
    for item in os.listdir(dataSaveDir):
        dirPath = os.path.join(dataSaveDir, item)
        if(os.path.isfile(dirPath)):
            continue
        gas = ct.Solution(rxnmech)
        gas.TP = Tin, Pin
        gas.transport_model = 'Mix'

        dataXmlPath = os.path.join(dirPath, 'data.xml')
        sim = ct.CounterflowDiffusionFlame(gas, width=width)
        sim.restore(filename=dataXmlPath, name='strPathu', loglevel=1)
        if np.max(sim.T) < 1200:
            continue
        T_max.append(np.max(sim.T))
        UFN = sim.fuel_inlet.mdot/rhoF
        UON = sim.oxidizer_inlet.mdot/rhoO
        UF.append(sim.fuel_inlet.mdot/rhoF)
        ag_max.append(2*(1 + UON/UFN * np.sqrt(rhoO/rhoF))* UFN/width )
        pos_f.append(sim.grid[np.argmax(sim.heat_release_rate)])
    np.savetxt(dataSaveDir + '/T_max_a_g2.csv', np.stack((ag_max, UF, T_max, pos_f), axis = -1), delimiter=',')
