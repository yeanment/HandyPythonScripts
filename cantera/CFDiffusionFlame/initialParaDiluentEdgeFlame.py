# -*- coding: utf-8 -*-
import cantera as ct
import numpy as np
import os
import matplotlib.pylab as plt 

plt.switch_backend('agg')

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
    # rxnmech = r'chem/H2He/H2He.xml'
    rxnmech = r'chem/Li_CXY/H2Li.xml'
    gas = ct.Solution(rxnmech)

    # Initial temperature and pressure
    Tin = 298
    Pin= 1.0 * ct.one_atm 

    # Initial fuel and oxidier profile
    fuel = 'H2'
    oxidizer = 'O2'
    diluent = 'N2'
    
    stoichOF = 0.5
    qDF = 4.76
    
    iFuel = gas.species_index(fuel)
    iOxid = gas.species_index(oxidizer)
    iDilu = gas.species_index(diluent)

    data_directory = 'EdgeFlameDNS/' + fuel + '-' + oxidizer + '-' + diluent + '/'
    if not os.path.exists(data_directory):
        os.makedirs(data_directory)
    
    # Obtain 
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
    print('// Fuel : Oxidizer : Diluent = {0} : {1} : {2} = {3:.2f} : {4:.2f} : {5:.2f}'.format( 
             fuel, oxidizer, diluent, 1, stoichOF, qDF))
    print('// LeF = {0:.3f}, LeO = {1:.3f}, phou = {2:.3f}, rhob = {3:.3f}, alpha = {4:.3f}, Tf = {5:.3f}.'.format(LeX[iFuel], LeX[iOxid], rhou, rhob, rhou/rhob, gas.T, qDF))

    gas.TPX = Tin, Pin, X
    initial_grid = 5*np.array([0.0, 0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005, 0.01],'d') # m
    width = 10.e-3 # 10mm wide
    tol_ss    = [1.0e-5, 1.0e-10]    # [rtol atol] for steady-state problem
    tol_ts    = [1.0e-5, 1.0e-10]    # [rtol atol] for time stepping
    loglevel  = 0                    # amount of diagnostic output (0 to 5)	
    refine_grid = True               # True to enable refinement, False to disable 
    sim = ct.FreeFlame(gas, initial_grid)
    sim.flame.set_steady_tolerances(default=tol_ss)
    sim.flame.set_transient_tolerances(default=tol_ts)
    sim.soret_enabled = False
    sim.inlet.X = X
    sim.inlet.T = Tin
    sim.energy_enabled = True
    sim.set_max_jac_age(50, 50)      # Max number of times the Jacobian before re-evaluated
    sim.set_time_step(0.1e-06, [2, 5,10, 20, 80]) # Time steps (s)
    sim.set_refine_criteria(ratio = 3.0, slope = 0.02, curve = 0.02, prune = 0.01)
    sim.solve(loglevel, refine_grid)
    sim.write_csv(data_directory + 'phi_{0:.2f}-qDF_{1:.2f}-Premixed.csv'.format(1, qDF))
    Sc = computeConsumptionSpeed(sim)
    deltaF = computeFlameThickness(sim)
    rho = np.zeros_like(sim.grid)
    for n in range(sim.flame.n_points):
        sim.set_gas_state(n)
        rho[n] = sim.gas.density_mass
    print('// SL = {0:.4f}, Sc = {1:.4f}, deltaF = {2:.4f} mm, rhou = {3:.4f}, rhob = {4:.4f}'.format(
              sim.velocity[0], Sc, deltaF*1000, rho[0], rho[-1]))
    # sim.show_stats
    fig=plt.figure(1)
    a=fig.add_subplot(221)
    a.plot(sim.flame.grid, sim.T)
    plt.title(r'T vs. Position')
    plt.xlabel(r'Position [m]', fontsize=12)
    plt.ylabel(r'Temperature [K]',fontsize=12)
    a.xaxis.set_major_locator(plt.MaxNLocator(4)) # this controls the number of tick marks on the axis
    b=fig.add_subplot(222)
    b.plot(sim.flame.grid, sim.velocity, label='$V$')
    plt.title(r'Velocity vs. Position')
    plt.xlabel(r'Position [m]', fontsize=12)
    plt.ylabel(r'Velocity [m/s]',fontsize=12)
    b.xaxis.set_major_locator(plt.MaxNLocator(4)) 
    c=fig.add_subplot(223)
    c.plot(sim.flame.grid, rho)
    plt.title(r'Rho vs. Position')
    plt.xlabel(r'Position [m]', fontsize=12)
    plt.ylabel(r'Rho [kg/m^3]',fontsize=12)
    c.xaxis.set_major_locator(plt.MaxNLocator(4)) 
    d=fig.add_subplot(224)
    YFuel = np.zeros(sim.flame.n_points,'d')
    for n in range(sim.flame.n_points):
        sim.set_gas_state(n)
        YFuel[n]= gas.Y[iFuel]
    d.plot(sim.flame.grid, YFuel)
    plt.title(r'FUEL vs. Position')
    plt.xlabel(r'Position [m]', fontsize=12)
    plt.ylabel(r'FUEL Mole Fraction',fontsize=12)
    d.xaxis.set_major_locator(plt.MaxNLocator(4))
    plt.tight_layout()
    plt.savefig(data_directory + 'phi_{0:.2f}-qDF_{1:.2f}-Premixed.png'.format(1, qDF))

    # Obtain Mass fraction for fuel and oxidizer
    phi = 0.0
    qDF = 9
    Zst = 1/(1 + phi)
    Zst = 0.9407
    phi = 1/Zst  -1
    UF = 1 # m/s
    MW = gas.molecular_weights
    MWF = MW[iFuel]
    MWO = MW[iOxid]
    MWD = MW[iDilu]
    nuOF = MWO*stoichOF/MWF  # nu = MWO * stochi / MWF
    YOO_YFF = nuOF/phi # Y_OO / Y_FF
    YFF = 1/(1 + (nuOF*(MWO + qDF*MWD) - phi*MWO) /(nuOF*MWF + phi*MWO))
    XFF = YFF/MWF/(YFF/MWF + (1-YFF)/MWD)
    YFF = 1/(Zst*(1 + nuOF + qDF*MWD/MWF))
    XFF = YFF/MWF/(YFF/MWF + (1-YFF)/MWD)
    YOO = YFF * YOO_YFF
    YFD = 1 - YFF
    YOD = 1 - YOO
    MW1_MW2 = (YOO/MWO + YOD/MWD)/(YFF/MWF + YFD/MWD)
    UO_UF = np.sqrt(MW1_MW2) # moment equal
    UO = UO_UF*UF


    gasF = ct.Solution(rxnmech)
    Y1 = np.zeros(gas.n_species)
    Y1[iFuel] = YFF
    Y1[iDilu] = YFD
    gasF.TPY = Tin, Pin, Y1
    gasF.transport_model = 'Multi'
    rhoF = gasF.density
    gasO = ct.Solution(rxnmech)
    Y2 = np.zeros(gas.n_species)
    Y2[iOxid] = YOO
    Y2[iDilu] = YOD
    gasO.TPY = Tin, Pin, Y2
    gasO.transport_model = 'Multi'
    rhoO = gasO.density
    print('// Fuel: YFF = {0:.4f}, YFD = {1:.4f}, rhoF = {2:.4f}, UF = {3:.4f}'.format(YFF, YFD, rhoF, UF))
    print('// Oxid: YOO = {0:.4f}, YOD = {1:.4f}, rhoO = {2:.4f}, UO = {3:.4f}'.format(YOO, YOD, rhoO, UO))

    # For Counterflow configuration
    width = 10.e-3 # 10mm wide
    gas = ct.Solution(rxnmech)
    gas.TP = Tin, Pin
    sim = ct.CounterflowDiffusionFlame(gas, width=width)
    sim.P = Pin  # 1 bar
    sim.fuel_inlet.mdot = rhoF*UF  # kg/m^2/s
    sim.fuel_inlet.Y = Y1
    sim.fuel_inlet.T = Tin  # K
    sim.oxidizer_inlet.mdot = rhoO*UO  # kg/m^2/s
    sim.oxidizer_inlet.Y = Y2
    sim.oxidizer_inlet.T = Tin  # K
    sim.set_refine_criteria(ratio=3.0, slope=0.1, curve=0.2, prune=0.03)
    sim.soret_enabled = False
    sim.energy_enabled = True
    sim.solve(loglevel=0, auto=True)
    file_name = 'phi_{0:.2f}-Zst_{2:05.3f}-qDF_{1:.2f}-UF_{3:06.3f}-UO_{4:06.3f}-Diffusion.xml'.format(phi, qDF, Zst, UF, UO)
    sim.save(data_directory + file_name, name='solution',
        description='Cantera version ' + ct.__version__ +
        ', reaction mechanism ' + rxnmech)
    sim.write_csv(data_directory + file_name[0:-4] + '.csv')
    print('// Tf = {0:.4f}, a_mean = {1:.4e}, a_max = {2:.4e}, a_fuel = {3:.4e}, a_oxid = {4:.4e}, a_stochi = {5:.5e}'.format(
              np.max(sim.T), sim.strain_rate('mean'), sim.strain_rate('max'), sim.strain_rate('potential_flow_fuel'), 
              sim.strain_rate('potential_flow_oxidizer'), sim.strain_rate('stoichiometric', fuel=fuel)))
    # sim.shostat
    fig=plt.figure(2)
    a=fig.add_subplot(221)
    a.plot(sim.flame.grid, sim.T)
    plt.title(r'T vs. Position')
    plt.xlabel(r'Position [m]', fontsize=12)
    plt.ylabel(r'Temperature [K]',fontsize=12)
    a.xaxis.set_major_locator(plt.MaxNLocator(4)) # this controls the number of tick marks on the axis
    b=fig.add_subplot(222)
    b.plot(sim.flame.grid, sim.velocity, label='$V$')
    plt.title(r'Velocity vs. Position')
    plt.xlabel(r'Position [m]', fontsize=12)
    plt.ylabel(r'Velocity [m/s]',fontsize=12)
    b.xaxis.set_major_locator(plt.MaxNLocator(4)) 
    rho = np.zeros_like(sim.grid)
    for n in range(sim.flame.n_points):
        sim.set_gas_state(n)
        rho[n] = sim.gas.density_mass
    c=fig.add_subplot(223)
    c.plot(sim.flame.grid, rho)
    plt.title(r'Rho vs. Position')
    plt.xlabel(r'Position [m]', fontsize=12)
    plt.ylabel(r'Rho [kg/m^3]',fontsize=12)
    c.xaxis.set_major_locator(plt.MaxNLocator(4)) 
    d=fig.add_subplot(224)
    YFuel = np.zeros(sim.flame.n_points,'d')
    for n in range(sim.flame.n_points):
        sim.set_gas_state(n)
        YFuel[n]= gas.Y[iFuel]
    d.plot(sim.flame.grid, YFuel)
    plt.title(r'FUEL vs. Position')
    plt.xlabel(r'Position [m]', fontsize=12)
    plt.ylabel(r'FUEL Mole Fraction',fontsize=12)
    d.xaxis.set_major_locator(plt.MaxNLocator(4))
    plt.tight_layout()
    file_name = 'phi_{0:.2f}-qDF_{1:.2f}-Zst_{2:05.3f}-UF_{3:06.3f}-UO_{4:06.3f}-Diffusion.png'.format(phi, qDF, Zst, UF, UO)
    plt.savefig(data_directory + file_name)

    

