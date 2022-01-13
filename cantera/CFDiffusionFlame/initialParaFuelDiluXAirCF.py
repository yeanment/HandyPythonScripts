#! /usr/bin/env python
# -*- coding: utf-8 -*-
import cantera as ct
import numpy as np
import sys, os
import matplotlib.pylab as plt
import argparse
# from .. import utils

plt.switch_backend('agg')
plt.style.use('/mnt/d/Work/pythonscripts/matplotlib/styles/science.mplstyle')

sys.path.append(os.getcwd() + "/..")
sys.path.append(os.getcwd() + "/.")
import utils


if __name__ == "__main__":
    parse = argparse.ArgumentParser( \
        description='Process some input parameters.')
    parse.add_argument('--Fuel', action='store', nargs=1, type=str, \
        default = ['H2'], help='Fuel for the flame.')
    parse.add_argument('--Diluent', action='store', nargs=1, type=str, \
        default = ['N2'], help='Diluent, e.g. HE, AR, N2.')
    parse.add_argument('--XFuel', action='store', nargs=1, type=float, \
        default = [0.2], help='XFuel, Fuel mole fraction in fuel stream.')
    # parse.add_argument('--Zst', action='store', nargs=1, type=float, \
    #     default = [0.5], help='Zst, i.e. stoichiometric mixture fraction.')
    parse.add_argument('--TransportModel', action='store', nargs=1, type=str, \
        default = ['Mix'], help='Transportmodel, e.g. Mix, Multi, UnityLeiws.')

    args = parse.parse_args()
    print(args)

    rxnmech = r'chem/H2_LiDryer2004/H2_LiDryer2004.xml'  # r'chem/H2He/H2He.xml'
    spec = ct.Species.listFromFile(rxnmech)
    rxns = ct.Reaction.listFromFile(rxnmech)
    gas = ct.Solution(thermo='IdealGas',
                      kinetics='GasKinetics',
                      species=spec,
                      reactions=rxns,
                      name='With reactions')
    gasCF = ct.Solution(thermo='IdealGas',
                        kinetics='GasKinetics',
                        species=spec,
                        reactions=(),
                        name='Without reactions')

    # Initial temperature and pressure
    Tin, Pin = 298.0, 1.0 * ct.one_atm

    # Initial fuel and oxidier profile
    fuel, oxidizer, diluent = args.Fuel[0], 'O2', args.Diluent[0]
    iFuel, iOxid = gas.species_index(fuel), gas.species_index(oxidizer)
    iDilu, iN2 = gas.species_index(diluent), gas.species_index('N2')

    MWF, MWO = gas.molecular_weights[iFuel], gas.molecular_weights[iOxid]
    MWD, MWN2 = gas.molecular_weights[iDilu], gas.molecular_weights[iN2]

    nuOF = gas.stoich_air_fuel_ratio(fuel, oxidizer,
                                     basis='mole')  # Mass based
    stoichOF = nuOF * MWF / MWO  # stoichOF = 0.5  # Mole based
    # airN2O2 = 3.76  # Mole based

    # Given XFuel in fuel stream, obtain corresponding properties
    XFF, XFD = args.XFuel[0], 1 - args.XFuel[0]
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

    # # Test for Zst
    Y1 = np.zeros(gas.n_species)
    Y1[[iFuel, iDilu]] = ('%.4f' % YFF), ('%.4f' % YFD)
    Y2 = np.zeros(gas.n_species)
    Y2[[iOxid, iN2]] = ('%.4f' % YOO), ('%.4f' % YON2)
    gas.set_mixture_fraction(Zst, Y1, Y2, basis='mass')
    # gas.mixture_fraction(Y1, Y2)
    zstmixture = gas.mass_fraction_dict()
    gas.Y = zstmixture
    gas.equivalence_ratio()

    dataSaveDir = os.path.join(
        'data/EdgeFlameFuelDiluXAir', '{0}-{1}-Air'.format(fuel, diluent),
        'XH2_{0:.2f}-{1}-Zst_{2:.2f}-Air-{3}-Premixed-CF'.format(
            XFF, diluent, Zst, args.TransportModel[0]))
    if not os.path.exists(dataSaveDir):
        os.makedirs(dataSaveDir)

    # For fuel and oxider stream
    gasF = ct.Solution(rxnmech)
    Y1 = np.zeros(gas.n_species)
    Y1[[iFuel, iDilu]] = ('%.4f' % YFF), ('%.4f' % YFD)
    gasF.TPY = Tin, Pin, Y1
    gasF.transport_model = args.TransportModel[0]
    rhoF = gasF.density
    LeF = gasF.thermal_conductivity / (gasF.density_mass * gasF.cp_mass *
                                       gasF.mix_diff_coeffs[iFuel])

    gasO = ct.Solution(rxnmech)
    Y2 = np.zeros(gas.n_species)
    Y2[[iOxid, iN2]] = ('%.4f' % YOO), ('%.4f' % YON2)
    gasO.TPY = Tin, Pin, Y2
    gasO.transport_model = args.TransportModel[0]
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
    width = 20.e-03
    agN = 40
    # Considering equal momentum, rho U ^2 = const
    UO = agN * width / 4
    UF = UO * np.sqrt(rhoO / rhoF)
    # UF = agN / (2 * (1 + np.sqrt(rhoO / rhoF)) / width), UO = UF
    print('// rho UF UF = rho UO UO, a_g = {0}, width = {1}'.format(
        agN, width))
    print('// UF = {0:.6f}, UO = {1:.6f}'.format(UF, UO))
    print('// rhoF = {0:.6f}, rhoO = {1:.6f}'.format(rhoF, rhoO))
    # exit()

    # For Counterflow configuration
    # Compute the mixing layer thickness
    UOArray = np.logspace(-2, 2, 51)
    agArray = np.zeros_like(UOArray)
    deltamArray = np.zeros_like(UOArray)
    uBalance = 'momentum'
    loglevel = 0

    for i in range(len(UOArray)):
        UO = UOArray[i]
        gasCF.TP = Tin, Pin
        gasCF.transport_model = args.TransportModel[0]
        # Compute the mixing layer thickness
        
        if (uBalance == 'velocity'):
            UF = UO
        elif (uBalance == 'mass'):
            UF = UO * rhoO / rhoF
        elif (uBalance == 'momentum'):
            UF = UO * np.sqrt(rhoO / rhoF)
        else:
            exit("The option is not avialable for uBalance.")

        ag = 2 * UO * (1 + UF / UO * np.sqrt(rhoF / rhoO)) / width
        agArray[i] = ag

        sim = ct.CounterflowDiffusionFlame(gasCF, width=width)
        sim.P = Pin  # 1 bar
        sim.fuel_inlet.mdot = rhoF * UF  # kg/m^2/s
        sim.fuel_inlet.T, sim.fuel_inlet.Y = Tin, Y1

        sim.oxidizer_inlet.mdot = rhoO * UO  # kg/m^2/s
        sim.oxidizer_inlet.T, sim.oxidizer_inlet.Y = Tin, Y2

        sim.max_grid_points = 1e5
        sim.max_time_step_count = 1E+6
        sim.soret_enabled = False
        sim.energy_enabled = True
        sim.transport_model = args.TransportModel[0]
        if (sim.transport_model == 'Multi'):
            sim.soret_enabled = True

        try:
            sim.set_refine_criteria(ratio=3.0, slope=0.05, curve=0.05, prune=0.01)
            sim.solve(loglevel=loglevel, auto=False)

            strPathu = uBalance.capitalize() + "-UO_%06.3f-CF" % (UO)
            pathSave = os.path.join(dataSaveDir, strPathu)
            if not os.path.exists(pathSave):
                os.makedirs(pathSave)
            sim.save(os.path.join(pathSave, 'data.xml'),
                    name='strPathu',
                    description='ag = ' + str(ag) + ', Cold flow.')
            
            # Compute mixture fraction
            Zmix = np.zeros_like(sim.grid)
            for n in range(sim.flame.n_points):
                sim.set_gas_state(n)
                Zmix[n] = sim.gas.mixture_fraction(Y1, Y2, basis='mass')
            
            fig_uZmix = plt.figure(facecolor='white', figsize=[12, 8])
            fig_uZmix.subplots_adjust(right=0.75, top=0.95)
            ax_u = fig_uZmix.subplots()
            ax_Zmix = ax_u.twinx()

            ax_u.plot(sim.grid, sim.velocity, 'b-', label='Axial Velocity')
            ax_Zmix.plot(sim.grid, Zmix, 'r-', label='Mixture Fraction')
            ax_u.set_xlim(sim.grid[0], sim.grid[-1])
            ax_u.set_xlabel('Distance [$m$]')
            ax_Zmix.yaxis.set_major_locator(plt.MaxNLocator(4))
            ax_Zmix.set_ylabel('Mixture Fraction')
            ax_u.set_ylabel('Axial Velocity [$m/s$]')
            ax_Zmix.set_ylabel('Heat Release Rate [$W/m^{{3}}$]')
            line1, label1 = ax_u.get_legend_handles_labels()
            line2, label2 = ax_Zmix.get_legend_handles_labels()
            ax_u.legend(line1 + line2, label1 + label2, loc='best')
            fig_uZmix.tight_layout()
            fig_uZmix.savefig(os.path.join(pathSave, "uZmix-Position.png"))
            plt.close(fig_uZmix)

            dZmix = np.gradient(Zmix, sim.grid)
            deltam = 1. / np.max(np.abs(dZmix))

            print('// velocityBalanceOption: ' + uBalance)
            print('// UF = {0:.6f}, UO = {1:.6f}'.format(UF, UO))
            print('// a_g = {0:.4e}, delta_m = {1:.6e}'.format(ag, deltam))

            agArray[i] = ag
            deltamArray[i] = deltam
        except ct.CanteraError as e:
            print('Error occurred while solving:', e)
            
    # Plot the maximum temperature over the maximum axial velocity gradient
    plt.figure()
    plt.semilogx(agArray, deltamArray * 1000)
    plt.xlabel(r'$a_{g}$ [1/s]')
    plt.ylabel(r'$\delta_{m}$ [mm]')
    plt.tight_layout()
    plt.savefig(dataSaveDir + '/ag_deltam.png')
    np.savetxt(dataSaveDir + '/ag_deltam.csv',
               np.stack((agArray, UOArray, deltamArray), axis=-1),
               delimiter=',')
