# -*- coding: utf-8 -*-
import cantera as ct
import numpy as np
import os

rxnmech = r'chem/H2He/H2He.xml'
gas = ct.Solution(rxnmech)
# Initial temperature and pressure
Tin, Pin = 298., 1.0 * ct.one_atm

# Initial fuel and oxidier profile
fuel, oxidizer, diluent = 'H2', 'O2', 'He'
iFuel, iOxid = gas.species_index(fuel), gas.species_index(oxidizer)
iDilu = gas.species_index(diluent)
MWF, MWO = gas.molecular_weights[iFuel], gas.molecular_weights[iOxid]
MWD = gas.molecular_weights[iDilu]

nuOF = gas.stoich_air_fuel_ratio(fuel, oxidizer, basis='mole')  # Mass based
stoichOF = nuOF * MWF / MWO  # stoichOF = 0.5  # Mole based

X = np.zeros(gas.n_species)
X[iFuel] = 1.
X[iOxid] = stoichOF
diluent_fuel_ratio = np.linspace(0, 50, 500)

print('qDF, LeF, LeO, phou, rhob, alpha, Tf')
for qDF in diluent_fuel_ratio:
    X[iDilu] = qDF
    gas.TPX = Tin, Pin, X
    gas.transport_model = 'Mix'
    LeX = gas.thermal_conductivity / (gas.density_mass * gas.cp_mass *
                                      gas.mix_diff_coeffs)
    rhou = gas.density
    gas.equilibrate('HP')
    rhob = gas.density
    print(
        '{0:.3f}, {1:.4f}, {2:.4f}, {3:.4f}, {4:.4f}, {5:.4f}, {6:.4f}'.format(
            qDF, LeX[iFuel], LeX[iOxid], rhou, rhob, rhou / rhob, gas.T, qDF))
