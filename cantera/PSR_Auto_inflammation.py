
import sys
from cantera import *
import numpy as np

cti = Solution('gri30.xml')
air = Solution('air.xml')


m=cti.n_species
fuel_species = 'CH4'
ifuel = cti.species_index(fuel_species)
io2 = cti.species_index('O2')
in2 = cti.species_index('N2')

if ifuel < 0:
    raise "fuel species "+fuel_species+" not present!"

if cti.n_atoms(fuel_species,'O') > 0 or  cti.n_atoms(fuel_species,'N') > 0:
    raise "Error: only hydrocarbon fuels are supported."

phi = input('Enter Stoichiometric ratio phi : ')
phi = float(phi)


#Air composition
air_N2_O2_ratio = 3.76   
stoich_O2 = cti.n_atoms(fuel_species,'C') + 0.25*cti.n_atoms(fuel_species,'H')

#Mass fraction vector
x = np.zeros(m,'d')
x[ifuel] = phi
x[io2] = stoich_O2
x[in2] = stoich_O2*air_N2_O2_ratio

# Specify intial pressures and temperature of the reactor
Ti = input('Enter temperature (in kelvin) : ')
Ti = float(Ti)       # Kelvin

Pi = input('Enter pressure (in bar) : ')
Pi = float(Pi)*1e5         # Pascal


#Set initial conditions
cti.TPX = Ti, Pi, x


#################################################################
# Program starts here
#################################################################
#Create the batch reactor
r = IdealGasReactor(cti)

##Specify the conditions: Pression or Volume constant
print "--------------------------------------------------- "
print "    Equilibirum conditions: "
print "--------------------------------------------------- "
print ""
print "For a constant volume equilibrium, enter :      UV "
print "For a constant pressure equilibrium, enter :    HP "
print ""
cond  = raw_input('Specify the equilibrium condition : ')
cond  = str(cond)

while cond != 'HP' and cond != 'UV':
    print ("You must choose between UV and HP !  ")
    cond  = raw_input('Specify the equilibrium condition : ')
    cond  = str(cond)

#Particular case of a constant-pressure reactor
if cond == 'HP':
	# Define a wall between the reactor and the environment, and
	# make it flexible, so that the pressure in the reactor is held
	# at the environment pressure.
    env = Reservoir(air)
    w = Wall(r,env)
    w.expansion_rate_coeff = 1.0e6 # set expansion parameter. dV/dt = KA(P_1 - P_2)
    w.area = 1.0 # set wall area


# Now create a reactor network consisting of the single batch reactor
# Reason: the only way to advance reactors in time is through a network
sim = ReactorNet([r])

#################
#Computational properties: we're going to advance the network in time
# Initial simulation time
time = 0.0e-0

# Specify the number of time steps
nt = input('Enter number of time steps: ')
nt = int(nt)

# Specify the time step
dtms = input('Enter the time step (in ms): ')
dtms = int(dtms)
dt = dtms * 1.0e-3 #s

#################
#Run the simulation

tim = np.zeros(nt,'d')
temp = np.zeros(nt,'d')
press = np.zeros(nt,'d')
mfrac = np.zeros([nt,m],'d')

#Loop for nt time steps of dt seconds. 
print ('time [s] ,   T [K] ,   p [Pa] ,   u [J/kg]')
for n in range(nt):
    time += dt
    sim.advance(time)
    tim[n] = time
    temp[n] = r.T
    press[n] = r.thermo.P
    for i in range(m):
        mfrac[n,i]=r.thermo[cti.species_name(i)].Y
    print ('%10.3e %10.3f %10.3f %14.6e' % (sim.time, r.T, 
                                           r.thermo.P, r.thermo.h))


##################################################################
## Save your results if needed
##################################################################
## write output CSV file for importing into Excel
#if cond == 'HP':
#     csvfile = 'Reactor_HP.csv'
#elif cond == 'UV':
#     csvfile = 'Reactor_UV.csv'
#
#csv_file = csvfile
#with open(csv_file, 'w') as outfile:
#	writer = csv.writer(outfile)
#	writer.writerow(['Time','Temperature','Pressure']+cti.species_names)
#	for n in range(nt):
#    		writer.writerow([tim[n], temp[n], press[n]]+list(mfrac[n,:]))
#print 'output written to '+csvfile


