import cantera as ct
import numpy as np
import matplotlib.pylab as plt 

gas1 = ct.Solution('GRI30.xml','gri30_mix')
Tin=300
Pin=1*101325
phi=1

ifuel = gas1.species_index('CH4')
io2 = gas1.species_index('O2')
in2 = gas1.species_index('N2')

air_N2_O2_ratio = 3.76
stoich_O2 = 2

X = np.zeros(gas1.n_species)
X[ifuel] = phi
X[io2] = stoich_O2
X[in2] = stoich_O2*air_N2_O2_ratio

gas1.TPX = Tin, Pin, X


# Refined grid at inlet and outlet, 6 points in x-direction :
initial_grid = 5*np.array([0.0, 0.001, 0.01, 0.02, 0.029, 0.03],'d')/3 # m
# Uniform grid, 6 points in x-direction (import numpy):
#initial_grid = 0.02*array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0],'d') # m		
# Uniform grid of 300 points using numpy :
#initial_grid = numpy.linspace(0,0.02 , 300)

#Set tolerance properties
tol_ss    = [1.0e-5, 1.0e-10]        # [rtol atol] for steady-state problem
tol_ts    = [1.0e-5, 1.0e-10]        # [rtol atol] for time stepping

loglevel  = 1                       # amount of diagnostic output (0 to 5)	    
refine_grid = True                  # True to enable refinement, False to disable 				   

# Creation de l'objet Flamme laminaire 1D
f = ct.FreeFlame(gas1, initial_grid)

f.flame.set_steady_tolerances(default=tol_ss)
f.flame.set_transient_tolerances(default=tol_ts)

f.inlet.X = X
f.inlet.T = Tin

# Première itération
# Pas d'équation d'énergie. A faire si schéma cinétique un peu trop costaud
f.energy_enabled = True
#Max number of times the Jacobian will be used before it must be re-evaluated
f.set_max_jac_age(50, 50)
#Set time steps whenever Newton convergence fails
f.set_time_step(0.1e-06, [2, 5,10, 20, 80]) #s

f.set_refine_criteria(ratio = 3.0, slope = 0.02, curve = 0.02, prune = 0.01)
f.solve(loglevel, refine_grid)

rnet = f.net_production_rates   # Unités Kmol/m3/s
hnet = f.heat_release_rate      # Unités W/m3

print ('Laminar flame speed = ',f.u[0],'m/s')

rnet = f.net_production_rates   # Unités Kmol/m3/s
hnet = f.heat_release_rate      # Unités W/m3

R11=gas1.reaction(11).equation
print(R11)
h_reaction_net_R11= - f.heat_production_rates[11] 

z = f.flame.grid
T = f.T
u = f.u

fig=plt.figure(1)

# create first subplot - adiabatic flame temperature
a=fig.add_subplot(221)
a.plot(z,T)
plt.title(r'T vs. Position')
plt.xlabel(r'Position [m]', fontsize=12)
plt.ylabel(r'Temperature [K]',fontsize=12)
a.xaxis.set_major_locator(plt.MaxNLocator(4)) # this controls the number of tick marks on the axis

# create second subplot - velocity
b=fig.add_subplot(222)
b.plot(z,u, label='$V$')
plt.title(r'Velocity vs. Position')
plt.xlabel(r'Position [m]', fontsize=12)
plt.ylabel(r'Velocity [m/s]',fontsize=12)
b.xaxis.set_major_locator(plt.MaxNLocator(4)) 

# create third subplot - rho
c=fig.add_subplot(223)
p = np.zeros(f.flame.n_points,'d')
for n in range(f.flame.n_points):
    f.set_gas_state(n)
    p[n]= gas1.density_mass
c.plot(z,p)
plt.title(r'Rho vs. Position')
plt.xlabel(r'Position [m]', fontsize=12)
plt.ylabel(r'Rho [kg/m^3]',fontsize=12)
c.xaxis.set_major_locator(plt.MaxNLocator(4)) 


# create fourth subplot - specie CH4
d=fig.add_subplot(224)
ch4 = np.zeros(f.flame.n_points,'d')
for n in range(f.flame.n_points):
    f.set_gas_state(n)
    ch4[n]= gas1.Y[ifuel]
d.plot(z,ch4)
plt.title(r'CH4 vs. Position')
plt.xlabel(r'Position [m]', fontsize=12)
plt.ylabel(r'CH4 Mole Fraction',fontsize=12)
d.xaxis.set_major_locator(plt.MaxNLocator(4))

#plt.subplots_adjust(left=0, right=1, wspace=0.5,top=0.9)
plt.tight_layout()
plt.savefig("test.png")
plt.show()

f.show_stats


