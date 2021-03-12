""" 
Script for creating an initial condition to load into Multi-BGK
edit stuff below

Generates Gaussian with mean 0 and variance sigma (in cm) 

f(x) = \frac{1}{\sigma(2\pi)^{1/2}} e^{-x^2 / 2 sigma^2}

input file layout  
Do it by species, e.g. 
species
0
n
1.0 1.0 1.0 ... 1.0
velocity
1.0 1.0 1.0 ... 1.0
temperature
1.0 1.0 1.0 ... 1.0

"""

import numpy as np

#!-------- EDIT THESE ---------!#

Lx = 2.0                 #size of domain in cm
Nx = 64
nspec = 2

outfile = 'Gaussian.init'  #Name of file to be written


# Data per species
dens_max = [1.0e16, 1.0e16]  #1/cm^3
Temps = [100.0, 100.0]       #eV
std_dev = [0.2, 0.2]

#!-----------------------------!#


#Create the grid on -Lx/2 : Lx/2
dx = float(Lx) / float(Nx)
cellcenters = np.linspace(-0.5*float(Lx) + 0.5*dx, 0.5*float(Lx) - 0.5*dx, Nx)

# Create the file
fout = open(outfile, 'w') 
    
for species in np.arange(0,nspec):
    fout.write("species\n")
    fout.write(str(species) + "\n")
    
    fout.write("dens\n")
    for x in cellcenters:
        # Calculate Maxwellian        
        dens = dens_max[species]*(1.0 / std_dev[species]) * (2.0 * np.pi) ** -0.5 * np.exp(-0.5 * (x / std_dev[species]) ** 2.0)
        fout.write(str(dens) + " ")
    fout.write("\n")
    
    fout.write("velocity\n")
    for x in cellcenters:
        # Calculate Maxwellian        
        fout.write('0.0' + " ")
    fout.write("\n")
    
    fout.write("temperature\n")
    for x in cellcenters:
        # Calculate Maxwellian        
        fout.write(str(Temps[species]) + " ")
    fout.write("\n")
    
