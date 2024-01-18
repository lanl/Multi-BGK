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

Lx = 2.574                 #size of domain in cm
Nx = 992
nspec = 1

infile = 'scott_gap_231204.txt'
outfile = 'scott_array_cutout_500.init'  #Name of file to be written

# Data per species
Temps = [8.61e-5]       #eV

#!-----------------------------!#


#Create the grid on -Lx/2 : Lx/2
dx = float(Lx) / float(Nx)
cellcenters = np.linspace(-0.5*float(Lx) + 0.5*dx, 0.5*float(Lx) - 0.5*dx, Nx)

# Create the file
fout = open(outfile, 'w') 

densvals_file = np.loadtxt(infile,delimiter=',')
print(densvals_file.shape)
shift = int(Nx / 4)
dens = np.zeros(Nx)

for species in range(0,nspec):
    fout.write("species\n")
    fout.write(str(species) + "\n")
    
    fout.write("dens\n")
    for i in range(Nx):
        # Calculate Maxwellian
        if(i < shift):
            dens[i] = densvals_file[0][2]
        elif(i > Nx - shift - 1):
            dens[i] = densvals_file[-1][2]
        else:
            dens[i] = densvals_file[i - shift][2]         

    for i in range(Nx):
        fout.write(str(dens[i]) + "\n")            
    
    fout.write("velocity\n")
    for x in range(Nx):
        # Calculate Maxwellian        
        fout.write('0.0' + "\n")
    
    fout.write("temperature\n")
    for x in range(Nx):
        # Calculate Maxwellian        
        fout.write(str(Temps[species]) + "\n")
    
fout.write("Stop\n")
fout.close()
