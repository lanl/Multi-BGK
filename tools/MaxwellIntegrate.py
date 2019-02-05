"""Script that calculates numerical error in calculating density and temperature for a Maxwellian 
using a uniform velocity grid

See __main__ below for setting up the species

"""



import numpy as np
import math

#boltzmann constant

#Use this if temperature units are K
#kb = 1.38064853e-23

#Use this if temperature units are ev
kb = 1


def buildMaxwellian(m,n,T,c2):
    """This function builds a Maxwellian with the given velocity grid c, mass m, 
    number density n, and temperature T. The bulk velocity is assumed to be zero.

    inputs [units]
    -------------
    m: mass [kg]
    n: number density [1/m^3]
    T: Temperature [K]
    c2: 3d array of velocities squared, flattened to 1D [m/s]

    output:
    -------

    M: flattened 3D array of distribution function values [1/m^3 / (m/s)^3]
    """
    
    Nv3 = len(c2)

    M = np.zeros(Nv3)

    kT = kb*T

    prefac = n * (m / (2.0*np.pi*kT))**1.5

    for i in range(Nv3):
        M[i] = prefac * math.exp(-0.5*m*c2[i]/kT)

    return M


def getVeloGrid(m,T,Nv,widthfactor):
    """This function returns 1d and flattened 3d velocity grids

    inputs [units]:
    --------------
    
    m: mass [kg]
    T: temperature [K]
    Nv: Number of velocity grid points in oen direction
    widthfactor: multiples of thermal speed for the width

    outputs:
    -------
    c: array of 1d velocities [m/s]
    (c3x, c3y,c3z): flattened 3d velocity grid [m/s]
    """

    kT = kb * T                    #J
    therm_speed = np.sqrt(kT/m)    #m/s
    Lv = widthfactor * therm_speed
    dv = 2.0*Lv/(Nv - 1.0)
    c = -Lv*np.ones(Nv) + dv*np.arange(Nv)

    c3 = np.zeros((3,Nv*Nv*Nv))

    for i in range(Nv):
        for j in range(Nv):
            for k in range(Nv):
                index = k + Nv*(j + Nv*i)
                c3[0,index] = c[i]
                c3[1,index] = c[j]
                c3[2,index] = c[k]

    return c, c3

def getc2(c3):
    """calculates c^2 at every point on the velocity grid"""

    Nv3 = len(c3[0])
    c2 = np.zeros(Nv3)

    for i in range(Nv3):
        c2[i] = c3[0,i]**2 + c3[1,i]**2 + c3[2,i]**2

    return c2

def getwts(c):
    """returns the flattened 3d array of integration weights for the uniform velocities c"""

    Nv = len(c)

    wts = np.zeros(Nv)
    wts3 = np.zeros(Nv*Nv*Nv)

    dv = c[1] - c[0]

    for i in range(Nv):
        wts[i] = dv
    wts[0] *= 0.5
    wts[-1] *= 0.5

    for i in range(Nv):
        for j in range(Nv):
            for k in range(Nv):
                index = k + Nv*(j + Nv*i)
                wts3[index] = wts[i]*wts[j]*wts[k]

    return wts3

if __name__ == "__main__":

    #THIS IS THE LIGHT SPECIES THAT I USED TO DEFINE THE GRID

    amu = 1.6605e-27
    m1 = 1.0*amu
    m2 = 12.0*amu
    n = 3.241e21                  #initial number density 1/cm^3
    T = 537.0                 
    T_low = 10.0

    Nv = 64                #points per direction in velocity space
    width = 8              #Multiple of thermal speeds to define edge of velocity space

    c, c3 = getVeloGrid(m1,T,Nv,width)

    c2 = getc2(c3)
    wts = getwts(c)

    dc = c[1] - c[0]
                   

    LightMaxwellian = buildMaxwellian(m1,n,T,c2)
    print "Max speed:", c[-1]

    #First, check to see how well different temperatures can be represented for species 1 on this grid

    print "\n\nOriginal Maxwellian that the grid is based on\n\n"

    n_approx = sum(wts*LightMaxwellian)
    T_approx = sum(wts*c2*LightMaxwellian)*m1/3.0/n/kb        

    # print n, n_approx, 100.*math.fabs(n_approx-n)/n

    print "Actual high temperature, light species: ", T
    print "Approximate high temperature, light species", T_approx
    print "Relative error", 100.*math.fabs(T_approx-T)/T


    print "\n\nSame grid, low temperature for species 1\n\n"

    LightMaxwellianLowTemp = buildMaxwellian(m1,n,T_low,c2)

    n_approx = sum(wts*LightMaxwellianLowTemp)
    T_approx = sum(wts*c2*LightMaxwellianLowTemp)*m1/3.0/n/kb

    #print n, n_approx, 100*math.fabs(n_approx-n)/n

    print "Actual high temperature, light species: ", T_low
    print "Approximate high temperature, light species", T_approx
    print "Relative error", 100.*math.fabs(T_approx-T_low)/T_low

    
    #Next, check to see how well the grid defined by species 1 can represent a distribution function of species 2

    print "\n \n Same grid, Heavy species \n\n"

    HeavyMaxwellian = buildMaxwellian(m2,n,T,c2)

    n_approx = sum(wts*HeavyMaxwellian)
    T_approx = sum(wts*c2*HeavyMaxwellian)*m2/3.0/n/kb

    #print n, n_approx, 100.*math.fabs(n_approx-n)/n
 
    print "Actual high temperature, heavy species: ", T
    print "Approximate high temperature, heavy species", T_approx
    print "Relative error", 100.*math.fabs(T_approx-T)/T


    print "\n \n Same grid, heavy species at low temperature \n\n"

    HeavyMaxwellian_low = buildMaxwellian(m2,n,T_low,c2)

    n_approx = sum(wts*HeavyMaxwellian_low)
    T_approx = sum(wts*c2*HeavyMaxwellian_low)*m2/3.0/n/kb

    #print n, n_approx, 100.*math.fabs(n_approx-n)/n

    print "Actual high temperature, heavy species: ", T_low
    print "Approximate high temperature, heavy species", T_approx
    print "Relative error", 100.*math.fabs(T_approx-T_low)/T_low


    beta1 = dc / np.sqrt(2.0*kb*T/m1)
    beta2 = dc / np.sqrt(2.0*kb*T/m2)
    beta1_low = dc / np.sqrt(2.0*kb*T_low/m1)
    beta2_low = dc / np.sqrt(2.0*kb*T_low/m2)

    print "\n\n \Delta v divided by thermal speeds for these four cases"
    print beta1, beta2, beta1_low, beta2_low, "\n"
