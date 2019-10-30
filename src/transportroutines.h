//Finite Volume code for LHS of 1D-3V Boltzmann

//Data structure for species information
//#include "species.h"

//Sets up everything used in transport side inputs are:
//  numV - Number of velocity points in each coordinate direction
//  numX - Number of physical points in x-direction. Note - if doing this with MPI, this is the number of points that are handled by the current process. Periodic + MPI not currently implemented
//  lv - Semi-length of velocity domain (v is assumed to be in [-lv,lv]
//  xnodes - ordered array of locations of the center of each node
//  dxnodes - ordered array with size of each node
//  Lx - width of x domain
//  vel - array of velocity grid points for one direction (full 3d velocity is just this tensored)
//  IC - flag that tells type of problem being solved (and thus, boundary condition) - currently just set to 6 for periodic BCs, other options from before are removed
//  TWall_in - specifies temperature of walls for certian BCs, can be ignored otherwise
void initialize_transport(int bcs, double ****f, int numV, int numX, int nspec, double *xnodes, double *dxnodes, double Lx, double **vel, int ord, double timestep);

/*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/
//these are what you actually call from the outside
//f is the input distribution - first coordinate is the x location, second is species, third is the v location
//f_conv is the output
//
//f array sizes: f[numX][nspec][numV*numV*numV]

//id is the process id, so it knows if it's at a boundary or interior of domain for MPI communication

void advectOne(double ***f, double *PoisPot, double **qm, double m, int sp);

void advectTwo(double ***f, double *PoisPot, double **qm, double m, int sp); 

void advectTwo_x(double ***f, double ***f_conv, int sp); 

void advectTwo_v(double ***f, double ***f_conv, double *PoisPot, double **qm, double m, int sp); 


//void upwindTwo(double ***f, double ***f_conv, double *PoisPot, double **qm, double m, int sp); 

/*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/

//deallocates stuff allocated in the initializer
void dealloc_trans();
