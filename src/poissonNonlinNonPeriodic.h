/*
 This solves the nonlinear Poisson equation
 -D_xx phi + g(phi) = s
 where g is a nonlinear function and s is a given source 

 The discretization is
 A_ij phi_j + g(phi_i) = s_i
 where A_ij is the usual discrete Laplcian on the 3-point stencil
 
*/
#include "gsl/gsl_linalg.h"


/**************************************************************/
/* PoissNonlinPeriodic1D                                      */
/*                                                            */
/* INPUTS                                                     */
/* N: number of cells in domain                               */
/* source: RHS vector, length N                               */
/* dx: size of each cell (fixed)                              */
/* Te: Electron temperature (eV)                              */
/*                                                            */
/* OUTPUT                                                     */
/* phi: solution to the equation                              */
/**************************************************************/

void PoissNonlinNonPeriodic1D(int N, int *order, double *source, double dx, double Lx, double *phi, double *Te);

void PoissLinNonPeriodic1D(int N, int *order, double *source, double dx, double Lx, double *phi, double *Te);

void PoissLinNonPeriodic1D_TF(int N, int *order, double *source, double dx, double Lx, double *phi, double *Te);

void PoissNonlinNonPeriodic1D_TF(int N, int *order, double *source, double dx, double Lx, double *phi, double *Te);

void simpleNonPeriodicPoisson(int N, double *source, double dx, double Lx, double *phi);
