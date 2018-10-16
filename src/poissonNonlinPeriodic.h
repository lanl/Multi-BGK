/*
 This solves the nonlinear Poisson equation
 -D_xx phi + g(phi) = s
 where g is a nonlinear function and s is a given source

 The discretization is
 A_ij phi_j + g(phi_i) = s_i
 where A_ij is the usual discrete Laplcian on the 3-point stencil
*/

#include <gsl/gsl_linalg.h>

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

void get_uniform_Te(double *Te, int Nx, double T0);

void get_ramp_Te(double *Te, int Nx, double T_start, double T_end, double t,
                 double tfinal);

void PoissNonlinPeriodic1D(int N, double *source, double dx, double Lx,
                           double *phi, double *Te);

void PoissLinPeriodic1D(int N, double *source, double dx, double Lx,
                        double *phi, double *Te);

void PoissLinPeriodic1D_TF(int N, double *source, double dx, double Lx,
                           double *phi, double *Te);

void PoissNonlinPeriodic1D_TF(int N, double *source, double dx, double Lx,
                              double *phi, double *Te);

void simplePoisson(int N, double *source, double dx, double Lx, double *phi);
