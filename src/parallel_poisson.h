#include "mpi.h"

#ifndef _PARALLEL_POISSON_H
#define _PARALLEL_POISSON_H

//parallel poisson solver that gives testament to the need to use structs to group variables together..
void poisson_solver(MPI_Comm *comm, int *rank, int *numRanks, MPI_Status *status, int *Nx, int **Nx_ranks, int *Nx_rank, 
                    int *ecouple, int *bcs, int *poissFlavor, int *ionFix, int *nspec, 
                    int **Z_max, double *Te_ref,  double *Te_start, int *order, double *dx, double *Lx, double *t, double *tfinal, double **Te_arr, 
                    double **Te_arr_allranks, double **T0_oned, double ***n_oned, double ***Z_oned, 
                    double **source, double **source_buf, double **source_allranks, double **PoisPot, double **PoisPot_allranks,
                    double **T0_bcs, double ***n_bcs, double ***Z_bcs, double **Te_bcs);

#endif //_PARALLEL_POISSON_H