#include "transportroutines.h"
#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

// Computes the transport term
static int BC;
static int N;
static int nX;
static int ns;
static int order;
static double Lx;
static double h_v;
static double dt;
static double *x, *dx, **c;
static double ***f_star;

/*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/
double min(double in1, double in2) { return in1 < in2 ? in1 : in2; }

/*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/
double max(double in1, double in2) { return in1 > in2 ? in1 : in2; }

/*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/
double minmod(double in1, double in2, double in3) {
  if ((in1 > 0) && (in2 > 0) && (in3 > 0)) {
    return min(min(in1, in2), in3);
  } else if ((in1 < 0) && (in2 < 0) && (in3 < 0)) {
    return max(max(in1, in2), in3);
  } else
    return 0;
}

void initialize_transport(int bcs, int numV, int numX, int nspec, double *xnodes,
                          double *dxnodes, double Lx_val, double **vel, int ord,
                          double timestep) {
  BC = bcs;
  N = numV;
  nX = numX; // Remember that this is nX in the rank
  x = xnodes;
  dx = dxnodes;
  ns = nspec;
  Lx = Lx_val;

  c = vel;

  order = ord;
  dt = timestep;

  int i, l;
  f_star = (double ***)malloc((nX + 2 * order) * sizeof(double **));
  for (l = 0; l < nX + 2 * order; l++) {
    f_star[l] = (double **)malloc(nspec * sizeof(double *));
    for (i = 0; i < nspec; i++)
      f_star[l][i] = malloc(N * N * N * sizeof(double));
  }

  x = xnodes;
  dx = dxnodes;
}

/*$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/
void fillGhostCells_firstorder(double ***f, int sp) {

  int rank, numRanks;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
  MPI_Status status;

  int i;

  int left_ghost = 0;
  int left_actual = 1;
  int right_actual = nX;
  int right_ghost = nX + 1;
  int up_nbr, down_nbr;

  double xlocal[N * N * N];
  double xremote[N * N * N];

  // If just one rank, no MPI communication is needed
  if (numRanks == 1) {
    if(BC == 0){
      for (i = 0; i < N * N * N; i++) {
        f[left_ghost][sp][i] = f[right_actual][sp][i];
        f[right_ghost][sp][i] = f[left_actual][sp][i];
      }
    }
  } else {

    /* Processors 0 and 1 exchange, 2 and 3 exchange, etc.  Then
       1 and 2 exchange, 3 and 4, etc.  The formula for this is
       if (even) exchng up else down
       if (odd)  exchng up else down
    */
    /* Note the use of xlocal[i] for &xlocal[i][0] */
    /* An edge case is when there are an odd number of ranks - then the last
     * (even) node is sending to 0, which is also an even node. */

    up_nbr = rank + 1;
    if (up_nbr >= numRanks)
      up_nbr = 0;
    down_nbr = rank - 1;
    if (down_nbr < 0)
      down_nbr = numRanks - 1;

    if ((rank % 2) == 0) {

      if (rank != numRanks - 1) {
        for (i = 0; i < N * N * N; i++){
          xlocal[i] = f[right_actual][sp][i];
        }

        /* exchange up */
        MPI_Send(&xlocal[0], N * N * N, MPI_DOUBLE, up_nbr, 1, MPI_COMM_WORLD);

        MPI_Recv(&xremote[0], N * N * N, MPI_DOUBLE, up_nbr, 2, MPI_COMM_WORLD,
                 &status);

        for (i = 0; i < N * N * N; i++){
          f[right_ghost][sp][i] = xremote[i];
        }
      }
    } else {

      for (i = 0; i < N * N * N; i++)
        xlocal[i] = f[left_actual][sp][i];

      /* exchange down */
      MPI_Recv(&xremote[0], N * N * N, MPI_DOUBLE, down_nbr, 1, MPI_COMM_WORLD,
               &status);

      MPI_Send(&xlocal[0], N * N * N, MPI_DOUBLE, down_nbr, 2, MPI_COMM_WORLD);

      for (i = 0; i < N * N * N; i++)
        f[left_ghost][sp][i] = xremote[i];
    }

    // Deal with edge case - boundary for odd numRanks case
    if (numRanks % 2 == 1) {
      if(BC == 0){
        if (rank == 0) {
          for (i = 0; i < N * N * N; i++)
            xlocal[i] = f[left_actual][sp][i];

          MPI_Recv(&xremote[0], N * N * N, MPI_DOUBLE, numRanks - 1, 1,
                   MPI_COMM_WORLD, &status);

          MPI_Send(&xlocal[0], N * N * N, MPI_DOUBLE, numRanks - 1, 2,
                   MPI_COMM_WORLD);
          for (i = 0; i < N * N * N; i++)
            f[left_ghost][sp][i] = xremote[i];
        } else if (rank == numRanks - 1) {

          for (i = 0; i < N * N * N; i++)
            xlocal[i] = f[right_actual][sp][i];

          MPI_Send(&xlocal[0], N * N * N, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);

          MPI_Recv(&xremote[0], N * N * N, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD,
                   &status);

          for (i = 0; i < N * N * N; i++)
            f[right_ghost][sp][i] = xremote[i];
        }
      }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Do the second set of exchanges
    // In this case, there will always be an even rank to the right (because 0
    // is even)
    // However we could again have the case where rank 0 and numRanks-1 are both
    // even. This is taken care of above.
    if ((rank % 2) == 1) {

      for (i = 0; i < N * N * N; i++)
        xlocal[i] = f[right_actual][sp][i];

      /* exchange up */
      MPI_Send(&xlocal[0], N * N * N, MPI_DOUBLE, up_nbr, 3, MPI_COMM_WORLD);

      MPI_Recv(&xremote[0], N * N * N, MPI_DOUBLE, up_nbr, 4, MPI_COMM_WORLD,
               &status);

      if(BC == 0 || rank != numRanks-1){
        for (i = 0; i < N * N * N; i++){
          f[right_ghost][sp][i] = xremote[i];
        }
      }
    } else {

      if (((numRanks % 2 == 0) ||
            (rank != 0))) { // this case was dealt with above
        for (i = 0; i < N * N * N; i++)
          xlocal[i] = f[left_actual][sp][i];

        /* exchange down */
        MPI_Recv(&xremote[0], N * N * N, MPI_DOUBLE, down_nbr, 3,
                 MPI_COMM_WORLD, &status);

        MPI_Send(&xlocal[0], N * N * N, MPI_DOUBLE, down_nbr, 4,
                 MPI_COMM_WORLD);

        if(BC == 0 || rank != 0){
          for (i = 0; i < N * N * N; i++){
            f[left_ghost][sp][i] = xremote[i];
          }
        }
      }
    }

    // Ghost cells complete, just wait for all ranks to catch up
    MPI_Barrier(MPI_COMM_WORLD);
  }
}

void fillGhostCells_secondorder(double ***f, int sp) {

  int rank, numRanks;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
  MPI_Status status;
  int i;

  // Indentify ghost and shared cells
  // Note that '1' and '2' refer to distance from the processor boundary
  int left_ghost_2 = 0;
  int left_ghost_1 = 1;
  int left_actual_1 = 2;
  int left_actual_2 = 3;
  int right_actual_2 = nX;
  int right_actual_1 = nX + 1;
  int right_ghost_1 = nX + 2;
  int right_ghost_2 = nX + 3;
  int up_nbr, down_nbr;

  int N3 = N * N * N;

  double xlocal[2 * N3];
  double xremote[2 * N3];

  // If just one rank, no MPI communication is needed
  if (numRanks == 1) {
    if(BC==0){
      for (i = 0; i < N3; i++) {
        f[left_ghost_1][sp][i] = f[right_actual_1][sp][i];
        f[left_ghost_2][sp][i] = f[right_actual_2][sp][i];
        f[right_ghost_1][sp][i] = f[left_actual_1][sp][i];
        f[right_ghost_2][sp][i] = f[left_actual_2][sp][i];
      }
    }
  } else {

    /* Processors 0 and 1 exchange, 2 and 3 exchange, etc.  Then
       1 and 2 exchange, 3 and 4, etc.  The formula for this is
       if (even) exchng up else down
       if (odd)  exchng up else down
    */
    /* Note the use of xlocal[i] for &xlocal[i][0] */
    /* An edge case is when there are an odd number of ranks - then the last
     * (even) node is sending to 0, which is also an even node. */

    up_nbr = rank + 1;
    if (up_nbr >= numRanks)
      up_nbr = 0;
    down_nbr = rank - 1;
    if (down_nbr < 0)
      down_nbr = numRanks - 1;

    if ((rank % 2) == 0) {

      if (rank != numRanks - 1) {
        for (i = 0; i < N3; i++) {
          xlocal[i] = f[right_actual_1][sp][i];
          xlocal[i + N3] = f[right_actual_2][sp][i];
        }

        /* exchange up */
        MPI_Send(&xlocal[0], 2 * N3, MPI_DOUBLE, up_nbr, 1, MPI_COMM_WORLD);

        MPI_Recv(&xremote[0], 2 * N3, MPI_DOUBLE, up_nbr, 2, MPI_COMM_WORLD,
                 &status);

        for (i = 0; i < N3; i++) {
          f[right_ghost_1][sp][i] = xremote[i];
          f[right_ghost_2][sp][i] = xremote[i + N3];
        }
      }
    } else {

      for (i = 0; i < N3; i++) {
        xlocal[i] = f[left_actual_1][sp][i];
        xlocal[i + N3] = f[left_actual_2][sp][i];
      }

      /* exchange down */
      MPI_Recv(&xremote[0], 2 * N3, MPI_DOUBLE, down_nbr, 1, MPI_COMM_WORLD,
               &status);

      MPI_Send(&xlocal[0], 2 * N3, MPI_DOUBLE, down_nbr, 2, MPI_COMM_WORLD);

      for (i = 0; i < N3; i++) {
        f[left_ghost_1][sp][i] = xremote[i];
        f[left_ghost_2][sp][i] = xremote[i + N3];
      }
    }

    // Deal with edge case - boundary for odd numRanks case
    if (numRanks % 2 == 1) {
      if(BC == 0){
        if (rank == 0) {
          for (i = 0; i < N3; i++) {
            xlocal[i] = f[left_actual_1][sp][i];
            xlocal[i + N3] = f[left_actual_2][sp][i];
          }

          MPI_Recv(&xremote[0], 2 * N3, MPI_DOUBLE, numRanks - 1, 1,
                   MPI_COMM_WORLD, &status);

          MPI_Send(&xlocal[0], 2 * N3, MPI_DOUBLE, numRanks - 1, 2,
                   MPI_COMM_WORLD);

          for (i = 0; i < N3; i++) {
            f[left_ghost_1][sp][i] = xremote[i];
            f[left_ghost_2][sp][i] = xremote[i + N3];
          }
        } else if (rank == numRanks - 1) {

          for (i = 0; i < N3; i++) {
            xlocal[i] = f[right_actual_1][sp][i];
            xlocal[i + N3] = f[right_actual_2][sp][i];
          }

          MPI_Send(&xlocal[0], 2 * N3, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);

          MPI_Recv(&xremote[0], 2 * N3, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD,
                   &status);

          for (i = 0; i < N3; i++) {
            f[right_ghost_1][sp][i] = xremote[i];
            f[right_ghost_2][sp][i] = xremote[i + N3];
          }
        }
      }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Do the second set of exchanges
    // In this case, there will always be an even rank to the right (because 0
    // is even)
    // However we could again have the case where rank 0 and numRanks-1 are both
    // even. This is taken care of above.
    if ((rank % 2) == 1) {

      for (i = 0; i < N3; i++) {
        xlocal[i] = f[right_actual_1][sp][i];
        xlocal[i + N3] = f[right_actual_2][sp][i];
      }

      /* exchange up */
      MPI_Send(&xlocal[0], 2 * N3, MPI_DOUBLE, up_nbr, 3, MPI_COMM_WORLD);

      MPI_Recv(&xremote[0], 2 * N3, MPI_DOUBLE, up_nbr, 4, MPI_COMM_WORLD,
               &status);

      if(BC == 0 || rank != numRanks-1){
        for (i = 0; i < N3; i++) {
          f[right_ghost_1][sp][i] = xremote[i];
          f[right_ghost_2][sp][i] = xremote[i + N3];
        }
      }
    } else {

      if (!((numRanks % 2 == 1) &&
            (rank == 0))) { // this case was dealt with above
        for (i = 0; i < N3; i++) {
          xlocal[i] = f[left_actual_1][sp][i];
          xlocal[i + N3] = f[left_actual_2][sp][i];
        }

        /* exchange down */
        MPI_Recv(&xremote[0], 2 * N3, MPI_DOUBLE, down_nbr, 3, MPI_COMM_WORLD,
                 &status);

        MPI_Send(&xlocal[0], 2 * N3, MPI_DOUBLE, down_nbr, 4, MPI_COMM_WORLD);

        if(BC == 0 || rank != 0){
          for (i = 0; i < N3; i++) {
            f[left_ghost_1][sp][i] = xremote[i];
            f[left_ghost_2][sp][i] = xremote[i + N3];
          }
        }
      }
    }

    // Ghost cells complete, just wait for all ranks to catch up
    MPI_Barrier(MPI_COMM_WORLD);
  }
}

// Computes the x direction first order upwind solution FOR A SINGLE SPECIES
void upwindOne_x(double ***f, double ***f_conv, double *v, int sp) {
  int i, j, k, l;
  int index;
  double CFL_NUM;

  // Update ghost with current solution, prior to transport
  fillGhostCells_firstorder(f, sp);

  // Note - the 'real' data lives in 1...Nx, the ghost points are 0 and Nx+1

  for (l = 1; l < nX + 1; l++) {
    //#pragma omp parallel for private(CFL_NUM,i,j,k,index)
    for (i = 0; i < N; i++) {
      CFL_NUM = dt * v[i] / dx[l];
      for (j = 0; j < N; j++)
        for (k = 0; k < N; k++) {
          index = k + N * (j + N * i);
          // the upwinding
          if (i < N / 2)
            f_conv[l][sp][index] =
                f[l][sp][index] -
                CFL_NUM * (f[l + 1][sp][index] - f[l][sp][index]);
          else
            f_conv[l][sp][index] =
                f[l][sp][index] -
                CFL_NUM * (f[l][sp][index] - f[l - 1][sp][index]);
        }
    }
  }

  // Get updated ghost solution
  fillGhostCells_firstorder(f_conv, sp);
}

// Computes the v direction first order upwind solution FOR A SINGLE SPECIES
void upwindOne_v(double ***f, double ***f_conv, double *PoisPot, double **qm,
                 double m, double *v, int sp) {
  int i, j, k, l;
  int index;
  double CFL_NUM, CFL_MAX = 0.0;
  h_v = v[1] - v[0];

  int rank, numNodes;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numNodes);
  // int numamt;

  // Set electric field acceleration term
  // This is given by
  // - Z (e_0 E) / m
  //
  // e_0 E = (e phi)_x is given in units of ergs
  // Z is the partial ionization, provided as qm in the input (no unit)
  // m is the mass of the species in grams
  //
  // Thus the acceleration term units are cm / s^2

  // Calculate the electric field from the potential.
  // Note that E is sized to nX+2 to ensure that l loop has enough memory
  // Also note that qm / Z is defined on nX, not nX+2
  double *E = malloc((nX + 2) * sizeof(double));
  for (l = 1; l < nX + 1; l++) {
    E[l] =
        -qm[l - 1][sp] / m * (PoisPot[l + 1] - PoisPot[l - 1]) / (2.0 * dx[l]);
  }

  for (l = 1; l < nX + 1; l++) {
    if (E[l] > 0) {
      CFL_NUM = dt * E[l] / h_v;

      if (CFL_NUM > CFL_MAX)
        CFL_MAX = CFL_NUM;

      //#pragma omp parallel for private(i,j,k,index)
      for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
          for (k = 0; k < N; k++) {
            index = k + N * (j + N * i);
            if (index >= N * N * N || index < 0) {
              printf("Index out of range");
              exit(1);
            }
            if (i != 0) {
              if (k + N * (j + N * (i - 1)) >= N * N * N ||
                  k + N * (j + N * (i - 1)) < 0) {
                printf("Index out of range");
                exit(1);
              }

              f_conv[l][sp][index] =
                  f[l][sp][index] -
                  CFL_NUM *
                      (f[l][sp][index] - f[l][sp][k + N * (j + N * (i - 1))]);
            } else {
              f_conv[l][sp][index] =
                  f[l][sp][index] - CFL_NUM * (f[l][sp][index]);
            }
          }
    } else {
      CFL_NUM = dt * E[l] / h_v;

      if (fabs(CFL_NUM) > CFL_MAX)
        CFL_MAX = fabs(CFL_NUM);

      //#pragma omp parallel for private(i,j,k,index)
      for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
          for (k = 0; k < N; k++) {
            index = k + N * (j + N * i);
            if (i != (N - 1))
              f_conv[l][sp][index] =
                  f[l][sp][index] -
                  CFL_NUM * (f[l][sp][k + N * (j + N * (i + 1))] -
                             f[l][sp][k + N * (j + N * i)]);
            else
              f_conv[l][sp][index] =
                  f[l][sp][index] + CFL_NUM * (f[l][sp][index]);
          }
    }
  }

  if (CFL_MAX > 1.0) {
    printf("\n***********\n");
    printf("Warning - Electric field CFL is %g\n", CFL_MAX);
    printf("Code may be unstable\n");
    printf("\n***********\n");
  }
  free(E);
}

// does a first order splitting, solving v first and returns the updated f
void advectOne(double ***f, double *PoisPot, double **qm, double m, int sp) {
  upwindOne_v(f, f_star, PoisPot, qm, m, c[sp], sp);
  upwindOne_x(f_star, f, c[sp], sp);
}

// Computes x portion of second order upwind solution, with minmod
void upwindTwo_x(double ***f, double ***f_conv, double *v, int sp) {
  int i, j, k, l;
  int index;
  double slope[3];
  double CFL_NUM;

  fillGhostCells_secondorder(f, sp);

  // main upwinding
  double f_l, f_r, f_ll, f_rr, x_l, x_r, x_ll, x_rr, dx_l, dx_r;

  for (l = 2; l < nX + 2; l++) {
    x_l = x[l - 1];
    x_ll = x[l - 2];
    dx_l = dx[l - 1];

    x_r = x[l + 1];
    x_rr = x[l + 2];
    dx_r = dx[l + 1];

    //#pragma omp parallel for private(i,j,k,index,CFL_NUM,f_l,f_ll,f_r,slope)

    for (i = N / 2; i < N; i++) {
      CFL_NUM = dt * v[i] / dx[l];

      for (j = 0; j < N; j++)
        for (k = 0; k < N; k++) {
          index = k + N * (j + N * i);
          // upwind coming from the left
          // generate the local slopes
          f_l = f[l - 1][sp][index];
          f_ll = f[l - 2][sp][index];

          f_r = f[l + 1][sp][index];

          slope[1] = minmod(2.0 * (f[l][sp][index] - f_l) / (x[l] - x_l),
                            2.0 * (f_r - f[l][sp][index]) / (x_r - x[l]),
                            (f_r - f_l) / (x_r - x_l));
          slope[0] = minmod(2.0 * (f_l - f_ll) / (x_l - x_ll),
                            2.0 * (f[l][sp][index] - f_l) / (x[l] - x_l),
                            (f[l][sp][index] - f_ll) / (x[l] - x_ll));

          f_conv[l][sp][index] =
              -CFL_NUM *
              (f[l][sp][index] + 0.5 * (dx[l] - v[i] * dt) * slope[1] -
               (f_l + 0.5 * (dx_l - v[i] * dt) * slope[0]));
        }
    }

    //#pragma omp parallel for private(i,j,k,index,CFL_NUM,f_l,f_rr,f_r,slope)
    for (i = 0; i < N / 2; i++) {
      CFL_NUM = dt * v[i] / dx[l];

      for (j = 0; j < N; j++)
        for (k = 0; k < N; k++) {
          index = k + N * (j + N * i);
          // upwind coming from the right
          // generate the local slopes
          f_l = f[l - 1][sp][index];

          f_r = f[l + 1][sp][index];
          f_rr = f[l + 2][sp][index];

          slope[2] = minmod(2.0 * (f_r - f[l][sp][index]) / (x_r - x[l]),
                            2.0 * (f_rr - f_r) / (x_rr - x_r),
                            (f_rr - f[l][sp][index]) / (x_rr - x[l]));

          slope[1] = minmod(2.0 * (f[l][sp][index] - f_l) / (x[l] - x_l),
                            2.0 * (f_r - f[l][sp][index]) / (x_r - x[l]),
                            (f_r - f_l) / (x_r - x_l));

          // the upwinding
          f_conv[l][sp][index] =
              -CFL_NUM *
              (f_r - 0.5 * (dx_r + v[i] * dt) * slope[2] -
               (f[l][sp][index] - 0.5 * (dx[l] + v[i] * dt) * slope[1]));
        }
    }
  }
}

// Computes second order upwind term for the field acceleration term
void upwindTwo_v(double ***f, double ***f_conv, double *PoisPot, double **qm,
                 double m, double *v, int sp) {
  int i, j, k, l;
  int index;
  double slope[3];
  double CFL_NUM, CFL_MAX = 0.0;

  h_v = v[1] - v[0];

  // Set electric field acceleration term
  // This is given by
  // - Z (e_0 E) / m
  //
  // e_0 E = (e phi)_x is given in units of ergs
  // Z is the partial ionization, provided as qm in the input (no unit)
  // m is the mass of the species in grams
  //
  // Thus the acceleration term units are cm / s^2

  // Note - dx and qm do not have ghost cells
  double *E = malloc((nX + 4) * sizeof(double));
  for (l = 2; l < nX + 2; l++) {
    E[l] = -qm[l - 2][sp] / m * (PoisPot[l + 1] - PoisPot[l - 1]) /
           (2.0 * dx[l - 2]);
  }

  if (fabs(v[2] - v[1] - h_v) > 1e-6) {
    printf("Error - you must use a uniform velocity grid for second order "
           "advection\n Please set discret=0 in your input file.\n");
    exit(1);
  }

  // main upwinding
  double up2, up1, neut, down1, down2;

  for (l = 2; l < nX + 2; l++) {
    CFL_NUM = dt * E[l] / h_v;

    if (fabs(CFL_NUM) > CFL_MAX)
      CFL_MAX = fabs(CFL_NUM);

    //#pragma omp parallel for
    // private(i,j,k,index,neut,up1,up2,down1,down2,slope)
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        for (k = 0; k < N; k++) {
          index = k + N * (j + N * i);
          // Poisson part
          if (CFL_NUM > 0) {
            neut = f[l][sp][index];
            if (i == 0) {
              up1 = f[l][sp][k + N * (j + N * (i + 1))];
              down1 = 0.0;
              down2 = 0.0;
            } else if (i == 1) {
              up1 = f[l][sp][k + N * (j + N * (i + 1))];
              down1 = f[l][sp][k + N * (j + N * (i - 1))];
              down2 = 0.0;
            } else if (i == N - 1) {
              up1 = 0.0;
              down1 = f[l][sp][k + N * (j + N * (i - 1))];
              down2 = f[l][sp][k + N * (j + N * (i - 2))];
            } else {
              up1 = f[l][sp][k + N * (j + N * (i + 1))];
              down1 = f[l][sp][k + N * (j + N * (i - 1))];
              down2 = f[l][sp][k + N * (j + N * (i - 2))];
            }

            slope[1] =
                minmod(2.0 * (up1 - neut) / h_v, 2.0 * (neut - down1) / h_v,
                       (up1 - down1) / (2.0 * h_v));
            slope[0] =
                minmod(2.0 * (neut - down1) / h_v, 2.0 * (down1 - down2) / h_v,
                       (neut - down2) / (2.0 * h_v));

            if (i != 0)
              f_conv[l][sp][index] =
                  -CFL_NUM *
                  (f[l][sp][index] + 0.5 * (h_v - E[l] * dt) * slope[1] -
                   (f[l][sp][k + N * (j + N * (i - 1))] +
                    0.5 * (h_v - E[l] * dt) * slope[0]));
            else
              f_conv[l][sp][index] =
                  -CFL_NUM *
                  (f[l][sp][index] + 0.5 * (h_v - E[l] * dt) * slope[1]);
          } else {
            neut = f[l][sp][index];
            if (i == 0) {
              up2 = f[l][sp][k + N * (j + N * (i + 2))];
              up1 = f[l][sp][k + N * (j + N * (i + 1))];
              down1 = 0.0;
            } else if (i == N - 1) {
              up2 = 0.0;
              up1 = 0.0;
              down1 = f[l][sp][k + N * (j + N * (i - 1))];
            } else if (i == N - 2) {
              up2 = 0.0;
              up1 = f[l][sp][k + N * (j + N * (i + 1))];
              down1 = f[l][sp][k + N * (j + N * (i - 1))];
            } else {
              up2 = f[l][sp][k + N * (j + N * (i + 2))];
              up1 = f[l][sp][k + N * (j + N * (i + 1))];
              down1 = f[l][sp][k + N * (j + N * (i - 1))];
            }

            slope[2] = minmod(2.0 * (up2 - up1) / h_v, 2.0 * (up1 - neut) / h_v,
                              (up2 - neut) / (2.0 * h_v));
            slope[1] =
                minmod(2.0 * (up1 - neut) / h_v, 2.0 * (neut - down1) / h_v,
                       (up1 - down1) / (2.0 * h_v));

            if (i != (N - 1))
              f_conv[l][sp][index] =
                  -CFL_NUM *
                  (f[l][sp][k + N * (j + N * (i + 1))] -
                   0.5 * (h_v + E[l] * dt) * slope[2] -
                   (f[l][sp][index] - 0.5 * (h_v + E[l] * dt) * slope[1]));
            else
              f_conv[l][sp][index] =
                  CFL_NUM *
                  (f[l][sp][index] - 0.5 * (h_v + E[l] * dt) * slope[1]);
          }
        }
  }

  if (CFL_MAX > 1.0) {
    printf("\n***********\n");
    printf("Warning - Electric field CFL is %g\n", CFL_MAX);
    printf("Code may be unstable\n");
    printf("\n***********\n");
  }

  free(E);
}

// Defunct, do not use. This is a placeholder for a first order time splitting,
// but with high-resolution methods for each piece
void advectTwo(double ***f, double *PoisPot, double **qm, double m, int sp) {
  upwindTwo_v(f, f_star, PoisPot, qm, m, c[sp], sp);
  upwindTwo_x(f_star, f, c[sp], sp);
}

// This returns the *update* to f as f_conv for the x transport piece
// used in second order time stepping method
void advectTwo_x(double ***f, double ***f_conv, int sp) {
  upwindTwo_x(f, f_conv, c[sp], sp);
}

// This returns the *update* to f as f_conv for the v transport piece
// used in second order time stepping method
void advectTwo_v(double ***f, double ***f_conv, double *PoisPot, double **qm,
                 double m, int sp) {
  upwindTwo_v(f, f_conv, PoisPot, qm, m, c[sp], sp);
}

void dealloc_trans() {
  int i, l;
  for (l = 0; l < nX; l++) {
    for (i = 0; i < ns; i++)
      free(f_star[l][i]);
    free(f_star[l]);
  }
  free(f_star);
}
