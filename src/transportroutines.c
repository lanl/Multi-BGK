#include "transportroutines.h"
#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

// Computes the transport term

static int N;
static int nX;
static int ns;
static int order;
static double Lx;
static double L_v;
static double h_v;
static double dt;
static int ICChoice;
static double TWall;
static double VWall;
static double T0, T1;
static double V0, V1;
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

// Computes the x direction first order upwind solution FOR A SINGLE SPECIES
void upwindOne_x(double ***f, double ***f_conv, double *v, int sp) {
  int i, j, k, l;
  int index;
  double CFL_NUM, CFL_NUM2;

  int rank, numNodes;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numNodes);
  MPI_Status status;

  for (l = 1; l < nX - 1; l++)
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

  // boundary cases
  //#pragma omp parallel for private(CFL_NUM,CFL_NUM2,i,j,k,index)
  for (i = 0; i < N; i++) {
    CFL_NUM = dt * v[i] / dx[0];
    CFL_NUM2 = dt * v[i] / dx[nX - 1];
    for (j = 0; j < N; j++)
      for (k = 0; k < N; k++) {
        index = k + N * (j + N * i);

        // the upwinding
        if (i < N / 2)
          f_conv[0][sp][index] =
              f[0][sp][index] - CFL_NUM * (f[1][sp][index] - f[0][sp][index]);
        else
          f_conv[0][sp][index] =
              f[0][sp][index] -
              CFL_NUM * (f[0][sp][index] - f[nX - 1][sp][index]);

        // the upwinding
        if (i < N / 2)
          f_conv[nX - 1][sp][index] =
              f[nX - 1][sp][index] -
              CFL_NUM2 * (f[0][sp][index] - f[nX - 1][sp][index]);
        else
          f_conv[nX - 1][sp][index] =
              f[nX - 1][sp][index] -
              CFL_NUM2 * (f[nX - 1][sp][index] - f[nX - 2][sp][index]);
      }
  }
}

// Computes the v direction first order upwind solution FOR A SINGLE SPECIES
void upwindOne_v(double ***f, double ***f_conv, double *PoisPot, double **qm,
                 double m, double *v, int sp) {
  int i, j, k, l;
  int index;
  double CFL_NUM;
  h_v = v[1] - v[0];

  int rank, numNodes;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numNodes);
  MPI_Status status;
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

  double *E = malloc(nX * sizeof(double));
  E[0] = -qm[0][sp] / m * (PoisPot[1] - PoisPot[nX - 1]) / (2.0 * dx[0]);
  for (i = 1; i < nX - 1; i++) {
    E[i] = -qm[i][sp] / m * (PoisPot[i + 1] - PoisPot[i - 1]) / (2.0 * dx[i]);
  }
  E[nX - 1] =
      -qm[nX - 1][sp] / m * (PoisPot[0] - PoisPot[nX - 2]) / (2.0 * dx[nX - 1]);

  for (l = 0; l < nX; l++)
    if (E[l] > 0) {
      CFL_NUM = dt * E[l] / h_v;
      //#pragma omp parallel for private(i,j,k,index)
      for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
          for (k = 0; k < N; k++) {
            index = k + N * (j + N * i);
            if (i != 0)
              f_conv[l][sp][index] =
                  f[l][sp][index] -
                  CFL_NUM *
                      (f[l][sp][index] - f[l][sp][k + N * (j + N * (i - 1))]);
            else
              f_conv[l][sp][index] =
                  f[l][sp][index] - CFL_NUM * (f[l][sp][index]);
          }
    } else {
      CFL_NUM = dt * E[l] / h_v;
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

  free(E);
}

// Computes x portion of second order upwind solution, with minmod
void upwindTwo_x(double ***f, double ***f_conv, double *v, int sp) {
  int i, j, k, l;
  int index;
  double slope[3];
  double CFL_NUM;

  int rank, numNodes;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numNodes);
  MPI_Status status;

  // main upwinding
  double f_l, f_r, f_ll, f_rr, x_l, x_r, x_ll, x_rr, dx_l, dx_r;

  for (l = 0; l < nX; l++) {
    if (l == 0) {
      x_l = x[nX - 1] - Lx;
      x_ll = x[nX - 2] - Lx;
      dx_l = dx[nX - 1];
    } else if (l == 1) {
      x_l = x[0];
      x_ll = x[nX - 1] - Lx;
      dx_l = dx[0];
    } else {
      x_l = x[l - 1];
      x_ll = x[l - 2];
      dx_l = dx[l - 1];
    }

    if (l == nX - 1) {
      x_r = x[0] + Lx;
      x_rr = x[1] + Lx;
      dx_r = dx[0];
    } else if (l == nX - 2) {
      x_r = x[nX - 1];
      x_rr = x[0] + Lx;
      dx_r = dx[nX - 1];
    } else {
      x_r = x[l + 1];
      x_rr = x[l + 2];
      dx_r = dx[l + 1];
    }
    //#pragma omp parallel for private(i,j,k,index,CFL_NUM,f_l,f_ll,f_r,slope)
    for (i = N / 2; i < N; i++) {
      CFL_NUM = dt * v[i] / dx[l];
      for (j = 0; j < N; j++)
        for (k = 0; k < N; k++) {
          index = k + N * (j + N * i);
          // upwind coming from the left
          // generate the local slopes
          if (l == 0) {
            f_l = f[nX - 1][sp][index];
            f_ll = f[nX - 2][sp][index];
          } else if (l == 1) {
            f_l = f[0][sp][index];
            f_ll = f[nX - 1][sp][index];
          } else {
            f_l = f[l - 1][sp][index];
            f_ll = f[l - 2][sp][index];
          }

          if (l == nX - 1)
            f_r = f[0][sp][index];
          else
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
          if (l == 0)
            f_l = f[nX - 1][sp][index];
          else
            f_l = f[l - 1][sp][index];

          if (l == nX - 1) {
            f_r = f[0][sp][index];
            f_rr = f[1][sp][index];
          } else if (l == nX - 2) {
            f_r = f[nX - 1][sp][index];
            f_rr = f[0][sp][index];
          } else {
            f_r = f[l + 1][sp][index];
            f_rr = f[l + 2][sp][index];
          }

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

// Computes seconrd order upwind term for the field acceleration term
void upwindTwo_v(double ***f, double ***f_conv, double *PoisPot, double **qm,
                 double m, double *v, int sp) {
  int i, j, k, l;
  int index;
  double slope[3];
  double CFL_NUM;

  int rank, numNodes;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numNodes);
  MPI_Status status;

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

  double *E = malloc(nX * sizeof(double));
  E[0] = -qm[0][sp] / m * (PoisPot[1] - PoisPot[nX - 1]) / (2.0 * dx[0]);
  for (i = 1; i < nX - 1; i++) {
    E[i] = -qm[i][sp] / m * (PoisPot[i + 1] - PoisPot[i - 1]) / (2.0 * dx[i]);
  }
  E[nX - 1] =
      -qm[nX - 1][sp] / m * (PoisPot[0] - PoisPot[nX - 2]) / (2.0 * dx[nX - 1]);

  if (fabs(v[2] - v[1] - h_v) > 1e-6) {
    printf("Error - you must use a uniform velocity grid for second order "
           "advection\n Please set discret=0 in your input file.\n");
    exit(1);
  }

  // main upwinding
  double f_l, f_r, f_ll, f_rr, x_l, x_r, x_ll, x_rr, dx_l, dx_r;
  double up2, up1, neut, down1, down2;

  for (l = 0; l < nX; l++) {
    CFL_NUM = dt * E[l] / h_v;
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

  free(E);
}

void initialize_transport(int numV, int numX, int nspec, double *xnodes,
                          double *dxnodes, double Lx_val, double **vel, int ord,
                          double timestep) {
  N = numV;
  nX = numX; // really NX NODE
  x = xnodes;
  dx = dxnodes;
  ns = nspec;
  Lx = Lx_val;

  c = vel;

  order = ord;
  dt = timestep;

  int i, l;
  f_star = malloc(nX * sizeof(double *));
  for (l = 0; l < nX; l++) {
    f_star[l] = malloc(nspec * sizeof(double *));
    for (i = 0; i < nspec; i++)
      f_star[l][i] = malloc(N * N * N * sizeof(double));
  }

  x = xnodes;
  dx = dxnodes;
}

// does a first order splitting, solving v first and returns the updated f
void advectOne(double ***f, double *PoisPot, double **qm, double m, int sp) {
  upwindOne_v(f, f_star, PoisPot, qm, m, c[sp], sp);
  upwindOne_x(f_star, f, c[sp], sp);
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
