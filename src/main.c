// Multi-BGK Code
// By Jeff Haack

// All units are CGS unless otherwise noted.

// C libraries
#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Utilities and setup stuff
#include "gauss_legendre.h"
#include "initialize_sol.h"
#include "momentRoutines.h"
#include "units/unit_data.c"
#include "zBar.h"

// I/O
#include "input.h"
#include "io.h"

// Vlasov/LHS packages
#include "mesh.h"
#include "poissonNonlinPeriodic.h"
#include "transportroutines.h"

// BGK/RHS packages
#include "BGK.h"

int main(int argc, char **argv) {

  MPI_Init(&argc, &argv);

  int rank, numRanks;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
  MPI_Status status;
  int rankCounter;
  double *momentBuffer;
  int *Nx_ranks;

  // get input information, set up the problem

  ////////////////
  // Declarations//
  ////////////////

  int nspec; // number of species
  int dims;  // flag for 0D or 1D

  int i, j, k, l, s;

  // species masses
  double *m; // in kg

  // species charges
  double **Z_oned, *Z_max, *Z_zerod;
  double E_max;

  // Flags for BGK collision rates
  int ecouple, CL_type, ion_type, MT_or_TR;
  int BGK_type;
  double beta;

  // distribution functions

  // For 0D - First dimension is species, second is velocity
  double **f_zerod;
  double **f_zerod_tmp;

  // For 1D - First dimension is species, second is position (x), third is
  // velocity
  double ***f, ***f_tmp, ***f_conv;

  // initial condition stuff - 0D
  double *v_val;

  double *n_zerod; // in cc
  double ntot, rhotot;

  double **v_zerod;
  double v0_zerod[3]; // mass-average velocity of the mixture

  double *T_zerod; // in eV

  double Htot, Htot_prev;
  double *H_spec, *H_spec_prev;
  double **BGK_f_minus_eq, **BGK_f_minus_eq_init;

  double T0;
  double Te;
  double *Te_arr;
  double Te_ref; // background electron temperature for interface problem
  double Te_start;

  // initial condition parameters - 1D
  int numint;
  double *intervalLimits;
  int input_file_data_flag = 0;

  double *vref;

  double *ndens_int;
  double **n_oned;

  double *velo_int;
  double ***v_oned;
  double **v0_oned;

  double *T_int;
  double **T_oned;
  double *T0_oned;
  double **T_for_zbar;
  double *T_max, T0_max;

  // Velocity grid setup - 0D and 1D
  int Nv;
  int discret;    // 0 - uniform, 1 - gauss
  double *Lv;     // semi-length, domain is [-Lv,Lv]
  double v_sigma; // width of velocity domains in terms of thermal speeds
  double vmax;
  double dvmin;
  double **c; // first dimension is species, second is 1D velo grid points. Full
              // velo pt is e.g. (v[0][i], v[0][j], v[0][k])
  double **wts; // first dimension is species, second is 1D velo grid points.
                // Full wts pt is e.g. (wts[0][i], wts[0][j], wts[0][k])
  int index;

  // Physical grid setup - 1D
  int order;   // spatial order of scheme, 1 or 2
  int Nx;      // Total number of physical grid points across all ranks
  int Nx_rank; // Physical grid points on the current rank (not including ghost
               // cells)
  double Lx;   // in cm
  double dx;
  double *x, *dxarray;

  double **dens; // for poisson
  double *PoisPot;
  double *source;
  int poissFlavor;

  // Time setup
  double dt;
  double tfinal;
  double t;
  int nT;
  double CFL1, CFL2, CFL3;
  double E_curr;
  double tau_min;
  int im_ex;
  int cfl_count = 0;
  int outcount = 0;
  int dataFreq;
  int outputDist;
  int restartFlag;
  int tauFlag;
  double RHS_tol;
  double RHS_min;
  double rel_BGK;

  int hydro_kinscheme_flag;

  /**********************************
         I/O Setup
  **********************************/

  char input_filename[100];
  char input_file_data_filename[100];
  strcpy(input_filename, argv[1]);

  read_input(&nspec, &dims, &Nx, &Lx, &Nv, &v_sigma, &discret, &poissFlavor, &m,
             &Z_max, &order, &im_ex, &dt, &tfinal, &numint, &intervalLimits,
             &ndens_int, &velo_int, &T_int, &ecouple, &Te_start, &Te_ref,
             &CL_type, &ion_type, &MT_or_TR, &n_zerod, &v_val, &T_zerod,
             &dataFreq, &outputDist, &RHS_tol, &BGK_type, &beta,
             &hydro_kinscheme_flag, &input_file_data_flag,
             input_file_data_filename, input_filename);

  char output_path[100] = {"./Data/"};
  strcat(output_path, argv[1]);

  if (argc < 3)
    restartFlag = 0;
  else
    restartFlag = atoi(argv[2]);

  if (argc < 4)
    tauFlag = 0;
  else
    tauFlag = atoi(argv[3]);

  if (argc > 2)
    printf("restart %d tau %d\n", restartFlag, tauFlag);

  // SHOULD FARM THIS OFF TO A ROUTINE
  char dens_path[100];
  char velo_path[100];
  char temp_path[100];

  char time_path[100];
  char H_path[100];
  char BGK_path[100];
  char poiss_path[100];
  char x_path[100];
  char name_tmp[100];

  FILE **outputFileDens = malloc(nspec * sizeof(FILE *));
  FILE **outputFileVelo = malloc(nspec * sizeof(FILE *));
  FILE **outputFileTemp = malloc(nspec * sizeof(FILE *));
  FILE *outputFileH;
  FILE *outputFileBGK;
  FILE *outputFileTime;
  FILE *outputFile_x;
  FILE *outputFilePoiss;

  H_spec = malloc(nspec * sizeof(double));
  H_spec_prev = malloc(nspec * sizeof(double));
  BGK_f_minus_eq = malloc(nspec * sizeof(double *));
  BGK_f_minus_eq_init = malloc(nspec * sizeof(double *));
  for (i = 0; i < nspec; i++) {
    BGK_f_minus_eq[i] = malloc(nspec * sizeof(double));
    BGK_f_minus_eq_init[i] = malloc(nspec * sizeof(double));
  }

  if ((dims == 0) || (rank == 0)) {
    for (i = 0; i < nspec; i++) {
      strcpy(dens_path, output_path);
      sprintf(name_tmp, "_dens%d", i);
      strcat(dens_path, name_tmp);
      if ((restartFlag == 2) || (restartFlag == 4))
        outputFileDens[i] = fopen(dens_path, "a");
      else
        outputFileDens[i] = fopen(dens_path, "w");

      strcpy(velo_path, output_path);
      sprintf(name_tmp, "_velo%d", i);
      strcat(velo_path, name_tmp);
      if ((restartFlag == 2) || (restartFlag == 4))
        outputFileVelo[i] = fopen(velo_path, "a");
      else
        outputFileVelo[i] = fopen(velo_path, "w");

      strcpy(temp_path, output_path);
      sprintf(name_tmp, "_temp%d", i);
      strcat(temp_path, name_tmp);
      if ((restartFlag == 2) || (restartFlag == 4))
        outputFileTemp[i] = fopen(temp_path, "a");
      else
        outputFileTemp[i] = fopen(temp_path, "w");
    }
  }

  // Set up file to store the times
  strcpy(time_path, output_path);
  sprintf(name_tmp, "_time");
  strcat(time_path, name_tmp);

  strcpy(H_path, output_path);
  sprintf(name_tmp, "_H");
  strcat(H_path, name_tmp);

  strcpy(BGK_path, output_path);
  sprintf(name_tmp, "_BGK");
  strcat(BGK_path, name_tmp);

  if ((dims == 0) || (rank == 0)) {

    if ((restartFlag == 2) || (restartFlag == 4))
      outputFileTime = fopen(time_path, "a");
    else
      outputFileTime = fopen(time_path, "w");

    if ((restartFlag == 2) || (restartFlag == 4)) {
      outputFileH = fopen(H_path, "a");
      outputFileBGK = fopen(BGK_path, "a");
    }

    else {
      outputFileTime = fopen(time_path, "w");
      outputFileH = fopen(H_path, "w");
      outputFileBGK = fopen(BGK_path, "w");
    }
  }

  if ((rank == 0) && (dims == 1)) {
    // Set up file to store E-field
    strcpy(poiss_path, output_path);
    strcat(poiss_path, "_poiss");

    if ((restartFlag == 2) || (restartFlag == 4))
      outputFilePoiss = fopen(poiss_path, "a");
    else
      outputFilePoiss = fopen(poiss_path, "w");

    // store physical mesh
    FILE *outputFile_x;
    strcpy(x_path, output_path);
    strcat(x_path, "_x");
    outputFile_x = fopen(x_path, "w");
  }

  printf("Input file: %s\n", input_filename);
  printf("Output filename: %s\n", output_path);

  ///////////////////
  // PARAMETER SETUP//
  ///////////////////

  //// 0D Case ////
  if (dims == 0) {

    if (numRanks > 1) {
      printf("Error - only use multiple MPI ranks if running a 1D problem");
      exit(37);
    }

    v0_zerod[0] = 0.0;
    v0_zerod[1] = 0.0;
    v0_zerod[2] = 0.0;
    v_zerod = malloc(nspec * sizeof(double *));
    rhotot = 0.0;
    ntot = 0.0;

    for (i = 0; i < nspec; i++) {
      v_zerod[i] = malloc(3 * sizeof(double));
      v_zerod[i][0] = v_val[i];
      v_zerod[i][1] = 0.0;
      v_zerod[i][2] = 0.0;
      v0_zerod[0] += m[i] * n_zerod[i] * v_zerod[i][0];
      ntot += n_zerod[i];
      rhotot += m[i] * n_zerod[i];
    }
    v0_zerod[0] = v0_zerod[0] / rhotot;

    T0 = 0.0;
    for (i = 0; i < nspec; i++) {
      if (n_zerod[i] != 0.0) {
        T0 += n_zerod[i] * T_zerod[i] / ntot;
        for (j = 0; j < 3; j++)
          T0 += m[i] * n_zerod[i] * (v_zerod[i][j] - v0_zerod[j]) *
                (v_zerod[i][j] - v0_zerod[j]) / (3.0 * ntot);
      }
    }
    if (ecouple == 1) {
      if (T0 < Te_ref)
        T0 = Te_ref;
    }

    Z_zerod = malloc(nspec * sizeof(double));

    // set the velocity grids
    vref = malloc(nspec * sizeof(double));
    Lv = malloc(nspec * sizeof(double));
    c = malloc(nspec * sizeof(double *));
    wts = malloc(nspec * sizeof(double *));

    for (i = 0; i < nspec; i++) {
      c[i] = malloc(Nv * sizeof(double));
      wts[i] = malloc(Nv * sizeof(double));
    }

    io_init_homog(Nv, nspec, c);

    double dv;

    if ((restartFlag == 2) || (restartFlag == 4)) {
      load_grid_restart(
          Lv, &t, input_filename); // note - this assumes uniform grid for now
      for (i = 0; i < nspec; i++) {
        dv = 2.0 * Lv[i] / (Nv - 1.0);
        for (j = 0; j < Nv; j++) {
          c[i][j] = -Lv[i] + dv * j;
          wts[i][j] = dv;
        }
        wts[i][0] *= 0.5;
        wts[i][Nv - 1] *= 0.5;
      }
    } else { // No restart, build the grid from scratch
      for (i = 0; i < nspec; i++) {
        if (T0 > T_zerod[i])
          vref[i] = sqrt((T0 / ERG_TO_EV_CGS) / m[i]);
        else
          vref[i] = sqrt((T_zerod[i] / ERG_TO_EV_CGS) / m[i]);
        Lv[i] = v_sigma * vref[i];
      }

      if (discret == 0) { // uniform velocity grid
        for (i = 0; i < nspec; i++) {
          if (rank == 0) {
            printf("Setting up uniform velocity grid for species %d with temps "
                   "%g %g \n",
                   i, T_zerod[i], T0);
          }
          dv = 2.0 * Lv[i] / (Nv - 1.0);
          for (j = 0; j < Nv; j++) {
            c[i][j] = -Lv[i] + dv * j;
            wts[i][j] = dv;
          }
          wts[i][0] *= 0.5;
          wts[i][Nv - 1] *= 0.5;
        }
      } else { // Using Gauss-Legendre velocity grid
        double *GLGrid = malloc(Nv * sizeof(double));
        double *GLWeights = malloc(Nv * sizeof(double));
        double A, B;
        int GLCenter;

        // set up gauss grid
        gauss_legendre_tbl(Nv, GLGrid, GLWeights, 1e-10);

        for (i = 0; i < nspec; i++) {
          if (rank == 0) {
            printf("Setting up gauss point velo grid for species %d with temps "
                   "%g %g \n",
                   i, T_zerod[i], T0);
          }

          A = Lv[i];
          B = 0.0;
          GLCenter = Nv / 2;
          if (2 * GLCenter < Nv) { // odd
            c[i][GLCenter] = B;
            wts[i][GLCenter] = A * GLWeights[0];
            for (j = 1; j <= GLCenter; j++) {
              c[i][GLCenter + j] = B + A * GLGrid[j];
              c[i][GLCenter - j] = B - A * GLGrid[j];
              wts[i][GLCenter + j] = A * GLWeights[j];
              wts[i][GLCenter - j] = A * GLWeights[j];
            }
          } else { // even
            for (j = 0; j < GLCenter; j++) {
              c[i][GLCenter + j] = B + A * GLGrid[j];
              c[i][GLCenter - j - 1] = B - A * GLGrid[j];
              wts[i][GLCenter + j] = A * GLWeights[j];
              wts[i][GLCenter - j - 1] = A * GLWeights[j];
            }
          }
        }
        free(GLGrid);
        free(GLWeights);
      }
    }

    f_zerod = malloc(nspec * sizeof(double *));
    f_zerod_tmp = malloc(nspec * sizeof(double *));
    for (i = 0; i < nspec; i++) {
      f_zerod[i] = malloc(Nv * Nv * Nv * sizeof(double));
      f_zerod_tmp[i] = malloc(Nv * Nv * Nv * sizeof(double));
    }
    if (rank == 0) {
      printf("Done allocating for 0D\n");
    }
  }

  if (dims == 1) {

    // Physical grid allocation and initialization
    make_mesh(Nx, Lx, order, &Nx_rank, &Nx_ranks, &x, &dxarray);
    momentBuffer = malloc(3 * (Nx_rank + 1) * sizeof(double));

    dx = Lx / Nx;

    /*********
    //Old style
    dx = Lx/Nx;
    x = malloc(Nx*sizeof(double));
    dxarray = malloc(Nx*sizeof(double));
    *********/

    // allocate rank-local moment arrays
    n_oned = malloc(Nx_rank * sizeof(double *));
    v_oned = malloc(Nx_rank * sizeof(double **));
    T_oned = malloc(Nx_rank * sizeof(double *));

    v0_oned = malloc(Nx_rank * sizeof(double *));
    T0_oned = malloc(Nx_rank * sizeof(double));

    T_max = malloc(nspec * sizeof(double));
    for (i = 0; i < nspec; i++)
      T_max[i] = 0.0;
    T0_max = 0.0;

    for (l = 0; l < Nx_rank; l++) {
      // fprintf(outputFile_x,"%+le ",x[l]);

      n_oned[l] = malloc(nspec * sizeof(double));
      v_oned[l] = malloc(nspec * sizeof(double *));
      for (i = 0; i < nspec; i++)
        v_oned[l][i] = malloc(3 * sizeof(double));
      T_oned[l] = malloc(nspec * sizeof(double));

      v0_oned[l] = malloc(3 * sizeof(double));

      v0_oned[l][0] = 0.0;
      v0_oned[l][1] = 0.0;
      v0_oned[l][2] = 0.0;

      T0_oned[l] = 0.0;
    }

    if (input_file_data_flag) {
      if (rank == 0) {
        printf("input_file_data_flag: %d, filename: %s\n", input_file_data_flag,
               input_file_data_filename);
      }

      initialize_sol_load_inhom_file(Nx, nspec, n_oned, v_oned, T_oned,
                                     input_file_data_filename);

      // Find maximum temperature
      for (l = 0; l < Nx_rank; l++)
        for (s = 0; s < nspec; s++) {
          T_max[s] = (T_oned[l][s] > T_max[s]) ? T_oned[l][s] : T_max[s];
          T0_max = (T_max[s] > T0_max) ? T_max[s] : T0_max;
        }

      if ((T0_max < Te_ref) && (ecouple == 1))
        T0_max = Te_ref;

    } else {
      for (i = 0; i < numint; i++) {
        for (j = 0; j < nspec; j++) {
          T_max[j] = (T_int[i * nspec + j] > T_max[j]) ? T_int[i * nspec + j]
                                                       : T_max[j];
          T0_max = (T_max[j] > T0_max) ? T_max[j] : T0_max;
        }
      }
      if ((T0_max < Te_ref) && (ecouple == 1))
        T0_max = Te_ref;
    }

    // Allocate for poisson calc
    PoisPot = malloc(Nx_rank * sizeof(double));
    source = malloc(Nx_rank * sizeof(double));
    T_for_zbar = malloc(Nx_rank * sizeof(double));
    Te_arr = malloc(Nx_rank * sizeof(double));

    Z_oned = malloc(Nx_rank * sizeof(double *));

    for (l = 0; l < Nx_rank; l++) {
      Z_oned[l] = malloc(nspec * sizeof(double));
      for (i = 0; i < nspec; i++)
        Z_oned[l][i] = Z_max[i]; // initialized at full ionization
    }

    // Velocity grid allocation and initialization

    // set the velocity grids
    vref = malloc(nspec * sizeof(double));
    Lv = malloc(nspec * sizeof(double));
    c = malloc(nspec * sizeof(double *));
    wts = malloc(nspec * sizeof(double *));

    for (i = 0; i < nspec; i++) {
      if (T0_max > T_max[i])
        vref[i] = sqrt((T0_max / ERG_TO_EV_CGS) / m[i]);
      else
        vref[i] = sqrt((T_max[i] / ERG_TO_EV_CGS) / m[i]);

      Lv[i] = v_sigma * vref[i];
      c[i] = malloc(Nv * sizeof(double));
      wts[i] = malloc(Nv * sizeof(double));
    }

    if (discret == 0) { // uniform velocity grid
      double dv;
      for (i = 0; i < nspec; i++) {
        if (rank == 0) {
          printf("Setting up uniform velocity grid for species %d with temps "
                 "%g %g \n",
                 i, T_max[i], T0_max);
        }
        dv = 2.0 * Lv[i] / (Nv - 1.0);
        for (j = 0; j < Nv; j++) {
          c[i][j] = -Lv[i] + dv * j;
          wts[i][j] = dv;
        }
        wts[i][0] *= 0.5;
        wts[i][Nv - 1] *= 0.5;
      }
    } else { // Using Gauss-Legendre velocity grid
      double *GLGrid = malloc(Nv * sizeof(double));
      double *GLWeights = malloc(Nv * sizeof(double));
      double A, B;
      int GLCenter;

      // set up gauss grid
      gauss_legendre_tbl(Nv, GLGrid, GLWeights, 1e-10);

      for (i = 0; i < nspec; i++) {
        if (rank == 0) {
          printf("Setting up gauss point velo grid for species %d with temps "
                 "%g %g \n",
                 i, T_max[i], T0_max);
        }

        A = Lv[i];
        B = 0.0;
        GLCenter = Nv / 2;
        if (2 * GLCenter < Nv) { // odd
          c[i][GLCenter] = B;
          wts[i][GLCenter] = A * GLWeights[0];
          for (j = 1; j <= GLCenter; j++) {
            c[i][GLCenter + j] = B + A * GLGrid[j];
            c[i][GLCenter - j] = B - A * GLGrid[j];
            wts[i][GLCenter + j] = A * GLWeights[j];
            wts[i][GLCenter - j] = A * GLWeights[j];
          }
        } else { // even
          for (j = 0; j < GLCenter; j++) {
            c[i][GLCenter + j] = B + A * GLGrid[j];
            c[i][GLCenter - j - 1] = B - A * GLGrid[j];
            wts[i][GLCenter + j] = A * GLWeights[j];
            wts[i][GLCenter - j - 1] = A * GLWeights[j];
          }
        }
      }
      free(GLGrid);
      free(GLWeights);
    }

    io_init_inhomog(Nx_rank, Nv, nspec, c);
    if (outputDist == 1) {
      store_grid(input_filename);
    }

    // Distribution function setup

    f = malloc((Nx_rank + 2 * order) * sizeof(double *));
    f_tmp = malloc((Nx_rank + 2 * order) * sizeof(double *));
    f_conv = malloc((Nx_rank + 2 * order) * sizeof(double *));

    for (l = 0; l < (Nx_rank + 2 * order); l++) {
      f[l] = malloc(nspec * sizeof(double *));
      f_tmp[l] = malloc(nspec * sizeof(double *));
      f_conv[l] = malloc(nspec * sizeof(double *));
      for (i = 0; i < nspec; i++) {
        f[l][i] = malloc(Nv * Nv * Nv * sizeof(double));
        f_tmp[l][i] = malloc(Nv * Nv * Nv * sizeof(double));
        f_conv[l][i] = malloc(Nv * Nv * Nv * sizeof(double));
      }
    }
  }

  //////////////////////
  // INITIAL CONDITIONS//
  //////////////////////

  // initial condition initialization

  if (dims == 0) {

    printf("Setting initial conditions for 0D\n");

    if ((restartFlag == 2) || (restartFlag == 4)) {
      load_distributions_homog(f_zerod, input_filename);
      printf("Loaded distribution from %s.dat\n", input_filename);
    } else {
      int index;
      for (l = 0; l < nspec; l++) {
        if (v_zerod[l][0] > Lv[l]) { // adjust velocity grid...
          printf("Shifting by %g\n", v_zerod[l][0]);
          for (i = 0; i < Nv; i++)
            c[l][i] += v_zerod[l][0];
        }
        for (i = 0; i < Nv; i++)
          for (j = 0; j < Nv; j++)
            for (k = 0; k < Nv; k++) {
              index = k + Nv * (j + Nv * i);
              f_zerod[l][index] =
                  n_zerod[l] *
                  pow(m[l] / (2.0 * M_PI * T_zerod[l] / ERG_TO_EV_CGS), 1.5) *
                  exp(-(0.5 * m[l] / (T_zerod[l] / ERG_TO_EV_CGS)) *
                      ((c[l][i] - v_zerod[l][0]) * (c[l][i] - v_zerod[l][0]) +
                       (c[l][j] - v_zerod[l][1]) * (c[l][j] - v_zerod[l][1]) +
                       (c[l][k] - v_zerod[l][2]) * (c[l][k] - v_zerod[l][2])));
            }
      }
    }

    printf("Done setting initial conditions for 0D\n");
  }

  if (dims == 1) {

    initialize_transport(Nv, Nx_rank, nspec, x, dxarray, Lx, c, order, dt);

    // check if we are loading moment data from a file
    // Note: This is not yet MPI-ified
    if (input_file_data_flag) {
      initialize_sol_inhom_file(f, Nx_rank, nspec, Nv, order, c, m, n_oned,
                                v_oned, T_oned);
    } else
      initialize_sol_inhom(f, numint, intervalLimits, ndens_int, velo_int,
                           T_int, Nx_rank, x, nspec, Nv, order, c, m, n_oned,
                           v_oned, T_oned);

    if (rank == 0) {
      printf("Initial condition setup complete\n");
      fflush(stdout);
    }
  }

  initialize_moments(Nv, nspec, c, wts);
  initialize_BGK(nspec, Nv, m, c, order, ecouple, CL_type, ion_type, MT_or_TR,
                 tauFlag, input_filename);

  if (!((restartFlag == 2) || (restartFlag == 4)))
    t = 0.0;

  if (restartFlag > 2) {
    // Check to see if we need to store initial BGK error data
    if (rank == 0) {
      printf("Setting initial BGK norm\n");
    }
    zBarFunc2(nspec, T0, Z_max, n_zerod, Z_zerod);
    BGK_norm(f_zerod, BGK_f_minus_eq_init, Z_zerod, dt, T0);
  }

  nT = 0;
  Htot_prev = 0.0;
  Htot = 0.0;
  for (i = 0; i < nspec; i++) {
    H_spec[i] = 0;
  }

  // Setup complete, begin main loop
  while (t < tfinal) {

    if (rank == 0) {
      printf("At time %g of %g\n", t, tfinal);
    }

    if (dims == 0) {
      // GET MOMENTS

      ntot = 0.0;
      rhotot = 0.0;
      Htot_prev = Htot;
      Htot = 0.0;

      for (i = 0; i < nspec; i++) {
        n_zerod[i] = getDensity(f_zerod[i], i);
        ntot += n_zerod[i];
        rhotot += m[i] * n_zerod[i];

        H_spec_prev[i] = H_spec[i];
        H_spec[i] = getH(n_zerod[i], f_zerod[i], i);
        Htot += H_spec[i];

        getBulkVel(f_zerod[i], v_zerod[i], n_zerod[i], i);
      }

      // get mixture mass avg velocity
      for (j = 0; j < 3; j++) {
        v0_zerod[j] = 0.0;
        for (i = 0; i < nspec; i++)
          v0_zerod[j] += m[i] * n_zerod[i] * v_zerod[i][j];
        v0_zerod[j] = v0_zerod[j] / rhotot;
      }

      // Find temperatures BASED ON INDIVIDUAL SPECIES VELOCITY. Note - result
      // is in eV
      for (i = 0; i < nspec; i++)
        T_zerod[i] = getTemp(m[i], n_zerod[i], v_zerod[i], f_zerod[i], i);

      // calculate mixture temperature (not the same as what goes into the
      // Maxwellians)
      T0 = 0.0;
      for (i = 0; i < nspec; i++) {
        if (n_zerod[i] != 0.0) {
          T0 += n_zerod[i] * T_zerod[i];
          for (j = 0; j < 3; j++)
            T0 += ERG_TO_EV_CGS * m[i] * n_zerod[i] *
                  (v_zerod[i][j] - v0_zerod[j]) *
                  (v_zerod[i][j] - v0_zerod[j]) / 3.0;
        }
      }
      T0 = T0 / ntot;

      // note - T0 is used as the electron temperature
      if (ecouple == 1)
        T0 = Te_ref;
      else if (ecouple == 2)
        T0 = T_zerod[0];

      // calc Zbar
      if (ecouple != 2)
        zBarFunc2(nspec, T0, Z_max, n_zerod, Z_zerod);

      // Output to files

      BGK_norm(f_zerod, BGK_f_minus_eq, Z_zerod, dt, T0);

      outcount += 1;
      if (outcount == dataFreq) {
        fprintf(outputFileTime, "%le\n", t);
        for (i = 0; i < nspec; i++) {
          fprintf(outputFileDens[i], "%le ", n_zerod[i]);
          fprintf(outputFileVelo[i], "%le ", v_zerod[i][0]);
          fprintf(outputFileTemp[i], "%le ", T_zerod[i]);
          fprintf(outputFileH, "%le ", (H_spec[i] - H_spec_prev[i]) / dt);
          fprintf(outputFileDens[i], "\n");
          fprintf(outputFileVelo[i], "\n");
          fprintf(outputFileTemp[i], "\n");
          outcount = 0;
        }
        fprintf(outputFileH, "%le\n", (Htot - Htot_prev) / dt);

        for (i = 0; i < nspec; i++)
          for (j = 0; j < nspec; j++)
            fprintf(outputFileBGK, "%le,", BGK_f_minus_eq[i][j]);
        fprintf(outputFileBGK, "\n");
      }

      // check to make sure that we don't need to stop
      RHS_min = 1.0e37;
      if (restartFlag > 2) {
        for (i = 0; i < nspec; i++)
          for (j = 0; j < nspec; j++) {
            if (i != j) {
              rel_BGK = BGK_f_minus_eq[i][j] / BGK_f_minus_eq_init[i][j];
              RHS_min = rel_BGK < RHS_min ? rel_BGK : RHS_min;
            }
          }

        // Stop computation
        if (RHS_min < RHS_tol) {
          printf("Stale collision rate detected\n");
          printf("RHS_min %g\n", RHS_min);
          for (i = 0; i < nspec; i++)
            for (j = 0; j < nspec; j++)
              printf("%d %d %g\n", i, j,
                     BGK_f_minus_eq[i][j] / BGK_f_minus_eq_init[i][j]);
          break;
        }
      }

      if (im_ex == 0) {

        // pick yer poison
        if (BGK_type == 0)
          BGK_ex(f_zerod, f_zerod_tmp, Z_zerod, dt, T0);
        else if (BGK_type == 1)
          BGK_Greene(f_zerod, f_zerod_tmp, Z_max, dt, beta, T0);
        else if (BGK_type == 2)
          BGK_NRL(f_zerod, f_zerod_tmp, Z_max, dt, T0);

        for (i = 0; i < nspec; i++)
          for (j = 0; j < Nv * Nv * Nv; j++)
            f_zerod[i][j] += dt * f_zerod_tmp[i][j];
      } else {
        BGK_im(f_zerod, Z_max, dt, T0);
      }
    } else if (dims == 1) {

      outcount += 1;

      // Calculate moment data in all cells
      for (l = 0; l < Nx_rank; l++) {
        ntot = 0.0;
        rhotot = 0.0;
        for (i = 0; i < nspec; i++) {
          n_oned[l][i] = getDensity(f[l + order][i], i);
          ntot += n_oned[l][i];
          rhotot += m[i] * n_oned[l][i];
          getBulkVel(f[l + order][i], v_oned[l][i], n_oned[l][i], i);
        }

        // get mixture mass avg velocity
        for (j = 0; j < 3; j++) {
          v0_oned[l][j] = 0.0;
          for (i = 0; i < nspec; i++)
            v0_oned[l][j] += m[i] * n_oned[l][i] * v_oned[l][i][j];
          v0_oned[l][j] = v0_oned[l][j] / rhotot;
        }
      }

      // Set T0 in all cells
      for (l = 0; l < Nx_rank; l++) {
        if (ecouple == 2)
          T0_oned[l] = T_oned[l][0];
        else {
          T0_oned[l] = 0.0;
          for (i = 0; i < nspec; i++) {
            if (n_oned[l][i] != 0.0) {
              T0_oned[l] += n_oned[l][i] * T_oned[l][i];
              for (j = 0; j < 3; j++)
                T0_oned[l] += ERG_TO_EV_CGS * m[i] * n_oned[l][i] *
                              (v_oned[l][i][j] - v0_oned[l][j]) *
                              (v_oned[l][i][j] - v0_oned[l][j]) / 3.0;
            }
          }
          T0_oned[l] = T0_oned[l] / ntot;
        }
      }

      // Flag - do we want to run this like the kinetic scheme for hydro
      if (hydro_kinscheme_flag == 1) {
        for (l = 0; l < Nx_rank; l++)
          for (i = 0; i < nspec; i++)
            GetMaxwell(m[i], n_oned[l][i], v_oned[l][i], T_oned[l][i],
                       f[l + order][i], i);
      }

      /**************
       ELECTRIC FIELD CALCULATION
       Note that the poisson solve gives e*phi (in ergs), not phi, more
      convenient for units issues
      ************/

      if (ecouple == 1) { // electrons only in background
        if (Te_start != Te_ref)
          get_ramp_Te(Te_arr, Nx_rank, Te_start, Te_ref, t, tfinal);
        else
          get_uniform_Te(Te_arr, Nx_rank,
                         Te_ref); // fixed background temperature
      } else
        Te_arr = T0_oned;

      if (ecouple == 2) {
        for (l = 0; l < Nx_rank; l++) {
          source[l] = 0.0;
          for (i = 1; i < nspec; i++) {
            source[l] +=
                Z_oned[l][i] *
                n_oned[l][i]; // total number of free electrons in each cell
          }
          source[l] -= n_oned[l][0]; // electrons
        }

      } else {
        for (l = 0; l < Nx_rank; l++) {
          source[l] = 0.0;
          // Get ionization
          zBarFunc2(nspec, Te_arr[l], Z_max, n_oned[l], Z_oned[l]);
          for (i = 0; i < nspec; i++) {
            source[l] +=
                Z_oned[l][i] *
                n_oned[l][i]; // total number of free electrons in each cell
          }
        }
      }

      if (ecouple == 2) {
        simplePoisson(Nx_rank, source, dx, Lx, PoisPot);
      } else {
        if (poissFlavor == 0) { // no E-field
          for (l = 0; l < Nx_rank; l++)
            PoisPot[l] = 0.0;
        } else if (poissFlavor == 11) // Linear Yukawa
          PoissLinPeriodic1D(Nx_rank, source, dx, Lx, PoisPot, Te_arr);
        else if (poissFlavor == 12) // Nonlinear Yukawa
          PoissNonlinPeriodic1D(Nx_rank, source, dx, Lx, PoisPot, Te_arr);
        else if (poissFlavor == 21) // Linear Thomas-Fermi
          PoissLinPeriodic1D_TF(Nx_rank, source, dx, Lx, PoisPot, Te_arr);
        else if (poissFlavor == 22) // Nonlinear Thomas-Fermi
          PoissNonlinPeriodic1D_TF(Nx_rank, source, dx, Lx, PoisPot, Te_arr);
      }

      if ((numRanks > 1) && (poissFlavor != 0)) {
        printf("Error - MPI not implemented for poisson solve");
        exit(1);
      }

      // Moments and initial electric field calculated - save if needed

      // MPI-ified output.
      if ((outcount == dataFreq) || t == 0.0) {
        if (rank == 0) {
          for (l = 0; l < Nx_rank; l++) {
            for (i = 0; i < nspec; i++) {
              T_oned[l][i] =
                  getTemp(m[i], n_oned[l][i], v_oned[l][i], f[l + order][i], i);
              if (outcount == dataFreq) {
                fprintf(outputFileDens[i], "%le ", n_oned[l][i]);
                fprintf(outputFileVelo[i], "%le ", v_oned[l][i][0]);
                fprintf(outputFileTemp[i], "%le ", T_oned[l][i]);
              }
            }
          }

          // get from other ranks
          for (rankCounter = 1; rankCounter < numRanks; rankCounter++) {
            printf("Calculating rank %d quantities\n", rankCounter);
            for (s = 0; s < nspec; s++) {
              MPI_Recv(momentBuffer, 3 * Nx_ranks[rankCounter], MPI_DOUBLE,
                       rankCounter, s, MPI_COMM_WORLD, &status);
              for (l = 0; l < Nx_ranks[rankCounter]; l++) {
                fprintf(outputFileDens[s], "%le ", momentBuffer[0 + 3 * l]);
                fprintf(outputFileVelo[s], "%le ", momentBuffer[1 + 3 * l]);
                fprintf(outputFileTemp[s], "%le ", momentBuffer[2 + 3 * l]);
              }
            }
          }

          // Close out this timestep
          fprintf(outputFileTime, "%le\n", t);
          for (i = 0; i < nspec; i++) {
            fprintf(outputFileDens[i], "\n");
            fprintf(outputFileVelo[i], "\n");
            fprintf(outputFileTemp[i], "\n");
            fprintf(outputFilePoiss, "\n");
          }

          if (outputDist == 1)
            store_distributions_inhomog(f, input_filename, nT);

          outcount = 0;
        }

        else { // send to rank 0 for output purposes
          for (s = 0; s < nspec; s++) {
            for (l = 0; l < Nx_rank; l++) {
              momentBuffer[0 + 3 * l] = n_oned[l][s];
              momentBuffer[1 + 3 * l] = v_oned[l][s][0];
              momentBuffer[2 + 3 * l] = T_oned[l][s];
            }
            MPI_Send(momentBuffer, 3 * Nx_rank, MPI_DOUBLE, 0, s,
                     MPI_COMM_WORLD);
          }

          outcount = 0;
        }
      }
      // IO done, advance to the actual solution...

      if (order == 1) {
        // ADVECT

        for (i = 0; i < nspec; i++) {
          advectOne(f, PoisPot, Z_oned, m[i], i);
        }

        // COLLIDE

        for (l = 0; l < Nx_rank; l++) {
          BGK_ex(f[l + order], f_tmp[l + order], Z_oned[l], dt, Te_arr[l]);
          for (i = 0; i < nspec; i++)
            for (j = 0; j < Nv * Nv * Nv; j++)
              f[l + order][i][j] += dt * f_tmp[l + order][i][j];
        }
      }

      if (order == 2) {
        // ADVECT

        // First strang step - V advection with timestep dt/2

        // SSP RK2 method:
        // u^(1) = u^n + dt*f(u^n)
        // u^(2) = 0.5*(u^n + u^(1)) + 0.5*dt*f(u^(1))
        // u^n+1 = u^(2)
        // but we are doing dt/2 steps due to strang

        // RK2 Step 1
        for (i = 0; i < nspec; i++) {
          advectTwo_v(f, f_conv, PoisPot, Z_oned, m[i], i);
        }

        for (l = 0; l < Nx_rank; l++)
          for (i = 0; i < nspec; i++)
            //#pragma omp parallel for private(j)
            for (j = 0; j < Nv * Nv * Nv; j++)
              f_tmp[l + order][i][j] =
                  f[l + order][i][j] + 0.5 * f_conv[l + order][i][j];

        if (ecouple == 1) { // electrons only in background

          if (Te_start != Te_ref)
            get_ramp_Te(Te_arr, Nx_rank, Te_start, Te_ref, t, tfinal);
          else
            get_uniform_Te(Te_arr, Nx_rank,
                           Te_ref); // fixed background temperature

        } else
          Te_arr = T0_oned;

        // Recalc Poiss
        if (ecouple == 2) {
          for (l = 0; l < Nx_rank; l++) {
            source[l] = 0.0;

            for (i = 0; i < nspec; i++)
              n_oned[l][i] = getDensity(f_tmp[l + order][i], i);

            for (i = 1; i < nspec; i++)
              source[l] +=
                  Z_oned[l][i] *
                  n_oned[l][i]; // total number of free electrons in each cell
            source[l] -= n_oned[l][0];
          }

        } else {
          for (l = 0; l < Nx_rank; l++) {
            for (i = 0; i < nspec; i++)
              n_oned[l][i] = getDensity(f_tmp[l + order][i], i);
            zBarFunc2(nspec, Te_arr[l], Z_max, n_oned[l], Z_oned[l]);

            source[l] = 0.0;
            for (i = 0; i < nspec; i++)
              source[l] +=
                  Z_oned[l][i] *
                  n_oned[l][i]; // total number of free electrons in each cell
          }
        }

        if (ecouple == 2) {
          simplePoisson(Nx_rank, source, dx, Lx, PoisPot);
        } else {
          if (poissFlavor == 0) { // no E-field
            for (l = 0; l < Nx_rank; l++)
              PoisPot[l] = 0.0;
          } else if (poissFlavor == 11) // Linear Yukawa
            PoissLinPeriodic1D(Nx_rank, source, dx, Lx, PoisPot, Te_arr);
          else if (poissFlavor == 12) // Nonlinear Yukawa
            PoissNonlinPeriodic1D(Nx_rank, source, dx, Lx, PoisPot, Te_arr);
          else if (poissFlavor == 21) // Linear Thomas-Fermi
            PoissLinPeriodic1D_TF(Nx_rank, source, dx, Lx, PoisPot, Te_arr);
          else if (poissFlavor == 22) // Nonlinear Thomas-Fermi
            PoissNonlinPeriodic1D_TF(Nx_rank, source, dx, Lx, PoisPot, Te_arr);
        }

        // RK2 Step 2
        for (i = 0; i < nspec; i++) {
          advectTwo_v(f_tmp, f_conv, PoisPot, Z_oned, m[i], i);
        }

        for (l = 0; l < Nx_rank; l++)
          for (i = 0; i < nspec; i++)
            //#pragma omp parallel for private(j)
            for (j = 0; j < Nv * Nv * Nv; j++)
              f[l + order][i][j] =
                  0.5 * (f[l + order][i][j] + f_tmp[l + order][i][j]) +
                  0.25 * f_conv[l + order][i][j];

        // Next strang step - x advection with timestep dt/2

        // RK2 Step 1
        for (i = 0; i < nspec; i++) {
          advectTwo_x(f, f_conv, i);
        }

        for (l = 0; l < Nx_rank; l++)
          for (i = 0; i < nspec; i++)
            //#pragma omp parallel for private(j)
            for (j = 0; j < Nv * Nv * Nv; j++)
              f_tmp[l + order][i][j] =
                  f[l + order][i][j] + 0.5 * f_conv[l + order][i][j];

        // RK2 Step 2
        for (i = 0; i < nspec; i++) {
          advectTwo_x(f_tmp, f_conv, i);
        }

        for (l = 0; l < Nx_rank; l++)
          for (i = 0; i < nspec; i++)
            //#pragma omp parallel for private(j)
            for (j = 0; j < Nv * Nv * Nv; j++)
              f[l + order][i][j] =
                  0.5 * (f[l + order][i][j] + f_tmp[l + order][i][j]) +
                  0.25 * f_conv[l + order][i][j];

        // Next Strang step - RK2 for collision with timstep dt

        for (l = 0; l < Nx_rank; l++) {

          // Step 1
          BGK_ex(f[l + order], f_conv[l + order], Z_oned[l], dt, Te_arr[l]);
          for (i = 0; i < nspec; i++)
            for (j = 0; j < Nv * Nv * Nv; j++)
              f_tmp[l + order][i][j] =
                  f[l + order][i][j] + dt * f_conv[l + order][i][j];

          // Step 2
          BGK_ex(f_tmp[l + order], f_conv[l + order], Z_oned[l], dt, Te_arr[l]);
          for (i = 0; i < nspec; i++)
            for (j = 0; j < Nv * Nv * Nv; j++)
              f[l + order][i][j] =
                  0.5 * (f[l + order][i][j] + f_tmp[l + order][i][j]) +
                  0.5 * dt * f_conv[l + order][i][j];
        }

        // Next strang step - x advection with timestep dt/2

        // RK2 Step 1
        for (i = 0; i < nspec; i++) {
          advectTwo_x(f, f_conv, i);
        }

        for (l = 0; l < Nx_rank; l++)
          for (i = 0; i < nspec; i++)
            //#pragma omp parallel for private(j)
            for (j = 0; j < Nv * Nv * Nv; j++)
              f_tmp[l + order][i][j] =
                  f[l + order][i][j] + 0.5 * f_conv[l + order][i][j];

        // RK2 Step 2
        for (i = 0; i < nspec; i++) {
          advectTwo_x(f_tmp, f_conv, i);
        }
        for (l = 0; l < Nx_rank; l++)
          for (i = 0; i < nspec; i++)
            //#pragma omp parallel for private(j)
            for (j = 0; j < Nv * Nv * Nv; j++)
              f[l + order][i][j] =
                  0.5 * (f[l + order][i][j] + f_tmp[l + order][i][j]) +
                  0.25 * f_conv[l + order][i][j];

        // Last strang step - V advection with timestep dt/2 (combine?)

        // Recalc Poiss

        if (ecouple == 1) { // electrons only in background
          if (Te_start != Te_ref)
            get_ramp_Te(Te_arr, Nx_rank, Te_start, Te_ref, t, tfinal);
          else
            get_uniform_Te(Te_arr, Nx_rank,
                           Te_ref); // fixed background temperature
        } else
          Te_arr = T0_oned;

        if (ecouple == 2) {
          for (l = 0; l < Nx_rank; l++) {
            source[l] = 0.0;
            for (i = 0; i < nspec; i++)
              n_oned[l][i] = getDensity(f[l + order][i], i);

            for (i = 1; i < nspec; i++)
              source[l] +=
                  Z_oned[l][i] *
                  n_oned[l][i]; // total number of free electrons in each cell
            source[l] -= n_oned[l][0];
          }
        } else {
          for (l = 0; l < Nx_rank; l++) {
            for (i = 0; i < nspec; i++)
              n_oned[l][i] = getDensity(f[l + order][i], i);
            zBarFunc2(nspec, Te_arr[l], Z_max, n_oned[l], Z_oned[l]);
            source[l] = 0.0;
            for (i = 0; i < nspec; i++)
              source[l] +=
                  Z_oned[l][i] *
                  n_oned[l][i]; // total number of free electrons in each cell
          }
        }

        if (ecouple == 2) {
          simplePoisson(Nx_rank, source, dx, Lx, PoisPot);
        } else {
          if (poissFlavor == 0) { // no E-field
            for (l = 0; l < Nx_rank; l++)
              PoisPot[l] = 0.0;
          } else if (poissFlavor == 11) // Linear Yukawa
            PoissLinPeriodic1D(Nx_rank, source, dx, Lx, PoisPot, Te_arr);
          else if (poissFlavor == 12) // Nonlinear Yukawa
            PoissNonlinPeriodic1D(Nx_rank, source, dx, Lx, PoisPot, Te_arr);
          else if (poissFlavor == 21) // Linear Thomas-Fermi
            PoissLinPeriodic1D_TF(Nx_rank, source, dx, Lx, PoisPot, Te_arr);
          else if (poissFlavor == 22) // Nonlinear Thomas-Fermi
            PoissNonlinPeriodic1D_TF(Nx_rank, source, dx, Lx, PoisPot, Te_arr);
        }

        // RK2 Step 1
        for (i = 0; i < nspec; i++) {
          advectTwo_v(f, f_conv, PoisPot, Z_oned, m[i], i);
        }

        for (l = 0; l < Nx_rank; l++)
          for (i = 0; i < nspec; i++)
            //#pragma omp parallel for private(j)
            for (j = 0; j < Nv * Nv * Nv; j++)
              f_tmp[l + order][i][j] =
                  f[l + order][i][j] + 0.5 * f_conv[l + order][i][j];

        // Recalc Poiss

        if (ecouple == 1) { // electrons only in background
          if (Te_start != Te_ref)
            get_ramp_Te(Te_arr, Nx_rank, Te_start, Te_ref, t, tfinal);
          else
            get_uniform_Te(Te_arr, Nx_rank,
                           Te_ref); // fixed background temperature
        } else
          Te_arr = T0_oned;

        if (ecouple == 2) {
          for (l = 0; l < Nx_rank; l++) {
            source[l] = 0.0;

            for (i = 0; i < nspec; i++)
              n_oned[l][i] = getDensity(f_tmp[l + order][i], i);
            for (i = 0; i < nspec; i++)
              source[l] +=
                  Z_oned[l][i] *
                  n_oned[l][i]; // total number of free electrons in each cell
            source[l] -= n_oned[l][0];
          }

        } else {
          for (l = 0; l < Nx_rank; l++) {
            for (i = 0; i < nspec; i++)
              n_oned[l][i] = getDensity(f_tmp[l + 1][i], i);
            zBarFunc2(nspec, Te_arr[l], Z_max, n_oned[l], Z_oned[l]);
            source[l] = 0.0;
            for (i = 0; i < nspec; i++)
              source[l] +=
                  Z_oned[l][i] *
                  n_oned[l][i]; // total number of free electrons in each cell
          }
        }

        if (ecouple == 2) {
          simplePoisson(Nx_rank, source, dx, Lx, PoisPot);
        } else {
          if (poissFlavor == 0) { // no E-field
            for (l = 0; l < Nx_rank; l++)
              PoisPot[l] = 0.0;
          } else if (poissFlavor == 11) // Linear Yukawa
            PoissLinPeriodic1D(Nx_rank, source, dx, Lx, PoisPot, Te_arr);
          else if (poissFlavor == 12) // Nonlinear Yukawa
            PoissNonlinPeriodic1D(Nx_rank, source, dx, Lx, PoisPot, Te_arr);
          else if (poissFlavor == 21) // Linear Thomas-Fermi
            PoissLinPeriodic1D_TF(Nx_rank, source, dx, Lx, PoisPot, Te_arr);
          else if (poissFlavor == 22) // Nonlinear Thomas-Fermi
            PoissNonlinPeriodic1D_TF(Nx_rank, source, dx, Lx, PoisPot, Te_arr);
        }

        // RK2 Step 2
        for (i = 0; i < nspec; i++) {
          advectTwo_v(f_tmp, f_conv, PoisPot, Z_oned, m[i], i);
        }
        for (l = 0; l < Nx_rank; l++)
          for (i = 0; i < nspec; i++)
            //#pragma omp parallel for private(j)
            for (j = 0; j < Nv * Nv * Nv; j++)
              f[l + order][i][j] =
                  0.5 * (f[l + order][i][j] + f_tmp[l + order][i][j]) +
                  0.25 * f_conv[l + order][i][j];
      }
    }
    t += dt;
    nT++;
  }

  if ((dims == 0) && (restartFlag > 0))
    store_distributions_homog(f_zerod, t, input_filename);

  // clean up

  // universal to both
  for (i = 0; i < nspec; i++) {
    free(c[i]);
    free(wts[i]);
  }
  free(c);
  free(wts);
  free(vref);
  free(Lv);

  for (i = 0; i < nspec; i++) {
    fclose(outputFileDens[i]);
    fclose(outputFileVelo[i]);
    fclose(outputFileTemp[i]);
  }
  free(outputFileDens);
  free(outputFileVelo);
  free(outputFileTemp);

  if (dims == 0) {
    for (i = 0; i < nspec; i++) {
      free(f_zerod[i]);
    }
    free(f_zerod);
  }

  if (dims == 1) {
    free(x);
    free(dxarray);

    fclose(outputFile_x);
    fclose(outputFilePoiss);

    for (l = 0; l < Nx_rank; l++) {
      free(n_oned[l]);
      free(T_oned[l]);
      free(v0_oned[l]);
      for (i = 0; i < nspec; i++) {
        free(v_oned[l][i]);
        free(f[l][i]);
        free(f_tmp[l][i]);
        free(f_conv[l][i]);
      }
      free(v_oned[l]);
      free(Z_oned[l]);

      free(f[l]);
      free(f_tmp[l]);
      free(f_conv[l]);
    }

    free(f);
    free(f_tmp);
    free(f_conv);
    free(n_oned);
    free(v_oned);
    free(v0_oned);
    free(T_oned);
    free(T0_oned);
    free(T_max);
    free(PoisPot);
    free(source);
    free(T_for_zbar);
    free(Z_oned);
    free(H_spec);
    free(H_spec_prev);
  }

  return 0;
}
