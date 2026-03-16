// BGK Code
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
#include "LB.h"
#include "TNB.h"
#include "fluid_electron.h"

int main(int argc, char **argv) {

  MPI_Init(&argc, &argv);

  int rank, numRanks;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
  MPI_Status status;
  int rankCounter;
  int rankOffset;
  double *momentBuffer, *momentBuffer2;
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

  // Flags for BGK collision rates
  int ecouple, CL_type, ion_type, MT_or_TR;
  int BGK_type;
  int TNBFlag;
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
  double *Te_arr, *Te_arr_allranks;
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
  double ***T_dir_oned;
  double *T0_oned;
  double **T_for_zbar;
  double *T_max, T0_max;

  double ***PTens_oned;

  double ***marginal_oned;

  // Velocity grid setup - 0D and 1D
  int Nv;
  int discret;    // 0 - uniform, 1 - gauss
  double *Lv;     // semi-length, domain is [-Lv,Lv]
  double v_sigma; // width of velocity domains in terms of thermal speeds
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

  double *PoisPot, *PoisPot_allranks;
  double *source, *source_buf, *source_allranks;
  int poissFlavor;
  int ionFix;

  // Time setup
  double dt;
  double tfinal;
  double t;
  int nT;
  int im_ex;
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
             &ndens_int, &velo_int, &T_int, &ecouple, &ionFix, &Te_start,
             &Te_ref, &CL_type, &ion_type, &MT_or_TR, &TNBFlag, &n_zerod,
             &v_val, &T_zerod, &dataFreq, &outputDist, &RHS_tol, &BGK_type,
             &beta, &hydro_kinscheme_flag, &input_file_data_flag,
             input_file_data_filename, input_filename);

  double ESBGK_nu = 0.5;
  if (rank == 0) {
    printf("NOTE: ESBGK_nu is hardcoded to 0.5 \n");
  }

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
  char phi_path[100];
  char x_path[100];
  char name_tmp[100];

  FILE **outputFileDens = malloc(nspec * sizeof(FILE *));
  FILE **outputFileVelo = malloc(nspec * sizeof(FILE *));
  FILE **outputFileTemp = malloc(nspec * sizeof(FILE *));
  FILE **outputFileDirTemp = malloc(nspec * sizeof(FILE *));
  FILE **outputFilePTens = malloc(nspec * sizeof(FILE *));
  FILE *outputFileH;
  FILE *outputFileBGK;
  FILE *outputFileTime;
  FILE *outputFile_x;
  FILE *outputFilePoiss;
  FILE *outputFilePhi;
  FILE **outputFileMarginal = malloc(nspec * sizeof(FILE *));

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

      strcpy(temp_path, output_path);
      sprintf(name_tmp, "_dirtemp%d", i);
      strcat(temp_path, name_tmp);
      outputFileDirTemp[i] = fopen(temp_path, "w");

      strcpy(temp_path, output_path);
      sprintf(name_tmp, "_PTens%d", i);
      strcat(temp_path, name_tmp);
      outputFilePTens[i] = fopen(temp_path, "w");

      strcpy(temp_path, output_path);
      sprintf(name_tmp, "_marginal%d", i);
      strcat(temp_path, name_tmp);
      if ((restartFlag == 2) || (restartFlag == 4))
        outputFileMarginal[i] = fopen(temp_path, "a");
      else
        outputFileMarginal[i] = fopen(temp_path, "w");
    }

    strcpy(temp_path, output_path);
    sprintf(name_tmp, "_tempe");
  }

  // Set up file to store the times
  strcpy(time_path, output_path);
  sprintf(name_tmp, "_time");
  strcat(time_path, name_tmp);

  strcpy(H_path, output_path);
  sprintf(name_tmp, "_H");
  strcat(H_path, name_tmp);

  strcpy(BGK_path, output_path);
  sprintf(name_tmp, "_BGK.%d", rank);
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
    strcpy(phi_path, output_path);
    strcat(phi_path, "_phi");

    if ((restartFlag == 2) || (restartFlag == 4)) {
      outputFilePoiss = fopen(poiss_path, "a");
      outputFilePhi = fopen(phi_path, "a");
    } else {
      outputFilePoiss = fopen(poiss_path, "w");
      outputFilePhi = fopen(phi_path, "w");
    }
    // store physical mesh
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
          Lv, &t, &nT,
          input_filename); // note - this assumes uniform grid for now
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

    if (restartFlag == 2) {
      load_time_inhomog(input_filename, &nT, &t);
      printf("Initial timestep and time %d %lg\n", nT, t);
    }

    // Physical grid allocation and initialization
    make_mesh(Nx, Lx, order, &Nx_rank, &Nx_ranks, &x, &dxarray);
    momentBuffer = malloc(15 * (Nx_rank + 1) * sizeof(double));
    momentBuffer2 = malloc(Nv * (Nx_rank + 1) * sizeof(double));

    if (input_file_data_flag) {
      intervalLimits =
          malloc(sizeof(double)); // we only care about the leftmost guy
      intervalLimits[0] = x[order];
    }

    dx = Lx / Nx;

    // Write the full mesh for storage purposes
    if (rank == 0) {
      for (l = 0; l < Nx; l++)
        fprintf(outputFile_x, "%+le ", (l + 0.5) * dx - 0.5 * Lx);

      fclose(outputFile_x);
    }

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
    T_dir_oned = malloc(Nx_rank * sizeof(double **));
    PTens_oned = malloc(Nx_rank * sizeof(double **));

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
      PTens_oned[l] = malloc(nspec * sizeof(double *));
      T_dir_oned[l] = malloc(nspec * sizeof(double *));
      for (i = 0; i < nspec; i++) {
        v_oned[l][i] = malloc(3 * sizeof(double));
        T_dir_oned[l][i] = malloc(3 * sizeof(double));
        PTens_oned[l][i] = malloc(9 * sizeof(double));
      }
      T_oned[l] = malloc(nspec * sizeof(double));

      v0_oned[l] = malloc(3 * sizeof(double));

      v0_oned[l][0] = 0.0;
      v0_oned[l][1] = 0.0;
      v0_oned[l][2] = 0.0;

      T0_oned[l] = 0.0;
    }

    printf("setting up marginals\n");
    if (outputDist == 2) {
      marginal_oned = malloc(nspec * sizeof(double **));
      for (s = 0; s < nspec; s++) {
        marginal_oned[s] = malloc(Nx_rank * sizeof(double *));
        for (j = 0; j < Nx_rank; j++)
          marginal_oned[s][j] = malloc(Nv * sizeof(double));
      }
    }
    printf("done setting up marginals\n");

    if (input_file_data_flag) {
      if (rank == 0) {
        printf("input_file_data_flag: %d, filename: %s\n", input_file_data_flag,
               input_file_data_filename);
      }

      initialize_sol_load_inhom_file(rank, Nx_ranks, Nx, nspec, n_oned, v_oned,
                                     T_oned, input_file_data_filename);

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

    printf("Setting up input done\n");

    // Allocate for poisson calc
    PoisPot = malloc((Nx_rank + 2 * order) * sizeof(double));
    PoisPot_allranks = malloc(Nx * sizeof(double));
    source = malloc(Nx_rank * sizeof(double));
    source_buf = malloc(2 * (Nx_rank + 1) * sizeof(double));
    source_allranks = malloc(Nx * sizeof(double));
    T_for_zbar = malloc(Nx_rank * sizeof(double));
    Te_arr = malloc(Nx_rank * sizeof(double));
    Te_arr_allranks = malloc(Nx * sizeof(double));

    // Initialize electron temperature for fluid electron case
    if (Te_ref == 0) {
      for (int Te_i = 0; Te_i < Nx_rank; Te_i++)
        Te_arr[Te_i] = Te_start;
    }

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

    double vmax = 0.0;

    for (i = 0; i < nspec; i++) {
      if (T0_max > T_max[i])
        vref[i] = sqrt((T0_max / ERG_TO_EV_CGS) / m[i]);
      else
        vref[i] = sqrt((T_max[i] / ERG_TO_EV_CGS) / m[i]);

      if (vmax < vref[i])
        vmax = vref[i];

      Lv[i] = v_sigma * vref[i];
      if (rank == 0)
        printf("Lv for species %d: %g\n", i, Lv[i]);
      c[i] = malloc(Nv * sizeof(double));
      wts[i] = malloc(Nv * sizeof(double));
    }

    double CFL = dt * vmax / dx;

    if (rank == 0) {
      printf("Streaming CFL: %g\n", CFL);

      if ((CFL > 1.0) || ((order == 2) && (CFL > 0.5))) {
        printf("\n\n******************************\n");
        printf("STREAMING COURANT CONDITION NOT SATISFIED, MAY BE UNSTABLE\n");
        printf("******************************\n\n\n");
      }
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

    io_init_inhomog(Nx_rank, Nv, nspec, order, c, m);
    if (outputDist == 1) {
      store_grid_inhomog(input_filename, rank);
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

    // if(rank == 0)
    //   printf("Warning: Maxwell BCs hacked in\n");
    // double n_bc[2] = {1.0622425812505312e-10, 6.547270642496101e-12};
    // double v_bc[2] = {0.0, 0.0};
    // double T_bc[2] = {8.61e-5, 8.61e-5};
    // initialize_transport_max_bc(n_bc, v_bc, T_bc, m[0]);

    initialize_transport(Nv, Nx_rank, nspec, x, dxarray, Lx, c, order, dt,
                         poissFlavor);

    // check if we are loading moment data from a file
    // Note: This is not yet MPI-ified
    if (input_file_data_flag) {
      initialize_sol_inhom_file(f, Nx_rank, nspec, Nv, order, c, m, n_oned,
                                v_oned, T_oned);
    } else
      initialize_sol_inhom(f, numint, intervalLimits, ndens_int, velo_int,
                           T_int, Nx_rank, x, nspec, Nv, order, c, m, n_oned,
                           v_oned, T_oned);

    if (restartFlag == 2) {
      // initialize distribution function with something non-physical
      for (l = 0; l < Nx_rank + 2 * order; l++) {
        for (s = 0; s < nspec; s++) {
          for (i = 0; i < Nv * Nv * Nv; i++) {
            f[l][s][i] = -1.0;
          }
        }
      }

      // Now load the distros
      load_distributions_inhomog(f, input_filename, nT, rank);
    }

    if (rank == 0) {
      printf("Initial condition setup complete\n");
      fflush(stdout);
    }
  }

  initialize_moments(Nv, nspec, c, wts);
  initialize_BGK(nspec, Nv, m, c, order, ecouple, CL_type, ion_type, MT_or_TR,
                 tauFlag, TNBFlag, input_filename);
  initialize_LB(nspec, Nv, m, c);
  initialize_fluid_electron(Nx_rank, nspec, dt);

  if (TNBFlag)
    initializeTNB(Nv, c, wts);

  if (!((restartFlag == 2) || (restartFlag == 4))) {
    t = 0.0;
    nT = 0;
  }
  if (restartFlag > 2) {
    // Check to see if we need to store initial BGK error data
    if (rank == 0) {
      printf("Setting initial BGK norm\n");
    }
    zBarFunc2(nspec, T0, Z_max, n_zerod, Z_zerod);
    BGK_norm(f_zerod, BGK_f_minus_eq_init, Z_zerod, dt, T0);
  }

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
      if ((ecouple != 2) && (ionFix != 1))
        zBarFunc2(nspec, T0, Z_max, n_zerod, Z_zerod);
      else
        Z_zerod = Z_max;

      // Output to files

      BGK_norm(f_zerod, BGK_f_minus_eq, Z_zerod, dt, T0);

      if (outcount == dataFreq || (t == 0)) {
        fprintf(outputFileTime, "%e\n", t);
        for (i = 0; i < nspec; i++) {
          fprintf(outputFileDens[i], "%e ", n_zerod[i]);
          fprintf(outputFileVelo[i], "%e ", v_zerod[i][0]);
          fprintf(outputFileTemp[i], "%e ", T_zerod[i]);
          fprintf(outputFileH, "%e ", (H_spec[i] - H_spec_prev[i]) / dt);
          fprintf(outputFileDens[i], "\n");
          fprintf(outputFileVelo[i], "\n");
          fprintf(outputFileTemp[i], "\n");
          outcount = 0;
        }
        fprintf(outputFileH, "%e\n", (Htot - Htot_prev) / dt);

        if (outputDist == 1)
          store_distributions_homog(f_zerod, t, nT, input_filename);

        for (i = 0; i < nspec; i++)
          for (j = 0; j < nspec; j++)
            fprintf(outputFileBGK, "%e,", BGK_f_minus_eq[i][j]);
        fprintf(outputFileBGK, "\n");
      }
      outcount += 1;

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
      } else if (im_ex == 1) {
        BGK_im_linear(f_zerod, f_zerod_tmp, Z_zerod, dt, T0);
        for (i = 0; i < nspec; i++)
          for (j = 0; j < Nv * Nv * Nv; j++)
            f_zerod[i][j] = f_zerod_tmp[i][j];
      } else if (im_ex == 2) {
        /*
        BGK_im_nonlinear();
        for (i = 0; i < nspec; i++)
          for (j = 0; j < Nv * Nv * Nv; j++)
            f_zerod[i][j] = f_zerod_tmp[i][j];
        */
      } else {
        printf("Error - please set im_ex = 0 (explicit), 1 (linear implicit), "
               "or 2 (nonlinear implicit) in your input file.\n");
        exit(1);
      }
    } else if (dims == 1) {

      Htot = 0.0;
      // Calculate moment data in all cells
      for (l = 0; l < Nx_rank; l++) {

        ntot = 0.0;
        rhotot = 0.0;
        for (i = 0; i < nspec; i++) {
          n_oned[l][i] = getDensity(f[l + order][i], i);
          ntot += n_oned[l][i];
          rhotot += m[i] * n_oned[l][i];
          getBulkVel(f[l + order][i], v_oned[l][i], n_oned[l][i], i);
          T_oned[l][i] =
              getTemp(m[i], n_oned[l][i], v_oned[l][i], f[l + order][i], i);
          for (int idir = 0; idir < 3; ++idir) {
            T_dir_oned[l][i][idir] = getDirTemp(
                m[i], n_oned[l][i], v_oned[l][i], f[l + order][i], i, idir);
          }
        }
	
        // calculate cross moment data
        if (BGK_type == 22) {
          for (int sp = 0; sp < nspec; sp++)
	    getStressTensor(m[sp], n_oned[l][sp], v_oned[l][sp],
                            f[l + order][sp], PTens_oned[l][sp], sp);
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
        for (k = 0; k < numRanks; k++) {
          if (k == rank) {
            for (i = 0; i < nspec; i++) {
              if (isnan(n_oned[l][i])) {
                printf("Something weird is going on: What did Jan Say? The "
                       "Michael Scott "
                       "Story. By Michael Scott. With Dwight Schrute.\n NaN "
                       "detected \n");
                for (j = 0; j < nspec; j++) {
                  printf("rank %d x %d i %d j %d n: %g v: %g T: %g Z: %g Te: "
                         "%g\n",
                         rank, l, i, j, n_oned[l][j], v0_oned[l][j],
                         T_oned[l][j], Z_oned[l][j], Te_arr[l]);
                }
                exit(1);
              }
            }
          }
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

      if (ecouple == 1) {    // electrons only in background
        if (Te_start == 0) { // cubic case
          get_ramp_Te_cubic(Te_arr, Nx_rank, 7.481, 1e-9, t);
        } else if (Te_ref == 0) { // coupled electron case, note only works for
                                  // 1 spec + ecouple 1
          // no-op, will update this later
          // fluid_electron_coll_only(n_oned, Te_arr, T_oned);
        } else if (Te_start != Te_ref)
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
            source[l] += Z_oned[l][i] * n_oned[l][i];
          }
          source[l] -= n_oned[l][0]; // electrons
        }
      } else {

        // Calculate ionization data
        for (l = 0; l < Nx_rank; l++) {
          // Get ionization

          if ((ecouple != 2) && (ionFix != 1))
            zBarFunc2(nspec, Te_arr[l], Z_max, n_oned[l], Z_oned[l]);
          else
            Z_oned[l] = Z_max;

          // Set up local Poisson RHS vector
          source[l] = 0.0;
          for (i = 0; i < nspec; i++) {
            source[l] +=
                Z_oned[l][i] *
                n_oned[l][i]; // total number of free electrons in each cell
          }
        }
      }

      // Set up the source/RHS array for the Poisson solve

      if (rank == 0) {

        for (l = 0; l < Nx_rank; l++) {
          source_allranks[l] = source[l];
          Te_arr_allranks[l] = Te_arr[l];
        }

        if (numRanks > 1) {
          rankOffset = Nx_rank;

          // Get the source/RHS value and electron temperatures from all the
          // other ranks.
          for (rankCounter = 1; rankCounter < numRanks; rankCounter++) {
            MPI_Recv(source_buf, 2 * Nx_ranks[rankCounter], MPI_DOUBLE,
                     rankCounter, 200 + rankCounter, MPI_COMM_WORLD, &status);
            for (l = 0; l < Nx_ranks[rankCounter]; l++) {
              source_allranks[l + rankOffset] = source_buf[0 + 2 * l];
              Te_arr_allranks[l + rankOffset] = source_buf[1 + 2 * l];
            }
            rankOffset += Nx_ranks[rankCounter];
          }
        }
      } else {
        // Send the RHS value to rank 0
        for (l = 0; l < Nx_rank; l++) {
          source_buf[0 + 2 * l] = source[l];
          source_buf[1 + 2 * l] = Te_arr[l];
        }
        MPI_Send(source_buf, 2 * Nx_rank, MPI_DOUBLE, 0, 200 + rank,
                 MPI_COMM_WORLD);
      }

      // Wait until all source stuff is sent
      MPI_Barrier(MPI_COMM_WORLD);

      if (rank == 0) { // Rank 0 performs the Poisson Solve

        if (ecouple == 2) {
          simplePoisson(Nx, source_allranks, dx, Lx, PoisPot_allranks);
        } else {

          if (poissFlavor == 0) { // no E-field
            for (l = 0; l < Nx; l++)
              PoisPot_allranks[l] = 0.0;
          } else if (poissFlavor == 11) // Linear Yukawa
            PoissLinPeriodic1D(Nx, source_allranks, dx, Lx, PoisPot_allranks,
                               Te_arr_allranks);
          else if (poissFlavor == 12) // Nonlinear Yukawa
            PoissNonlinPeriodic1D(Nx, source_allranks, dx, Lx, PoisPot_allranks,
                                  Te_arr_allranks);
          else if (poissFlavor == 21) // Linear Thomas-Fermi
            PoissLinPeriodic1D_TF(Nx, source_allranks, dx, Lx, PoisPot_allranks,
                                  Te_arr_allranks);
          else if (poissFlavor == 22) // Nonlinear Thomas-Fermi
            PoissNonlinPeriodic1D_TF(Nx, source_allranks, dx, Lx,
                                     PoisPot_allranks, Te_arr_allranks);
          else if (poissFlavor == 42) // ohm law
            fluid_electron_get_eE(Nx, source_allranks, dx, PoisPot_allranks,
                                  Te_arr_allranks);
        }
      }

      // Distribute back to the other ranks
      if (rank == 0) {
        rankOffset = 0;

        // SET LOCAL GHOST CELLS
        if (order == 1) {
          PoisPot[0] = PoisPot_allranks[0]; // MODIFIED
          if (numRanks == 1)
            PoisPot[Nx_rank + 1] = PoisPot_allranks[0];
          else
            PoisPot[Nx_rank + 1] = PoisPot_allranks[Nx_rank + 1];
        } else if (order == 2) {
          PoisPot[0] = PoisPot_allranks[Nx - 2];
          PoisPot[1] = PoisPot_allranks[Nx - 1];
          if (numRanks == 1) {
            PoisPot[Nx_rank + 2] = PoisPot_allranks[0];
            PoisPot[Nx_rank + 3] = PoisPot_allranks[1];
          } else {
            PoisPot[Nx_rank + 2] = PoisPot_allranks[Nx_rank];
            PoisPot[Nx_rank + 3] = PoisPot_allranks[Nx_rank + 1];
          }
        }

        // Set main body of PoisPot
        for (l = 0; l < Nx_rank; l++) {
          PoisPot[l + order] = PoisPot_allranks[l];
        }

        if (numRanks > 1) {

          rankOffset = Nx_rank;

          for (rankCounter = 1; rankCounter < numRanks - 1; rankCounter++) {

            // Send local ghost cell info
            if (order == 1) {
              source_buf[0] = PoisPot_allranks[rankOffset - 1];
              source_buf[Nx_ranks[rankCounter] + 1] =
                  PoisPot_allranks[rankOffset + Nx_ranks[rankCounter]];
            } else if (order == 2) {
              source_buf[0] = PoisPot_allranks[rankOffset - 2];
              source_buf[1] = PoisPot_allranks[rankOffset - 1];
              source_buf[Nx_ranks[rankCounter] + 2] =
                  PoisPot_allranks[rankOffset + Nx_ranks[rankCounter]];
              source_buf[Nx_ranks[rankCounter] + 3] =
                  PoisPot_allranks[rankOffset + Nx_ranks[rankCounter] + 1];
            }

            // Set main body of PoisPot on rankCounter
            for (l = 0; l < Nx_ranks[rankCounter]; l++) {
              source_buf[l + order] = PoisPot_allranks[rankOffset + l];
            }

            MPI_Send(source_buf, Nx_ranks[rankCounter] + 2 * order, MPI_DOUBLE,
                     rankCounter, rankCounter, MPI_COMM_WORLD);

            rankOffset += Nx_ranks[rankCounter];
          }

          // Deal with periodic BC for rightmost rank

          // SET LOCAL GHOST CELLS
          if (order == 1) {
            source_buf[0] = PoisPot_allranks[rankOffset - 1];
            source_buf[Nx_ranks[numRanks - 1] + 1] = PoisPot_allranks[Nx - 1];
          } else if (order == 2) {
            source_buf[0] = PoisPot_allranks[rankOffset - 2];
            source_buf[1] = PoisPot_allranks[rankOffset - 1];
            source_buf[Nx_ranks[numRanks - 1] + 2] = PoisPot_allranks[0];
            source_buf[Nx_ranks[numRanks - 1] + 3] = PoisPot_allranks[1];
          }

          for (l = 0; l < Nx_ranks[numRanks - 1]; l++) {
            source_buf[l + order] = PoisPot_allranks[rankOffset + l];
          }

          MPI_Send(source_buf, Nx_ranks[numRanks - 1] + 2 * order, MPI_DOUBLE,
                   numRanks - 1, numRanks - 1, MPI_COMM_WORLD);
        }
      } else {
        // Get potential from rank 0

        MPI_Recv(source_buf, Nx_rank + 2 * order, MPI_DOUBLE, 0, rank,
                 MPI_COMM_WORLD, &status);

        // Set Poispot
        for (l = 0; l < Nx_rank + 2 * order; l++) {
          PoisPot[l] = source_buf[l];
        }
      }

      // Moments and initial electric field calculated - save if needed

      if (rank == 0) {
        printf("Calculating moments for output\n");
	fflush(stdout);
      }

      // MPI-ified output.
      if (outcount == dataFreq) {
        if (rank == 0) {
          for (l = 0; l < Nx_rank; l++) {
            for (i = 0; i < nspec; i++) {
              fprintf(outputFileDens[i], "%e ", n_oned[l][i]);
              fprintf(outputFileVelo[i], "%e ", v_oned[l][i][0]);
              fprintf(outputFileTemp[i], "%e ", T_oned[l][i]);
              fprintf(outputFileDirTemp[i], "%e %e %e ", T_dir_oned[l][i][0],
                      T_dir_oned[l][i][1], T_dir_oned[l][i][2]);
              fprintf(outputFilePTens[i], "%e %e %e %e %e %e %e %e %e ",
                      PTens_oned[l][i][0], PTens_oned[l][i][1],
                      PTens_oned[l][i][2], PTens_oned[l][i][3],
                      PTens_oned[l][i][4], PTens_oned[l][i][5],
                      PTens_oned[l][i][6], PTens_oned[l][i][7],
                      PTens_oned[l][i][8]);	    
              if (outputDist == 2) {
                getMarginal(f[l + order][i], marginal_oned[i][l], i);
                for (j = 0; j < Nv; ++j) {
                  fprintf(outputFileMarginal[i], "%e ", marginal_oned[i][l][j]);
                }
                fprintf(outputFileMarginal[i], "\n");
              }
            }
          }

          // get from other ranks
          for (rankCounter = 1; rankCounter < numRanks; rankCounter++) {
            for (s = 0; s < nspec; s++) {
              MPI_Recv(momentBuffer, 15 * Nx_ranks[rankCounter], MPI_DOUBLE,
                       rankCounter, 100 + s, MPI_COMM_WORLD, &status);
              for (l = 0; l < Nx_ranks[rankCounter]; l++) {
                fprintf(outputFileDens[s], "%e ", momentBuffer[0 + 15 * l]);
                fprintf(outputFileVelo[s], "%e ", momentBuffer[1 + 15 * l]);
                fprintf(outputFileTemp[s], "%e ", momentBuffer[2 + 15 * l]);
                fprintf(outputFileDirTemp[s], "%e %e %e ",
                        momentBuffer[3 + 15 * l], momentBuffer[4 + 15 * l],
                        momentBuffer[5 + 15 * l]);

                fprintf(outputFilePTens[s], "%e %e %e %e %e %e %e %e %e",
                        momentBuffer[6 + 15 * l], momentBuffer[7 + 15 * l],
                        momentBuffer[8 + 15 * l], momentBuffer[9 + 15 * l],
                        momentBuffer[10 + 15 * l], momentBuffer[11 + 15 * l],
                        momentBuffer[12 + 15 * l], momentBuffer[13 + 15 * l],
                        momentBuffer[14 + 15 * l]);
              }
            }
          }

          if (outputDist == 2) {
            for (rankCounter = 1; rankCounter < numRanks; rankCounter++) {
              for (s = 0; s < nspec; s++) {
                MPI_Recv(momentBuffer2, Nv * Nx_ranks[rankCounter], MPI_DOUBLE,
                         rankCounter, 200 + s, MPI_COMM_WORLD, &status);
              }
              printf("setting up marginal output, rank %d \n", rank);
              for (l = 0; l < Nx_ranks[rankCounter]; l++) {
                for (j = 0; j < Nv; j++) {
                  fprintf(outputFileMarginal[s], "%e ",
                          momentBuffer2[j + Nv * l]);
                }
                fprintf(outputFileMarginal[s], "\n");
              }
            }
          }

          // Print poisson information - this produces eE in V / m
          if (poissFlavor != 42) {
            fprintf(outputFilePoiss, "%e ",
                    -50 * (PoisPot_allranks[1] - PoisPot_allranks[Nx - 1]) / dx *
                        ERG_TO_EV_CGS);
            fprintf(outputFilePhi, "%e ", Te_arr_allranks[0]);
            for (l = 1; l < Nx - 1; l++) {
              fprintf(outputFilePoiss, "%e ",
                      -50 * (PoisPot_allranks[l + 1] - PoisPot_allranks[l - 1]) / dx *
                          ERG_TO_EV_CGS);

              fprintf(outputFilePhi, "%e ", Te_arr_allranks[l]);
            }
            fprintf(outputFilePoiss, "%e ",
                    -50 * (PoisPot_allranks[0] - PoisPot_allranks[Nx - 2]) / dx *
                        ERG_TO_EV_CGS);
            fprintf(outputFilePhi, "%e ", Te_arr_allranks[Nx - 1]);
          } else {
            for (l = 0; l < Nx; l++) {
              fprintf(outputFilePoiss, "%e ",
                     50 *  PoisPot_allranks[l] * ERG_TO_EV_CGS);
              fprintf(outputFilePhi, "%e ", Te_arr_allranks[l]);
            }
          }
          fprintf(outputFilePoiss, "\n");
          fprintf(outputFilePhi, "\n");

          // Close out this timestep
          fprintf(outputFileTime, "%e\n", t);
          for (i = 0; i < nspec; i++) {
            fprintf(outputFileDens[i], "\n");
            fprintf(outputFileVelo[i], "\n");
            fprintf(outputFileTemp[i], "\n");
            fprintf(outputFileDirTemp[i], "\n");
            fprintf(outputFilePTens[i], "\n");
          }
          if (outputDist > 1)
            store_distributions_inhomog(f, input_filename, nT, t, rank);
          if (outputDist == 3) {
            char BGK_name_buffer[256];
            for (int s = 0; s < nspec; ++s) {
              sprintf(BGK_name_buffer, "Data/%s_spec%d.step%d_rank%d_BGK.dat",
                      input_filename, s, nT, rank);
              FILE *fid_store = fopen(BGK_name_buffer, "w");
              for (int xj = 0; xj < Nx_rank; ++xj) {
                fwrite(f_tmp[xj + order][s], sizeof(double), Nv * Nv * Nv,
                       fid_store);
              }
              fclose(fid_store);
            }
          }

          outcount = 0;
        } else { // send to rank 0 for output purposes
          for (s = 0; s < nspec; s++) {
            for (l = 0; l < Nx_rank; l++) {
              T_oned[l][s] =
                  getTemp(m[s], n_oned[l][s], v_oned[l][s], f[l + order][s], s);
              for (int idir = 0; idir < 3; ++idir)
                T_dir_oned[l][s][idir] = getDirTemp(
                    m[s], n_oned[l][s], v_oned[l][s], f[l + order][s], s, idir);

              momentBuffer[0 + 15 * l] = n_oned[l][s];
              momentBuffer[1 + 15 * l] = v_oned[l][s][0];
              momentBuffer[2 + 15 * l] = T_oned[l][s];
              momentBuffer[3 + 15 * l] = T_dir_oned[l][s][0];
              momentBuffer[4 + 15 * l] = T_dir_oned[l][s][1];
              momentBuffer[5 + 15 * l] = T_dir_oned[l][s][2];
              momentBuffer[6 + 15 * l] = PTens_oned[l][s][0];
              momentBuffer[7 + 15 * l] = PTens_oned[l][s][1];
              momentBuffer[8 + 15 * l] = PTens_oned[l][s][2];
              momentBuffer[9 + 15 * l] = PTens_oned[l][s][3];
              momentBuffer[10 + 15 * l] = PTens_oned[l][s][4];
              momentBuffer[11 + 15 * l] = PTens_oned[l][s][5];
              momentBuffer[12 + 15 * l] = PTens_oned[l][s][6];
              momentBuffer[13 + 15 * l] = PTens_oned[l][s][7];
              momentBuffer[14 + 15 * l] = PTens_oned[l][s][8];
            }
            MPI_Send(momentBuffer, 15 * Nx_rank, MPI_DOUBLE, 0, 100 + s,
                     MPI_COMM_WORLD);
            if (outputDist == 2) {
              for (s = 0; s < nspec; ++s) {
                for (l = 0; l < Nx_rank; ++l) {
                  getMarginal(f[l + order][s], marginal_oned[s][l], i);
                  for (j = 0; j < Nv; ++j) {
                    momentBuffer2[j + Nv * l] = marginal_oned[s][l][j];
                  }
                }
                MPI_Send(momentBuffer2, Nv * Nx_rank, MPI_DOUBLE, 0, 200 + s,
                         MPI_COMM_WORLD);
              }
            }
          }
          if (outputDist > 0)
            store_distributions_inhomog(f, input_filename, nT, t, rank);
          if (outputDist == 3) {
            char BGK_name_buffer[256];
            for (int s = 0; s < nspec; ++s) {
              sprintf(BGK_name_buffer, "Data/%s_spec%d.step%d_rank%d_BGK.dat",
                      input_filename, s, nT, rank);
              FILE *fid_store = fopen(BGK_name_buffer, "w");
              for (int xj = 0; xj < Nx_rank; ++xj) {
                fwrite(f_tmp[xj + order][s], sizeof(double), Nv * Nv * Nv,
                       fid_store);
              }
              fclose(fid_store);
            }
          }
          outcount = 0;
        }
      }
      outcount += 1;

      // IO done, advance to the actual solution...
      if(rank == 0) {
	printf("IO done, now solving\n");
	fflush(stdout);
      }
      
      MPI_Barrier(MPI_COMM_WORLD);

      if (order == 1) {
        // ADVECT

        for (i = 0; i < nspec; i++) {
          advectOne(f, PoisPot, Z_oned, m[i], i);
        }

        if (ecouple == 1 && Te_ref == 0) {
          fluid_electron_ohmic_term(n_oned, v_oned, Te_arr, PoisPot, Z_oned);
        }

        if (!(BGK_type < 0)) {
          // COLLIDE
          if (im_ex == 0) {
            for (l = 0; l < Nx_rank; l++) {
              if (BGK_type == 0)
                BGK_ex(f[l + order], f_tmp[l + order], Z_oned[l], dt,
                       Te_arr[l]);
              else if (BGK_type == 11)
                LB_ex(f[l + order], f_tmp[l + order], Z_oned[l], dt, Te_arr[l]);
              else if (BGK_type == 22)
                ESBGK_ex(f[l + order], f_tmp[l + order], Z_oned[l], dt,
                         Te_arr[l], ESBGK_nu);

              for (i = 0; i < nspec; i++)
                for (j = 0; j < Nv * Nv * Nv; j++)
                  f[l + order][i][j] += dt * f_tmp[l + order][i][j];
            }
          } else if (im_ex == 1) {
            for (l = 0; l < Nx_rank; l++) {
              BGK_im_linear(f[l + order], f_tmp[l + order], Z_oned[l], dt,
                            Te_arr[l]);
              for (i = 0; i < nspec; i++)
                for (j = 0; j < Nv * Nv * Nv; j++)
                  f[l + order][i][j] = f_tmp[l + order][i][j];
            }
          } else if (im_ex == 2) {
            printf("Error - nonlinear implicit solve not implemented for 1D\n");
            exit(1);
            /*
              BGK_im_nonlinear();
              for (i = 0; i < nspec; i++)
              for (j = 0; j < Nv * Nv * Nv; j++)
              f_zerod[i][j] = f_zerod_tmp[i][j];
            */
          } else {
            printf("Error - please set Imp_exp = 0 (explicit), 1 (linear "
                   "implicit), or 2 (nonlinear implicit) in your input "
                   "file.\n");
            exit(1);
          }
          if (ecouple == 1 && Te_ref == 0) {
            fluid_electron_coll_only(n_oned, Te_arr, T_oned);
          }
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
            // #pragma omp parallel for private(j)
            for (j = 0; j < Nv * Nv * Nv; j++)
              f_tmp[l + order][i][j] =
                  f[l + order][i][j] + 0.5 * f_conv[l + order][i][j];

        // Do second step of Poisson solve - TODO refactor this into its own
        // routine

        if (ecouple == 1) { // electrons only in background

          if (Te_start == 0) { // cubic case
            get_ramp_Te_cubic(Te_arr, Nx_rank, 7.481, 1e-9, t);
          } else if (Te_ref ==
                     0) { // coupled electron case, note only works for 2 spec
            fluid_electron_coll_only(n_oned, Te_arr, T_oned);
          } else if (Te_start != Te_ref)
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

            if ((ecouple != 2) && (ionFix != 1))
              zBarFunc2(nspec, Te_arr[l], Z_max, n_oned[l], Z_oned[l]);
            else
              Z_oned[l] = Z_max;

            source[l] = 0.0;
            for (i = 0; i < nspec; i++)
              source[l] +=
                  Z_oned[l][i] *
                  n_oned[l][i]; // total number of free electrons in each cell
          }
        }

        // Set up the source/RHS array for the Poisson solve

        if (rank == 0) {

          for (l = 0; l < Nx_rank; l++) {
            source_allranks[l] = source[l];
            Te_arr_allranks[l] = Te_arr[l];
          }

          if (numRanks > 1) {
            rankOffset = Nx_rank;

            // Get the source/RHS value and electron temperatures from all the
            // other ranks.
            for (rankCounter = 1; rankCounter < numRanks; rankCounter++) {
              MPI_Recv(source_buf, 2 * Nx_ranks[rankCounter], MPI_DOUBLE,
                       rankCounter, 200 + rankCounter, MPI_COMM_WORLD, &status);
              for (l = 0; l < Nx_ranks[rankCounter]; l++) {
                source_allranks[l + rankOffset] = source_buf[0 + 2 * l];
                Te_arr_allranks[l + rankOffset] = source_buf[1 + 2 * l];
              }
              rankOffset += Nx_ranks[rankCounter];
            }
          }
        } else {
          // Send the RHS value to rank 0
          for (l = 0; l < Nx_rank; l++) {
            source_buf[0 + 2 * l] = source[l];
            source_buf[1 + 2 * l] = Te_arr[l];
          }
          MPI_Send(source_buf, 2 * Nx_rank, MPI_DOUBLE, 0, 200 + rank,
                   MPI_COMM_WORLD);
        }

        // Wait until all source stuff is sent
        MPI_Barrier(MPI_COMM_WORLD);

        if (rank == 0) { // Rank 0 performs the Poisson Solve

          if (ecouple == 2) {
            simplePoisson(Nx_rank, source, dx, Lx, PoisPot);
          } else {

            if (poissFlavor == 0) { // no E-field
              for (l = 0; l < Nx; l++)
                PoisPot_allranks[l] = 0.0;
            } else if (poissFlavor == 11) // Linear Yukawa
              PoissLinPeriodic1D(Nx, source_allranks, dx, Lx, PoisPot_allranks,
                                 Te_arr_allranks);
            else if (poissFlavor == 12) // Nonlinear Yukawa
              PoissNonlinPeriodic1D(Nx, source_allranks, dx, Lx,
                                    PoisPot_allranks, Te_arr_allranks);
            else if (poissFlavor == 21) // Linear Thomas-Fermi
              PoissLinPeriodic1D_TF(Nx, source_allranks, dx, Lx,
                                    PoisPot_allranks, Te_arr_allranks);
            else if (poissFlavor == 22) // Nonlinear Thomas-Fermi
              PoissNonlinPeriodic1D_TF(Nx, source_allranks, dx, Lx,
                                       PoisPot_allranks, Te_arr_allranks);
            else if (poissFlavor == 42) // ohm law
              fluid_electron_get_eE(Nx, source_allranks, dx, PoisPot_allranks,
                                    Te_arr_allranks);
          }
        }

        // Distribute back to the other ranks
        if (rank == 0) {
          rankOffset = 0;

          // SET LOCAL GHOST CELLS
          if (order == 1) {
            PoisPot[0] = PoisPot_allranks[0];
            if (numRanks == 1)
              PoisPot[Nx_rank + 1] = PoisPot_allranks[0];
            else
              PoisPot[Nx_rank + 1] = PoisPot_allranks[Nx_rank + 1];
          } else if (order == 2) {
            PoisPot[0] = PoisPot_allranks[Nx - 2];
            PoisPot[1] = PoisPot_allranks[Nx - 1];
            if (numRanks == 1) {
              PoisPot[Nx_rank + 2] = PoisPot_allranks[0];
              PoisPot[Nx_rank + 3] = PoisPot_allranks[1];
            } else {
              PoisPot[Nx_rank + 2] = PoisPot_allranks[Nx_rank];
              PoisPot[Nx_rank + 3] = PoisPot_allranks[Nx_rank + 1];
            }
          }

          // Set main body of PoisPot
          for (l = 0; l < Nx_rank; l++) {
            PoisPot[l + order] = PoisPot_allranks[l];
          }

          if (numRanks > 1) {

            rankOffset = Nx_rank;

            for (rankCounter = 1; rankCounter < numRanks - 1; rankCounter++) {

              // Send local ghost cell info
              if (order == 1) {
                source_buf[0] = PoisPot_allranks[rankOffset - 1];
                source_buf[Nx_ranks[rankCounter] + 1] =
                    PoisPot_allranks[rankOffset + Nx_ranks[rankCounter]];
              } else if (order == 2) {
                source_buf[0] = PoisPot_allranks[rankOffset - 2];
                source_buf[1] = PoisPot_allranks[rankOffset - 1];
                source_buf[Nx_ranks[rankCounter] + 2] =
                    PoisPot_allranks[rankOffset + Nx_ranks[rankCounter]];
                source_buf[Nx_ranks[rankCounter] + 3] =
                    PoisPot_allranks[rankOffset + Nx_ranks[rankCounter] + 1];
              }

              // Set main body of PoisPot on rankCounter
              for (l = 0; l < Nx_ranks[rankCounter]; l++) {
                source_buf[l + order] = PoisPot_allranks[rankOffset + l];
              }

              MPI_Send(source_buf, Nx_ranks[rankCounter] + 2 * order,
                       MPI_DOUBLE, rankCounter, rankCounter, MPI_COMM_WORLD);

              rankOffset += Nx_ranks[rankCounter];
            }

            // Deal with periodic BC for rightmost rank

            // SET LOCAL GHOST CELLS
            if (order == 1) {
              source_buf[0] = PoisPot_allranks[rankOffset - 1];
              source_buf[Nx_ranks[numRanks - 1] + 1] =
                  PoisPot_allranks[Nx - 1]; // modified
            } else if (order == 2) {
              source_buf[0] = PoisPot_allranks[rankOffset - 2];
              source_buf[1] = PoisPot_allranks[rankOffset - 1];
              source_buf[Nx_ranks[numRanks - 1] + 2] = PoisPot_allranks[0];
              source_buf[Nx_ranks[numRanks - 1] + 3] = PoisPot_allranks[1];
            }

            for (l = 0; l < Nx_ranks[numRanks - 1]; l++) {
              source_buf[l + order] = PoisPot_allranks[rankOffset + l];
            }

            MPI_Send(source_buf, Nx_ranks[numRanks - 1] + 2 * order, MPI_DOUBLE,
                     numRanks - 1, numRanks - 1, MPI_COMM_WORLD);
          }
        } else {
          // Get potential from rank 0

          MPI_Recv(source_buf, Nx_rank + 2 * order, MPI_DOUBLE, 0, rank,
                   MPI_COMM_WORLD, &status);

          // Set Poispot
          for (l = 0; l < Nx_rank + 2 * order; l++) {
            PoisPot[l] = source_buf[l];
          }
        }

        // Finished with determining poisson solve, now advect

        // RK2 Step 2
        for (i = 0; i < nspec; i++) {
          advectTwo_v(f_tmp, f_conv, PoisPot, Z_oned, m[i], i);
        }

        for (l = 0; l < Nx_rank; l++)
          for (i = 0; i < nspec; i++)
            // #pragma omp parallel for private(j)
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
            // #pragma omp parallel for private(j)
            for (j = 0; j < Nv * Nv * Nv; j++)
              f_tmp[l + order][i][j] =
                  f[l + order][i][j] + 0.5 * f_conv[l + order][i][j];

        // RK2 Step 2
        for (i = 0; i < nspec; i++) {
          advectTwo_x(f_tmp, f_conv, i);
        }

        for (l = 0; l < Nx_rank; l++)
          for (i = 0; i < nspec; i++)
            // #pragma omp parallel for private(j)
            for (j = 0; j < Nv * Nv * Nv; j++)
              f[l + order][i][j] =
                  0.5 * (f[l + order][i][j] + f_tmp[l + order][i][j]) +
                  0.25 * f_conv[l + order][i][j];

        if (ecouple == 1 && Te_ref == 0) {
          fluid_electron_ohmic_term(n_oned, v_oned, Te_arr, PoisPot, Z_oned);
        }

        if (!(BGK_type == -1)) {
          // Next Strang step - RK2 for collision with timstep dt
          for (l = 0; l < Nx_rank; l++) {

            if (im_ex == 0) {
              // Step 1
              BGK_ex(f[l + order], f_conv[l + order], Z_oned[l], dt, Te_arr[l]);
              for (i = 0; i < nspec; i++)
                for (j = 0; j < Nv * Nv * Nv; j++)
                  f_tmp[l + order][i][j] =
                      f[l + order][i][j] + dt * f_conv[l + order][i][j];

              // Step 2
              BGK_ex(f_tmp[l + order], f_conv[l + order], Z_oned[l], dt,
                     Te_arr[l]);
              for (i = 0; i < nspec; i++)
                for (j = 0; j < Nv * Nv * Nv; j++)
                  f[l + order][i][j] =
                      0.5 * (f[l + order][i][j] + f_tmp[l + order][i][j]) +
                      0.5 * dt * f_conv[l + order][i][j];
            }

            else if (im_ex == 1) {
              printf("Error - implicit solve in 1D only implemented for first "
                     "oder solve \n");
            } else if (im_ex == 2) {
              /*
                BGK_im_nonlinear();
                for (i = 0; i < nspec; i++)
                for (j = 0; j < Nv * Nv * Nv; j++)
                f_zerod[i][j] = f_zerod_tmp[i][j];
              */
            } else {
              printf("Error - please set im_ex = 0 (explicit), 1 (linear "
                     "implicit), or 2 (nonlinear implicit) in your input "
                     "file.\n");
              exit(1);
            }
          }
        }

        // Next strang step - x advection with timestep dt/2

        // RK2 Step 1
        for (i = 0; i < nspec; i++) {
          advectTwo_x(f, f_conv, i);
        }

        for (l = 0; l < Nx_rank; l++)
          for (i = 0; i < nspec; i++)
            // #pragma omp parallel for private(j)
            for (j = 0; j < Nv * Nv * Nv; j++)
              f_tmp[l + order][i][j] =
                  f[l + order][i][j] + 0.5 * f_conv[l + order][i][j];

        // RK2 Step 2
        for (i = 0; i < nspec; i++) {
          advectTwo_x(f_tmp, f_conv, i);
        }
        for (l = 0; l < Nx_rank; l++)
          for (i = 0; i < nspec; i++)
            // #pragma omp parallel for private(j)
            for (j = 0; j < Nv * Nv * Nv; j++)
              f[l + order][i][j] =
                  0.5 * (f[l + order][i][j] + f_tmp[l + order][i][j]) +
                  0.25 * f_conv[l + order][i][j];

        // Last strang step - V advection with timestep dt/2 (combine?)

        // Recalc Poiss

        if (ecouple == 1) {    // electrons only in background
          if (Te_start == 0) { // cubic case
            get_ramp_Te_cubic(Te_arr, Nx_rank, 7.481, 1e-9, t);
          } else if (Te_ref ==
                     0) { // coupled electron case, note only works for 2 spec
            fluid_electron_coll_only(n_oned, Te_arr, T_oned);
          } else if (Te_start != Te_ref)
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

            if ((ecouple != 2) && (ionFix != 1))
              zBarFunc2(nspec, Te_arr[l], Z_max, n_oned[l], Z_oned[l]);
            else
              Z_oned[l] = Z_max;

            source[l] = 0.0;
            for (i = 0; i < nspec; i++)
              source[l] +=
                  Z_oned[l][i] *
                  n_oned[l][i]; // total number of free electrons in each cell
          }
        }

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

            if ((ecouple != 2) && (ionFix != 1))
              zBarFunc2(nspec, Te_arr[l], Z_max, n_oned[l], Z_oned[l]);
            else
              Z_oned[l] = Z_max;

            source[l] = 0.0;
            for (i = 0; i < nspec; i++)
              source[l] +=
                  Z_oned[l][i] *
                  n_oned[l][i]; // total number of free electrons in each cell
          }
        }

        // Set up the source/RHS array for the Poisson solve

        if (rank == 0) {

          for (l = 0; l < Nx_rank; l++) {
            source_allranks[l] = source[l];
            Te_arr_allranks[l] = Te_arr[l];
          }

          if (numRanks > 1) {
            rankOffset = Nx_rank;

            // Get the source/RHS value and electron temperatures from all the
            // other ranks.
            for (rankCounter = 1; rankCounter < numRanks; rankCounter++) {
              MPI_Recv(source_buf, 2 * Nx_ranks[rankCounter], MPI_DOUBLE,
                       rankCounter, 200 + rankCounter, MPI_COMM_WORLD, &status);
              for (l = 0; l < Nx_ranks[rankCounter]; l++) {
                source_allranks[l + rankOffset] = source_buf[0 + 2 * l];
                Te_arr_allranks[l + rankOffset] = source_buf[1 + 2 * l];
              }
              rankOffset += Nx_ranks[rankCounter];
            }
          }
        } else {
          // Send the RHS value to rank 0
          for (l = 0; l < Nx_rank; l++) {
            source_buf[0 + 2 * l] = source[l];
            source_buf[1 + 2 * l] = Te_arr[l];
          }
          MPI_Send(source_buf, 2 * Nx_rank, MPI_DOUBLE, 0, 200 + rank,
                   MPI_COMM_WORLD);
        }

        // Wait until all source stuff is sent
        MPI_Barrier(MPI_COMM_WORLD);

        if (rank == 0) { // Rank 0 performs the Poisson Solve

          if (ecouple == 2) {
            simplePoisson(Nx_rank, source, dx, Lx, PoisPot);
          } else {

            if (poissFlavor == 0) { // no E-field
              for (l = 0; l < Nx; l++)
                PoisPot_allranks[l] = 0.0;
            } else if (poissFlavor == 11) // Linear Yukawa
              PoissLinPeriodic1D(Nx, source_allranks, dx, Lx, PoisPot_allranks,
                                 Te_arr_allranks);
            else if (poissFlavor == 12) // Nonlinear Yukawa
              PoissNonlinPeriodic1D(Nx, source_allranks, dx, Lx,
                                    PoisPot_allranks, Te_arr_allranks);
            else if (poissFlavor == 21) // Linear Thomas-Fermi
              PoissLinPeriodic1D_TF(Nx, source_allranks, dx, Lx,
                                    PoisPot_allranks, Te_arr_allranks);
            else if (poissFlavor == 22) // Nonlinear Thomas-Fermi
              PoissNonlinPeriodic1D_TF(Nx, source_allranks, dx, Lx,
                                       PoisPot_allranks, Te_arr_allranks);
            else if (poissFlavor == 42) // ohm law
              fluid_electron_get_eE(Nx, source_allranks, dx, PoisPot_allranks,
                                    Te_arr_allranks);
          }
        }

        // Distribute back to the other ranks
        if (rank == 0) {
          rankOffset = 0;

          // SET LOCAL GHOST CELLS
          if (order == 1) {
            PoisPot[0] = PoisPot_allranks[0];
            if (numRanks == 1)
              PoisPot[Nx_rank + 1] = PoisPot_allranks[0];
            else
              PoisPot[Nx_rank + 1] = PoisPot_allranks[Nx_rank + 1];
          } else if (order == 2) {
            PoisPot[0] = PoisPot_allranks[Nx - 2];
            PoisPot[1] = PoisPot_allranks[Nx - 1];
            if (numRanks == 1) {
              PoisPot[Nx_rank + 2] = PoisPot_allranks[0];
              PoisPot[Nx_rank + 3] = PoisPot_allranks[1];
            } else {
              PoisPot[Nx_rank + 2] = PoisPot_allranks[Nx_rank];
              PoisPot[Nx_rank + 3] = PoisPot_allranks[Nx_rank + 1];
            }
          }

          // Set main body of PoisPot
          for (l = 0; l < Nx_rank; l++) {
            PoisPot[l + order] = PoisPot_allranks[l];
          }

          if (numRanks > 1) {

            rankOffset = Nx_rank;

            for (rankCounter = 1; rankCounter < numRanks - 1; rankCounter++) {

              // Send local ghost cell info
              if (order == 1) {
                source_buf[0] = PoisPot_allranks[rankOffset - 1];
                source_buf[Nx_ranks[rankCounter] + 1] =
                    PoisPot_allranks[rankOffset + Nx_ranks[rankCounter]];
              } else if (order == 2) {
                source_buf[0] = PoisPot_allranks[rankOffset - 2];
                source_buf[1] = PoisPot_allranks[rankOffset - 1];
                source_buf[Nx_ranks[rankCounter] + 2] =
                    PoisPot_allranks[rankOffset + Nx_ranks[rankCounter]];
                source_buf[Nx_ranks[rankCounter] + 3] =
                    PoisPot_allranks[rankOffset + Nx_ranks[rankCounter] + 1];
              }

              // Set main body of PoisPot on rankCounter
              for (l = 0; l < Nx_ranks[rankCounter]; l++) {
                source_buf[l + order] = PoisPot_allranks[rankOffset + l];
              }

              MPI_Send(source_buf, Nx_ranks[rankCounter] + 2 * order,
                       MPI_DOUBLE, rankCounter, rankCounter, MPI_COMM_WORLD);

              rankOffset += Nx_ranks[rankCounter];
            }

            // Deal with periodic BC for rightmost rank

            // SET LOCAL GHOST CELLS
            if (order == 1) {
              source_buf[0] = PoisPot_allranks[rankOffset - 1];
              source_buf[Nx_ranks[numRanks - 1] + 1] = PoisPot_allranks[Nx - 1];
            } else if (order == 2) {
              source_buf[0] = PoisPot_allranks[rankOffset - 2];
              source_buf[1] = PoisPot_allranks[rankOffset - 1];
              source_buf[Nx_ranks[numRanks - 1] + 2] = PoisPot_allranks[0];
              source_buf[Nx_ranks[numRanks - 1] + 3] = PoisPot_allranks[1];
            }

            for (l = 0; l < Nx_ranks[numRanks - 1]; l++) {
              source_buf[l + order] = PoisPot_allranks[rankOffset + l];
            }

            MPI_Send(source_buf, Nx_ranks[numRanks - 1] + 2 * order, MPI_DOUBLE,
                     numRanks - 1, numRanks - 1, MPI_COMM_WORLD);
          }
        } else {
          // Get potential from rank 0

          MPI_Recv(source_buf, Nx_rank + 2 * order, MPI_DOUBLE, 0, rank,
                   MPI_COMM_WORLD, &status);

          // Set Poispot
          for (l = 0; l < Nx_rank + 2 * order; l++) {
            PoisPot[l] = source_buf[l];
          }
        }

        // Finished with determining poisson solve, now advect

        // RK2 Step 1
        for (i = 0; i < nspec; i++) {
          advectTwo_v(f, f_conv, PoisPot, Z_oned, m[i], i);
        }

        for (l = 0; l < Nx_rank; l++)
          for (i = 0; i < nspec; i++)
            // #pragma omp parallel for private(j)
            for (j = 0; j < Nv * Nv * Nv; j++)
              f_tmp[l + order][i][j] =
                  f[l + order][i][j] + 0.5 * f_conv[l + order][i][j];

        // Recalc Poiss

        if (ecouple == 1) {    // electrons only in background
          if (Te_start == 0) { // cubic case
            get_ramp_Te_cubic(Te_arr, Nx_rank, 7.481, 1e-9, t);
          } else if (Te_ref ==
                     0) { // coupled electron case, note only works for 2 spec
            fluid_electron_coll_only(n_oned, Te_arr, T_oned);
          } else if (Te_start != Te_ref)
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

            if ((ecouple != 2) && (ionFix != 1))
              zBarFunc2(nspec, Te_arr[l], Z_max, n_oned[l], Z_oned[l]);
            else
              Z_oned[l] = Z_max;

            source[l] = 0.0;
            for (i = 0; i < nspec; i++)
              source[l] +=
                  Z_oned[l][i] *
                  n_oned[l][i]; // total number of free electrons in each cell
          }
        }

        // Set up the source/RHS array for the Poisson solve

        if (rank == 0) {

          for (l = 0; l < Nx_rank; l++) {
            source_allranks[l] = source[l];
            Te_arr_allranks[l] = Te_arr[l];
          }

          if (numRanks > 1) {
            rankOffset = Nx_rank;

            // Get the source/RHS value and electron temperatures from all
            // the other ranks.
            for (rankCounter = 1; rankCounter < numRanks; rankCounter++) {
              MPI_Recv(source_buf, 2 * Nx_ranks[rankCounter], MPI_DOUBLE,
                       rankCounter, 200 + rankCounter, MPI_COMM_WORLD, &status);
              for (l = 0; l < Nx_ranks[rankCounter]; l++) {
                source_allranks[l + rankOffset] = source_buf[0 + 2 * l];
                Te_arr_allranks[l + rankOffset] = source_buf[1 + 2 * l];
              }
              rankOffset += Nx_ranks[rankCounter];
            }
          }
        } else {
          // Send the RHS value to rank 0
          for (l = 0; l < Nx_rank; l++) {
            source_buf[0 + 2 * l] = source[l];
            source_buf[1 + 2 * l] = Te_arr[l];
          }
          MPI_Send(source_buf, 2 * Nx_rank, MPI_DOUBLE, 0, 200 + rank,
                   MPI_COMM_WORLD);
        }

        // Wait until all source stuff is sent
        MPI_Barrier(MPI_COMM_WORLD);

        if (rank == 0) { // Rank 0 performs the Poisson Solve

          if (ecouple == 2) {
            simplePoisson(Nx_rank, source, dx, Lx, PoisPot);
          } else {

            if (poissFlavor == 0) { // no E-field
              for (l = 0; l < Nx; l++)
                PoisPot_allranks[l] = 0.0;
            } else if (poissFlavor == 11) // Linear Yukawa
              PoissLinPeriodic1D(Nx, source_allranks, dx, Lx, PoisPot_allranks,
                                 Te_arr_allranks);
            else if (poissFlavor == 12) // Nonlinear Yukawa
              PoissNonlinPeriodic1D(Nx, source_allranks, dx, Lx,
                                    PoisPot_allranks, Te_arr_allranks);
            else if (poissFlavor == 21) // Linear Thomas-Fermi
              PoissLinPeriodic1D_TF(Nx, source_allranks, dx, Lx,
                                    PoisPot_allranks, Te_arr_allranks);
            else if (poissFlavor == 22) // Nonlinear Thomas-Fermi
              PoissNonlinPeriodic1D_TF(Nx, source_allranks, dx, Lx,
                                       PoisPot_allranks, Te_arr_allranks);
            else if (poissFlavor == 42) // ohm law
              fluid_electron_get_eE(Nx, source_allranks, dx, PoisPot_allranks,
                                    Te_arr_allranks);
          }
        }

        // Distribute back to the other ranks
        if (rank == 0) {
          rankOffset = 0;

          // SET LOCAL GHOST CELLS
          if (order == 1) {
            PoisPot[0] = PoisPot_allranks[0];
            if (numRanks == 1)
              PoisPot[Nx_rank + 1] = PoisPot_allranks[0];
            else
              PoisPot[Nx_rank + 1] = PoisPot_allranks[Nx_rank + 1];
          } else if (order == 2) {
            PoisPot[0] = PoisPot_allranks[Nx - 2];
            PoisPot[1] = PoisPot_allranks[Nx - 1];
            if (numRanks == 1) {
              PoisPot[Nx_rank + 2] = PoisPot_allranks[0];
              PoisPot[Nx_rank + 3] = PoisPot_allranks[1];
            } else {
              PoisPot[Nx_rank + 2] = PoisPot_allranks[Nx_rank];
              PoisPot[Nx_rank + 3] = PoisPot_allranks[Nx_rank + 1];
            }
          }

          // Set main body of PoisPot
          for (l = 0; l < Nx_rank; l++) {
            PoisPot[l + order] = PoisPot_allranks[l];
          }

          if (numRanks > 1) {

            rankOffset = Nx_rank;

            for (rankCounter = 1; rankCounter < numRanks - 1; rankCounter++) {

              // Send local ghost cell info
              if (order == 1) {
                source_buf[0] = PoisPot_allranks[rankOffset - 1];
                source_buf[Nx_ranks[rankCounter] + 1] =
                    PoisPot_allranks[rankOffset + Nx_ranks[rankCounter]];
              } else if (order == 2) {
                source_buf[0] = PoisPot_allranks[rankOffset - 2];
                source_buf[1] = PoisPot_allranks[rankOffset - 1];
                source_buf[Nx_ranks[rankCounter] + 2] =
                    PoisPot_allranks[rankOffset + Nx_ranks[rankCounter]];
                source_buf[Nx_ranks[rankCounter] + 3] =
                    PoisPot_allranks[rankOffset + Nx_ranks[rankCounter] + 1];
              }

              // Set main body of PoisPot on rankCounter
              for (l = 0; l < Nx_ranks[rankCounter]; l++) {
                source_buf[l + order] = PoisPot_allranks[rankOffset + l];
              }

              MPI_Send(source_buf, Nx_ranks[rankCounter] + 2 * order,
                       MPI_DOUBLE, rankCounter, rankCounter, MPI_COMM_WORLD);

              rankOffset += Nx_ranks[rankCounter];
            }

            // Deal with periodic BC for rightmost rank

            // SET LOCAL GHOST CELLS
            if (order == 1) {
              source_buf[0] = PoisPot_allranks[rankOffset - 1];
              source_buf[Nx_ranks[numRanks - 1] + 1] = PoisPot_allranks[Nx - 1];
            } else if (order == 2) {
              source_buf[0] = PoisPot_allranks[rankOffset - 2];
              source_buf[1] = PoisPot_allranks[rankOffset - 1];
              source_buf[Nx_ranks[numRanks - 1] + 2] = PoisPot_allranks[0];
              source_buf[Nx_ranks[numRanks - 1] + 3] = PoisPot_allranks[1];
            }

            for (l = 0; l < Nx_ranks[numRanks - 1]; l++) {
              source_buf[l + order] = PoisPot_allranks[rankOffset + l];
            }

            MPI_Send(source_buf, Nx_ranks[numRanks - 1] + 2 * order, MPI_DOUBLE,
                     numRanks - 1, numRanks - 1, MPI_COMM_WORLD);
          }
        } else {
          // Get potential from rank 0

          MPI_Recv(source_buf, Nx_rank + 2 * order, MPI_DOUBLE, 0, rank,
                   MPI_COMM_WORLD, &status);

          // Set Poispot
          for (l = 0; l < Nx_rank + 2 * order; l++) {
            PoisPot[l] = source_buf[l];
          }
        }

        // Finished with determining poisson solve, now advect

        // RK2 Step 2
        for (i = 0; i < nspec; i++) {
          advectTwo_v(f_tmp, f_conv, PoisPot, Z_oned, m[i], i);
        }
        for (l = 0; l < Nx_rank; l++)
          for (i = 0; i < nspec; i++)
            // #pragma omp parallel for private(j)
            for (j = 0; j < Nv * Nv * Nv; j++)
              f[l + order][i][j] =
                  0.5 * (f[l + order][i][j] + f_tmp[l + order][i][j]) +
                  0.25 * f_conv[l + order][i][j];
      }
    }
    t += dt;
    nT++;
  }

  // Store final timstep data

  if (dims == 1) {
    // Calculate moment data in all cells
    for (l = 0; l < Nx_rank; l++) {

      ntot = 0.0;
      rhotot = 0.0;
      for (i = 0; i < nspec; i++) {
        n_oned[l][i] = getDensity(f[l + order][i], i);
        ntot += n_oned[l][i];
        rhotot += m[i] * n_oned[l][i];
        getBulkVel(f[l + order][i], v_oned[l][i], n_oned[l][i], i);
        T_oned[l][i] =
            getTemp(m[i], n_oned[l][i], v_oned[l][i], f[l + order][i], i);
        for (int idim = 0; idim < 3; ++idim) {
          T_dir_oned[l][i][idim] = getDirTemp(m[i], n_oned[l][i], v_oned[l][i],
                                              f[l + order][i], i, idim);
        }
      }
    }

    // moments calculated, now store final step. We are ignoring the e field
    // here

    if (rank == 0) {
      for (l = 0; l < Nx_rank; l++) {
        for (i = 0; i < nspec; i++) {
          fprintf(outputFileDens[i], "%e ", n_oned[l][i]);
          fprintf(outputFileVelo[i], "%e ", v_oned[l][i][0]);
          fprintf(outputFileTemp[i], "%e ", T_oned[l][i]);
          fprintf(outputFileDirTemp[i], "%e %e %e ", T_dir_oned[l][i][0],
                  T_dir_oned[l][i][1], T_dir_oned[l][i][2]);
        }
      }

      // get from other ranks
      for (rankCounter = 1; rankCounter < numRanks; rankCounter++) {
        for (s = 0; s < nspec; s++) {
          MPI_Recv(momentBuffer, 6 * Nx_ranks[rankCounter], MPI_DOUBLE,
                   rankCounter, 100 + s, MPI_COMM_WORLD, &status);
          for (l = 0; l < Nx_ranks[rankCounter]; l++) {
            fprintf(outputFileDens[s], "%e ", momentBuffer[0 + 6 * l]);
            fprintf(outputFileVelo[s], "%e ", momentBuffer[1 + 6 * l]);
            fprintf(outputFileTemp[s], "%e ", momentBuffer[2 + 6 * l]);
            fprintf(outputFileDirTemp[s], "%e %e %e ", momentBuffer[3 + 6 * l],
                    momentBuffer[4 + 6 * l], momentBuffer[5 + 6 * l]);
          }
        }
      }

      // Close out this timestep
      fprintf(outputFileTime, "%e\n", t);
      for (i = 0; i < nspec; i++) {
        fprintf(outputFileDens[i], "\n");
        fprintf(outputFileVelo[i], "\n");
        fprintf(outputFileTemp[i], "\n");
        fprintf(outputFileDirTemp[i], "\n");
      }

      if (outputDist > 0) {
        store_distributions_inhomog(f, input_filename, nT, t, rank);
      }
      if (outputDist == 3) {
        char BGK_name_buffer[256];
        for (int s = 0; s < nspec; ++s) {
          sprintf(BGK_name_buffer, "Data/%s_spec%d.step%d_rank%d_BGK.dat",
                  input_filename, s, nT, rank);
          FILE *fid_store = fopen(BGK_name_buffer, "w");
          for (int xj = 0; xj < Nx_rank; ++xj) {
            fwrite(f_tmp[xj + order][s], sizeof(double), Nv * Nv * Nv,
                   fid_store);
          }
          fclose(fid_store);
        }
      }
    } else { // send to rank 0 for output purposes
      for (s = 0; s < nspec; s++) {
        for (l = 0; l < Nx_rank; l++) {
          momentBuffer[0 + 6 * l] = n_oned[l][s];
          momentBuffer[1 + 6 * l] = v_oned[l][s][0];
          momentBuffer[2 + 6 * l] = T_oned[l][s];
          momentBuffer[3 + 6 * l] = T_dir_oned[l][s][0];
          momentBuffer[4 + 6 * l] = T_dir_oned[l][s][1];
          momentBuffer[5 + 6 * l] = T_dir_oned[l][s][2];
        }
        MPI_Send(momentBuffer, 6 * Nx_rank, MPI_DOUBLE, 0, 100 + s,
                 MPI_COMM_WORLD);
      }
      if (outputDist > 0)
        store_distributions_inhomog(f, input_filename, nT, t, rank);
      if (outputDist == 3) {
        char BGK_name_buffer[256];
        for (s = 0; s < nspec; ++s) {
          sprintf(BGK_name_buffer, "Data/%s_spec%d.step%d_rank%d_BGK.dat",
                  input_filename, s, nT, rank);
          FILE *fid_store = fopen(BGK_name_buffer, "w");
          for (int xj = 0; xj < Nx_rank; ++xj) {
            fwrite(f_tmp[xj + order][s], sizeof(double), Nv * Nv * Nv,
                   fid_store);
          }
        }
      }
    }

    MPI_Barrier(MPI_COMM_WORLD);
  }

  if (restartFlag > 0) {
    if (dims == 0)
      store_distributions_homog(f_zerod, t, -1 * nT, input_filename);
    else
      store_distributions_inhomog(f, input_filename, nT, t, rank);
  }

  // clean up

  // universal to both 0D and 1D
  for (i = 0; i < nspec; i++) {
    free(c[i]);
    free(wts[i]);
  }
  free(c);
  free(wts);
  free(vref);
  free(Lv);

  if (rank == 0) {
    for (i = 0; i < nspec; i++) {
      fclose(outputFileDens[i]);
      fclose(outputFileVelo[i]);
      fclose(outputFileTemp[i]);
      fclose(outputFileDirTemp[i]);
    }
    free(outputFileDens);
    free(outputFileVelo);
    free(outputFileTemp);
    free(outputFileDirTemp);
  }

  if (dims == 0) {
    for (i = 0; i < nspec; i++) {
      free(f_zerod[i]);
    }
    free(f_zerod);
  }

  if (dims == 1) {
    free(x);
    free(dxarray);

    if (rank == 0) {
      fclose(outputFilePoiss);
      fclose(outputFilePhi);
    }

    for (l = 0; l < Nx_rank; l++) {
      free(n_oned[l]);
      free(T_oned[l]);
      free(v0_oned[l]);
      for (i = 0; i < nspec; i++) {
        free(v_oned[l][i]);
        free(PTens_oned[l][i]);
        free(f[l][i]);
        free(f_tmp[l][i]);
        free(f_conv[l][i]);
      }
      free(v_oned[l]);
      free(PTens_oned[l]);

      free(f[l]);
      free(f_tmp[l]);
      free(f_conv[l]);
    }

    free(f);
    free(f_tmp);
    free(f_conv);
    free(n_oned);
    free(v_oned);
    free(PTens_oned);
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

  MPI_Finalize();

  return 0;
}
