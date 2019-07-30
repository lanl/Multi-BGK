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
#include "parallel_poisson.h"
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
  int rankOffset;
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
  double *T0_oned;
  double **T_for_zbar;
  double *T_max, T0_max;

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
  int bcs; //boundary conditions. 1 if fixed, 0 if periodic.
  //parameters needed for fixed bcs poisson solver
  //This could actually be improved by only tracking the total charge density....
  double *T0_bcs; //index for left or right
  double **n_bcs; //n_bcs[left/right][spec]
  double **Z_bcs; //Z_bcs[left/right][spec]
  double *Te_bcs; //Te_bcs[left/right]

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
  double eq_rtol; 

  int hydro_kinscheme_flag, zerot_flag, eq_flag;

  /**********************************
         I/O Setup
  **********************************/

  char input_filename[100];
  char input_file_data_filename[100];
  strcpy(input_filename, argv[1]);

  read_input(&nspec, &dims, &Nx, &Lx, &bcs, &Nv, &v_sigma, &discret, &poissFlavor, &m,
             &Z_max, &order, &im_ex, &dt, &tfinal, &numint, &intervalLimits,
             &ndens_int, &velo_int, &T_int, &ecouple, &ionFix, &Te_start,
             &Te_ref, &CL_type, &ion_type, &MT_or_TR, &n_zerod, &v_val,
             &T_zerod, &dataFreq, &outputDist, &RHS_tol, &BGK_type, &beta,
             &hydro_kinscheme_flag, &eq_flag, &eq_rtol, &input_file_data_flag,
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

  if(argc < 5){
    zerot_flag = 0;
  } else {
    zerot_flag = atoi(argv[4]);
  }

  if (argc > 2){
    printf("restart %d tau %d zerot %d\n", restartFlag, tauFlag, zerot_flag);
  }

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
    //store physical mesh
    //strcpy(x_path, output_path);
    //strcat(x_path, "_x");
 
    if ((restartFlag == 2) || (restartFlag == 4)){
      outputFilePoiss = fopen(poiss_path, "a");
      //outputFile_x = fopen(x_path, "a");
    } else {
      outputFilePoiss = fopen(poiss_path, "w");
      //outputFile_x = fopen(x_path, 'w');
    }
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

    // Physical grid allocation and initialization
    make_mesh(Nx, Lx, order, &Nx_rank, &Nx_ranks, &x, &dxarray);
    momentBuffer = malloc(3 * (Nx_rank + 1) * sizeof(double));

    dx = Lx / Nx;

    // Write the full mesh for storage purposes
    if (rank == 0) {
      for (l = 0; l < Nx; l++){
        fprintf(outputFile_x, "%+le ", (l + 0.5) * dx - 0.5 * Lx);
      }
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
    source = malloc(Nx_rank * sizeof(double));
    source_buf = malloc(2 * (Nx_rank + 1) * sizeof(double));
    PoisPot = malloc((Nx_rank + 2 * order) * sizeof(double));
    Te_arr = malloc(Nx_rank * sizeof(double));
    Z_oned = malloc(Nx_rank * sizeof(double*));

    for(l = 0; l < Nx_rank; l++){
      Z_oned[l] = malloc(nspec * sizeof(double*));
      for(i=0; i < nspec; i++){
        //initialize at full ionization
        Z_oned[l][i] = Z_max[i];
      }
    }


    if(rank==0){
      if(bcs == 0){
        source_allranks = malloc(Nx * sizeof(double));
        Te_arr_allranks = malloc(Nx * sizeof(double));
        PoisPot_allranks = malloc(Nx * sizeof(double));
      } else if(bcs == 1){
        //all ghost cells on the left & right have the same source/potential value,
        //i.e. you are now doing kinetic theory in between two infinite square plates,
        //connected to two proportionally unlimited voltage sources. FUN!
        source_allranks = malloc((Nx + 2) * sizeof(double));
        Te_arr_allranks = malloc((Nx + 2) * sizeof(double));
        PoisPot_allranks = malloc((Nx + 2) * sizeof(double));
        T0_bcs = malloc(2 * sizeof(double));
        Te_bcs = malloc(2 * sizeof(double));
        n_bcs = malloc(2 * sizeof(double*));
        Z_bcs = malloc(2*sizeof(double*));
        for(i=0; i<2; i++){
          n_bcs[i] = malloc(nspec*sizeof(double));
          Z_bcs[i] = malloc(nspec*sizeof(double));
        }
      }
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

    initialize_transport(bcs, Nv, Nx_rank, nspec, x, dxarray, Lx, c, order, dt);

    // check if we are loading moment data from a file
    // Note: This is not yet MPI-ified
    if (input_file_data_flag) {
      initialize_sol_inhom_file(f, Nx_rank, nspec, Nv, order, c, m, n_oned,
                                v_oned, T_oned);
    } else {
      initialize_sol_inhom(&rank, &numRanks,f, numint, intervalLimits, ndens_int, velo_int,
                           T_int, Nx_rank, x, nspec, Nv, order, c, m, n_oned,
                           v_oned, T_oned);
      if(rank==0){
        //set left and right params. 
        if(bcs==1){
          double tmp_mixv = 0;
          double tmp_totndens, tmp_totdens;
          //set left bcs (the first interval)
          for(i=0; i<nspec; i++){
            n_bcs[0][i] = ndens_int[i];
            Z_bcs[0][i] = Z_max[i];
            tmp_totndens += ndens_int[i];
            tmp_totdens += m[i]*ndens_int[i];
            tmp_mixv += m[i]*ndens_int[i]*velo_int[i];
          }
          tmp_mixv = tmp_mixv/tmp_totdens;
          for(i=0; i<nspec; i++){
            T0_bcs[0] += ndens_int[i]*(T_int[i] + m[i]*(velo_int[i]-tmp_mixv)*(velo_int[i]-tmp_mixv)/3.0);
          }
          T0_bcs[0] = T0_bcs[0]/tmp_totndens;
          //if only one interval, set the same bounds on each end.
          if(numint == 1){
            T0_bcs[1] = T0_bcs[0];
            for(i=0; i<nspec; i++){
              n_bcs[1][i] = n_bcs[0][i];
              Z_bcs[1][i] = n_bcs[0][i];
            }
          }else{
            tmp_totndens = 0.0;
            tmp_totdens = 0.0;
            tmp_mixv = 0.0;
            l = (numint-1)*nspec;
            for(i = 0; i<nspec; i++){
              n_bcs[1][i] = ndens_int[l + i];
              Z_bcs[1][i] = Z_max[i];
              tmp_totndens += ndens_int[l+i];
              tmp_totdens += m[i]*ndens_int[l+i];
              tmp_mixv += m[i]*ndens_int[l+i]*velo_int[l+i];
            }
            tmp_mixv = tmp_mixv/tmp_totdens;
            for(i=0; i<nspec; i++){
              T0_bcs[1] += n_bcs[1][i]*(T_int[l+i] + m[i]*(velo_int[l+i]-tmp_mixv)*(velo_int[l+i]-tmp_mixv)/3.0);
            }
            T0_bcs[1] = T0_bcs[1]/tmp_totndens;
          }
          printf("T0_left:%g, T0_right:%g\n", T0_bcs[0], T0_bcs[1]);
        }
      }
    }
    if (rank == 0) {
      printf("Initial condition setup complete\n");
      fflush(stdout);
    }
  }

  initialize_moments(Nv, nspec, c, wts);
  initialize_BGK(nspec, Nv, m, c, order, ecouple, CL_type, ion_type, MT_or_TR,
                 tauFlag, input_filename);

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
  if(eq_flag == 0){
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

        //output if time to or initial data.
        if (outcount == dataFreq || nT==0) {
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

       /*
       Regarding the number of input arguments:
       "It is true, we shall be monsters, cut off from the world; 
       but on that account we shall be more attached to one another" - MS, F
       */
       poisson_solver(MPI_COMM_WORLD, &rank, &numRanks, &status, &Nx, &Nx_ranks, &Nx_rank,
                     &ecouple, &bcs, &poissFlavor, &ionFix, &nspec, &Z_max, &Te_ref, 
                     &Te_start, &order, &dx, &Lx, &t, &tfinal, &Te_arr, &Te_arr_allranks,
                     &T0_oned, &n_oned, &Z_oned, &source, &source_buf, &source_allranks,
                     &PoisPot, &PoisPot_allranks, &T0_bcs, &n_bcs, &Z_bcs, &Te_bcs);

        // Moments and initial electric field calculated - save if needed

        // MPI-ified output.
        if (outcount == dataFreq) {
          if (rank == 0) {
            for (l = 0; l < Nx_rank; l++) {
              for (i = 0; i < nspec; i++) {
                T_oned[l][i] =
                    getTemp(m[i], n_oned[l][i], v_oned[l][i], f[l + order][i], i);
                fprintf(outputFileDens[i], "%e ", n_oned[l][i]);
                fprintf(outputFileVelo[i], "%e ", v_oned[l][i][0]);
                fprintf(outputFileTemp[i], "%e ", T_oned[l][i]);
              }
            }

            // get from other ranks
            for (rankCounter = 1; rankCounter < numRanks; rankCounter++) {
              for (s = 0; s < nspec; s++) {
                MPI_Recv(momentBuffer, 3 * Nx_ranks[rankCounter], MPI_DOUBLE,
                         rankCounter, 100 + s, MPI_COMM_WORLD, &status);
                for (l = 0; l < Nx_ranks[rankCounter]; l++) {
                  fprintf(outputFileDens[s], "%e ", momentBuffer[0 + 3 * l]);
                  fprintf(outputFileVelo[s], "%e ", momentBuffer[1 + 3 * l]);
                  fprintf(outputFileTemp[s], "%e ", momentBuffer[2 + 3 * l]);
                }
              }
            }

            // Print poisson information - TODO double check units
            if(bcs == 0){
              fprintf(outputFilePoiss, "%e ",
                      0.5*(PoisPot_allranks[1] - PoisPot_allranks[Nx - 1]) / dx);
              for (l = 1; l < Nx - 1; l++) {
                fprintf(outputFilePoiss, "%e ",
                        0.5*(PoisPot_allranks[l + 1] - PoisPot_allranks[l - 1]) / dx);
              }
              fprintf(outputFilePoiss, "%e ",
                      0.5*(PoisPot_allranks[0] - PoisPot_allranks[Nx - 2]) / dx);
              fprintf(outputFilePoiss, "\n");
            }else{
              for(l = 1; l < Nx+1; l++){
                fprintf(outputFilePoiss, "%e ",
                        0.5*(PoisPot_allranks[l + 1] - PoisPot_allranks[l - 1]) / dx);
              }
            }

            // Close out this timestep
            fprintf(outputFileTime, "%e\n", t);
            for (i = 0; i < nspec; i++) {
              fprintf(outputFileDens[i], "\n");
              fprintf(outputFileVelo[i], "\n");
              fprintf(outputFileTemp[i], "\n");
            }
            outcount = 0;
          } else { // send to rank 0 for output purposes
            for (s = 0; s < nspec; s++) {
              for (l = 0; l < Nx_rank; l++) {
                T_oned[l][s] =
                    getTemp(m[s], n_oned[l][s], v_oned[l][s], f[l + order][s], s);
                momentBuffer[0 + 3 * l] = n_oned[l][s];
                momentBuffer[1 + 3 * l] = v_oned[l][s][0];
                momentBuffer[2 + 3 * l] = T_oned[l][s];
              }
              MPI_Send(momentBuffer, 3 * Nx_rank, MPI_DOUBLE, 0, 100 + s,
                       MPI_COMM_WORLD);
            }
            outcount = 0;
          }

          MPI_Barrier(MPI_COMM_WORLD);
          if(outputDist==1){
            store_distributions_inhomog(&numRanks, &rank, &order, f, input_filename, nT);
          }
        }
        outcount += 1;

        // IO done, advance to the actual solution...

        MPI_Barrier(MPI_COMM_WORLD);

        if (order == 1) {
          // ADVECT

          for (i = 0; i < nspec; i++) {
            advectOne(f, PoisPot, Z_oned, m[i], i);
          }

          if (!(BGK_type < 0)) {
            // COLLIDE
            if (im_ex == 0) {
              for (l = 0; l < Nx_rank; l++) {
                BGK_ex(f[l + order], f_tmp[l + order], Z_oned[l], dt, Te_arr[l]);
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

          for (l = 0; l < Nx_rank; l++){
            for (i = 0; i < nspec; i++){
              //#pragma omp parallel for private(j)
              for (j = 0; j < Nv * Nv * Nv; j++){
                f_tmp[l + order][i][j] =
                    f[l + order][i][j] + 0.5 * f_conv[l + order][i][j];
              }
              n_oned[l][i] = getDensity(f_tmp[l+order][i], i);
            }
          }
          // Do second step of Poisson solve 

          //CONCERN: WE DO A SECOND POISSON SOLVE WITHOUT FIRST UPDATING THE MIXTURE TEMPERATURE
          //CHANGE RESULTING FROM THE ABOVE ADVECTION.

          //see line 1044
          poisson_solver(MPI_COMM_WORLD, &rank, &numRanks, &status, &Nx, &Nx_ranks, &Nx_rank,
                         &ecouple, &bcs, &poissFlavor, &ionFix, &nspec, &Z_max, &Te_ref, 
                         &Te_start, &order, &dx, &Lx, &t, &tfinal, &Te_arr, &Te_arr_allranks,
                         &T0_oned, &n_oned, &Z_oned, &source, &source_buf, &source_allranks,
                         &PoisPot, &PoisPot_allranks, &T0_bcs, &n_bcs, &Z_bcs, &Te_bcs);
          // Finished with determining poisson solve, now advect

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

          if (!(BGK_type == -1)) {
            // Next Strang step - RK2 for collision with timstep dt
            for (l = 0; l < Nx_rank; l++) {

              if (im_ex == 0) {
                // Step 1
                BGK_ex(f[l + order], f_conv[l + order], Z_oned[l], dt, Te_arr[l]);
                for (i = 0; i < nspec; i++)
                  for (j = 0; j < Nv * Nv * Nv; j++)
                    f_conv[l+order][i][j] = dt*f_conv[l+order][i][j];
                    f_tmp[l + order][i][j] += (f_conv[l+order][i][j] < -1*f_tmp[l+order][i][j]) ? 0.0 : f_conv[l+order][i][j];

                // Step 2
                BGK_ex(f_tmp[l + order], f_conv[l + order], Z_oned[l], dt,
                       Te_arr[l]);
                for (i = 0; i < nspec; i++)
                  for (j = 0; j < Nv * Nv * Nv; j++)
                    f_conv[l+order][i][j] = 0.5*dt*f_conv[l+order][i][j];
                    f[l+order][i][j] = 0.5*(f[l+order][i][j] + f_tmp[l+order][i][j]);
                    f_tmp[l + order][i][j] += (f_conv[l+order][i][j] < -1*f_tmp[l+order][i][j]) ? 0.0 : f_conv[l+order][i][j];
                    f[l+order][i][j] += (f_conv[l+order][i][j] < -1*f[l+order][i][j]) ? 0.0 : f_conv[l+order][i][j];
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
              //#pragma omp parallel for private(j)
              for (j = 0; j < Nv * Nv * Nv; j++)
                f_tmp[l + order][i][j] =
                    f[l + order][i][j] + 0.5 * f_conv[l + order][i][j];

          // RK2 Step 2
          for (i = 0; i < nspec; i++) {
            advectTwo_x(f_tmp, f_conv, i);
          }
          for (l = 0; l < Nx_rank; l++){
            for (i = 0; i < nspec; i++){
              //#pragma omp parallel for private(j)
              for (j = 0; j < Nv * Nv * Nv; j++){
                f[l + order][i][j] =
                    0.5 * (f[l + order][i][j] + f_tmp[l + order][i][j]) +
                    0.25 * f_conv[l + order][i][j];
              }
              n_oned[l][i] = getDensity(f[l + order][i], i);
            }
          }
          // Last strang step - V advection with timestep dt/2 (combine?)

          // Recalc Poiss

          //CONCERN: WE DO A SECOND POISSON SOLVE WITHOUT FIRST UPDATING THE MIXTURE TEMPERATURE
          //CHANGE RESULTING FROM THE ABOVE ADVECTION.

          //see line 1044
          poisson_solver(MPI_COMM_WORLD, &rank, &numRanks, &status, &Nx, &Nx_ranks, &Nx_rank,
                         &ecouple, &bcs, &poissFlavor, &ionFix, &nspec, &Z_max, &Te_ref, 
                         &Te_start, &order, &dx, &Lx, &t, &tfinal, &Te_arr, &Te_arr_allranks,
                         &T0_oned, &n_oned, &Z_oned, &source, &source_buf, &source_allranks,
                         &PoisPot, &PoisPot_allranks, &T0_bcs, &n_bcs, &Z_bcs, &Te_bcs);

          // Finished with determining poisson solve, now advect

          // RK2 Step 1
          for (i = 0; i < nspec; i++) {
            advectTwo_v(f, f_conv, PoisPot, Z_oned, m[i], i);
          }

          for (l = 0; l < Nx_rank; l++){
              for (i = 0; i < nspec; i++){
                //#pragma omp parallel for private(j)
                for (j = 0; j < Nv * Nv * Nv; j++){
                  f_tmp[l + order][i][j] =
                      f[l + order][i][j] + 0.5 * f_conv[l + order][i][j];
                }
                n_oned[l][i] = getDensity(f[l + order][i], i);
              }
            } 

            // Recalc Poiss

            //CONCERN: WE DO A SECOND POISSON SOLVE WITHOUT FIRST UPDATING THE MIXTURE TEMPERATURE
            //CHANGE RESULTING FROM THE ABOVE ADVECTION.

            //see line 1044
            poisson_solver(MPI_COMM_WORLD, &rank, &numRanks, &status, &Nx, &Nx_ranks, &Nx_rank,
                           &ecouple, &bcs, &poissFlavor, &ionFix, &nspec, &Z_max, &Te_ref, 
                           &Te_start, &order, &dx, &Lx, &t, &tfinal, &Te_arr, &Te_arr_allranks,
                           &T0_oned, &n_oned, &Z_oned, &source, &source_buf, &source_allranks,
                           &PoisPot, &PoisPot_allranks, &T0_bcs, &n_bcs, &Z_bcs, &Te_bcs); 

            // Finished with determining poisson solve, now advect

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
  } else {
    //ensure we step at least once
    double h = 1.0;
    double dh = eq_rtol*h + 1.0;
    while (dh > eq_rtol*h) {

      if (rank == 0) {
        printf("At time %g\n", t);
      }

      if (dims == 0) {
        //if dims == 0, the stagnation of the entropy is used rather than the grid norm
        //to determine equilibration.
        // GET MOMENTS

        ntot = 0.0;
        rhotot = 0.0;
        Htot_prev = Htot;
        Htot = 0.0;
        h = fabs(Htot_prev);

        for (i = 0; i < nspec; i++) {
          n_zerod[i] = getDensity(f_zerod[i], i);
          ntot += n_zerod[i];
          rhotot += m[i] * n_zerod[i];

          H_spec_prev[i] = H_spec[i];
          H_spec[i] = getH(n_zerod[i], f_zerod[i], i);
          Htot += H_spec[i];

          getBulkVel(f_zerod[i], v_zerod[i], n_zerod[i], i);
        }
        dh = abs(Htot - Htot_prev);

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


        BGK_norm(f_zerod, BGK_f_minus_eq, Z_zerod, dt, T0);

        // check to make sure that we don't need to stop based on collision rate.
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

         for (i = 0; i < nspec; i++){
            for (j = 0; j < Nv * Nv * Nv; j++){
              f_zerod[i][j] += dt * f_zerod_tmp[i][j];
            }
          }
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
        
        //put value of previous global H in h,
        h = Htot;
        Htot = 0.0;
        // Calculate moment data in all cells
        for (l = 0; l < Nx_rank; l++) {
          ntot = 0.0;
          rhotot = 0.0;
          for (i = 0; i < nspec; i++) {
            n_oned[l][i] = getDensity(f[l + order][i], i);
            ntot += n_oned[l][i];
            rhotot += m[i] * n_oned[l][i];
            Htot += getH(n_oned[l][i], f[l+order][i], i);
            getBulkVel(f[l + order][i], v_oned[l][i], n_oned[l][i], i);
            T_oned[l][i] = getTemp(m[i], n_oned[l][i], v_oned[l][i], f[l+order][i], i);
          }

          // get mixture mass avg velocity
          for (j = 0; j < 3; j++) {
            v0_oned[l][j] = 0.0;
            for (i = 0; i < nspec; i++)
              v0_oned[l][j] += m[i] * n_oned[l][i] * v_oned[l][i][j];
            v0_oned[l][j] = v0_oned[l][j] / rhotot;
          }
          //set T0 
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
        //put the current global entropy from every rank into dh.
        MPI_Allreduce(&Htot, &dh, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        dh = fabs(dh - h);
        h = fabs(h);
        
        // Flag - do we want to run this like the kinetic scheme for hydro
        if (hydro_kinscheme_flag == 1) {
          for (l = 0; l < Nx_rank; l++)
            for (i = 0; i < nspec; i++)
              GetMaxwell(m[i], n_oned[l][i], v_oned[l][i], T_oned[l][i],
                         f[l + order][i], i);
        }

       /*
       Regarding the number of input arguments:
       "It is true, we shall be monsters, cut off from the world; 
       but on that account we shall be more attached to one another" - MS, F
       */
       poisson_solver(MPI_COMM_WORLD, &rank, &numRanks, &status, &Nx, &Nx_ranks, &Nx_rank,
                     &ecouple, &bcs, &poissFlavor, &ionFix, &nspec, &Z_max, &Te_ref, 
                     &Te_start, &order, &dx, &Lx, &t, &tfinal, &Te_arr, &Te_arr_allranks,
                     &T0_oned, &n_oned, &Z_oned, &source, &source_buf, &source_allranks,
                     &PoisPot, &PoisPot_allranks, &T0_bcs, &n_bcs, &Z_bcs, &Te_bcs);

        // Moments and initial electric field calculated - don't save since equilibrating.
       // IO done, advance to the actual solution...

        MPI_Barrier(MPI_COMM_WORLD);

        if (order == 1) {
          // ADVECT

          for (i = 0; i < nspec; i++) {
            advectOne(f, PoisPot, Z_oned, m[i], i);
          }

          if (!(BGK_type < 0)) {
            // COLLIDE
            if (im_ex == 0) {
              for (l = 0; l < Nx_rank; l++) {
                BGK_ex(f[l + order], f_tmp[l + order], Z_oned[l], dt, Te_arr[l]);
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
  
          for (l = 0; l < Nx_rank; l++){
            for (i = 0; i < nspec; i++){
              //#pragma omp parallel for private(j)
              for (j = 0; j < Nv * Nv * Nv; j++){
                f_tmp[l + order][i][j] =
                    f[l + order][i][j] + 0.5 * f_conv[l + order][i][j];
              }
              n_oned[l][i] = getDensity(f_tmp[l+order][i], i);
            }
          }
          // Do second step of Poisson solve 

          //CONCERN: WE DO A SECOND POISSON SOLVE WITHOUT FIRST UPDATING THE MIXTURE TEMPERATURE
          //CHANGE RESULTING FROM THE ABOVE ADVECTION.

          //see line 1044
          poisson_solver(MPI_COMM_WORLD, &rank, &numRanks, &status, &Nx, &Nx_ranks, &Nx_rank,
                         &ecouple, &bcs, &poissFlavor, &ionFix, &nspec, &Z_max, &Te_ref, 
                         &Te_start, &order, &dx, &Lx, &t, &tfinal, &Te_arr, &Te_arr_allranks,
                         &T0_oned, &n_oned, &Z_oned, &source, &source_buf, &source_allranks,
                         &PoisPot, &PoisPot_allranks, &T0_bcs, &n_bcs, &Z_bcs, &Te_bcs);
          // Finished with determining poisson solve, now advect

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

          if (!(BGK_type == -1)) {
            // Next Strang step - RK2 for collision with timstep dt
            for (l = 0; l < Nx_rank; l++) {

              if (im_ex == 0) {
                // Step 1
                BGK_ex(f[l + order], f_conv[l + order], Z_oned[l], dt, Te_arr[l]);
                for (i = 0; i < nspec; i++){
                  for (j = 0; j < Nv * Nv * Nv; j++){
                    f_conv[l+order][i][j] = dt*f_conv[l+order][i][j];
                    if(f_conv[l+order][i][j] >= -1.0*(f[l+order][i][j])){
		      f_tmp[l+order][i][j] = f[l+order][i][j] + f_conv[l+order][i][j];
		    }else{
                      f_tmp[l+order][i][j] = f[l+order][i][j];
		    }
		  }
                }

                // Step 2
                BGK_ex(f_tmp[l + order], f_conv[l + order], Z_oned[l], dt,
                       Te_arr[l]);
                for (i = 0; i < nspec; i++){
                  for (j = 0; j < Nv * Nv * Nv; j++){
                    f_conv[l+order][i][j] = 0.5*dt*f_conv[l+order][i][j];
                    f[l+order][i][j] = 0.5*(f[l+order][i][j] + f_tmp[l+order][i][j]);
                    if(f_conv[l+order][i][j] > -1.0*f[l+order][i][j]){
		      f[l+order][i][j] = f[l+order][i][j] + f_conv[l+order][i][j];
		    }
		  }
                }
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
              //#pragma omp parallel for private(j)
              for (j = 0; j < Nv * Nv * Nv; j++)
                f_tmp[l + order][i][j] =
                    f[l + order][i][j] + 0.5 * f_conv[l + order][i][j];

          // RK2 Step 2
          for (i = 0; i < nspec; i++) {
            advectTwo_x(f_tmp, f_conv, i);
          }
          for (l = 0; l < Nx_rank; l++){
            for (i = 0; i < nspec; i++){
              //#pragma omp parallel for private(j)
              for (j = 0; j < Nv * Nv * Nv; j++){
                f[l + order][i][j] =
                    0.5 * (f[l + order][i][j] + f_tmp[l + order][i][j]) +
                    0.25 * f_conv[l + order][i][j];
              }
              n_oned[l][i] = getDensity(f[l + order][i], i);
            }
          }
          // Last strang step - V advection with timestep dt/2 (combine?)

          // Recalc Poiss

          //CONCERN: WE DO A SECOND POISSON SOLVE WITHOUT FIRST UPDATING THE MIXTURE TEMPERATURE
          //CHANGE RESULTING FROM THE ABOVE ADVECTION.

          //see line 1044
          poisson_solver(MPI_COMM_WORLD, &rank, &numRanks, &status, &Nx, &Nx_ranks, &Nx_rank,
                         &ecouple, &bcs, &poissFlavor, &ionFix, &nspec, &Z_max, &Te_ref, 
                         &Te_start, &order, &dx, &Lx, &t, &tfinal, &Te_arr, &Te_arr_allranks,
                         &T0_oned, &n_oned, &Z_oned, &source, &source_buf, &source_allranks,
                         &PoisPot, &PoisPot_allranks, &T0_bcs, &n_bcs, &Z_bcs, &Te_bcs);

          // Finished with determining poisson solve, now advect

          // RK2 Step 1
          for (i = 0; i < nspec; i++) {
            advectTwo_v(f, f_conv, PoisPot, Z_oned, m[i], i);
          }

          for (l = 0; l < Nx_rank; l++){
              for (i = 0; i < nspec; i++){
                //#pragma omp parallel for private(j)
                for (j = 0; j < Nv * Nv * Nv; j++){
                  f_tmp[l + order][i][j] =
                      f[l + order][i][j] + 0.5 * f_conv[l + order][i][j];
                }
                n_oned[l][i] = getDensity(f[l + order][i], i);
              }
            } 

            // Recalc Poiss

            //CONCERN: WE DO A SECOND POISSON SOLVE WITHOUT FIRST UPDATING THE MIXTURE TEMPERATURE
            //CHANGE RESULTING FROM THE ABOVE ADVECTION.

            //see line 1044
            poisson_solver(MPI_COMM_WORLD, &rank, &numRanks, &status, &Nx, &Nx_ranks, &Nx_rank,
                           &ecouple, &bcs, &poissFlavor, &ionFix, &nspec, &Z_max, &Te_ref, 
                           &Te_start, &order, &dx, &Lx, &t, &tfinal, &Te_arr, &Te_arr_allranks,
                           &T0_oned, &n_oned, &Z_oned, &source, &source_buf, &source_allranks,
                           &PoisPot, &PoisPot_allranks, &T0_bcs, &n_bcs, &Z_bcs, &Te_bcs); 

            // Finished with determining poisson solve, now advect

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
    }
  }
  // Store final timstep data
  
  if(dims == 1) {
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

      }
    }
   
    //moments calculated, now store final step. We aren't ignoring the e field here

    if(rank == 0) {
      for (l = 0; l < Nx_rank; l++) {
        for (i = 0; i < nspec; i++) {
          fprintf(outputFileDens[i], "%e ", n_oned[l][i]);
          fprintf(outputFileVelo[i], "%e ", v_oned[l][i][0]);
          fprintf(outputFileTemp[i], "%e ", T_oned[l][i]);
        }
      }
      
      // get from other ranks
      for (rankCounter = 1; rankCounter < numRanks; rankCounter++) {
        for (s = 0; s < nspec; s++) {
          MPI_Recv(momentBuffer, 3 * Nx_ranks[rankCounter], MPI_DOUBLE,
                   rankCounter, 100 + s, MPI_COMM_WORLD, &status);
          for (l = 0; l < Nx_ranks[rankCounter]; l++) {
            fprintf(outputFileDens[s], "%e ", momentBuffer[0 + 3 * l]);
            fprintf(outputFileVelo[s], "%e ", momentBuffer[1 + 3 * l]);
            fprintf(outputFileTemp[s], "%e ", momentBuffer[2 + 3 * l]);
          }
        }
      }

      //output e field.
      if(bcs == 0){
        fprintf(outputFilePoiss, "%e ",
              0.5*(PoisPot_allranks[1] - PoisPot_allranks[Nx - 1]) / dx);
        for (l = 1; l < Nx - 1; l++) {
          fprintf(outputFilePoiss, "%e ",
                  0.5*(PoisPot_allranks[l + 1] - PoisPot_allranks[l - 1]) / dx);
        }
        fprintf(outputFilePoiss, "%e ",
                  0.5*(PoisPot_allranks[0] - PoisPot_allranks[Nx - 2]) / dx);
        fprintf(outputFilePoiss, "\n");
      }else{
        for(l = 1; l < Nx+1; l++){
          fprintf(outputFilePoiss, "%e ",
                  0.5*(PoisPot_allranks[l + 1] - PoisPot_allranks[l - 1]) / dx);
        }
      }


      
      // Close out this timestep
      fprintf(outputFileTime, "%e\n", t);
      for (i = 0; i < nspec; i++) {
        fprintf(outputFileDens[i], "\n");
        fprintf(outputFileVelo[i], "\n");
        fprintf(outputFileTemp[i], "\n");
      }
    } else { // send to rank 0 for output purposes
      for (s = 0; s < nspec; s++) {
        for (l = 0; l < Nx_rank; l++) {
          momentBuffer[0 + 3 * l] = n_oned[l][s];
          momentBuffer[1 + 3 * l] = v_oned[l][s][0];
          momentBuffer[2 + 3 * l] = T_oned[l][s];
        }
        MPI_Send(momentBuffer, 3 * Nx_rank, MPI_DOUBLE, 0, 100 + s,
                 MPI_COMM_WORLD);
      }
      
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(outputDist==1){
      store_distributions_inhomog(&numRanks, &rank, &order, f, input_filename, nT);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
   


  if ((dims == 0) && (restartFlag > 0))
    store_distributions_homog(f_zerod, t, -1 * nT, input_filename);

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
    }
    free(outputFileDens);
    free(outputFileVelo);
    free(outputFileTemp);
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
    }

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

  MPI_Finalize();

  return 0;
}
