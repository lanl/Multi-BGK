#include "input.h"
#include "units/unit_data.c"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*!
 *  This function reads the input file
 */

void read_input(int *nspec, int *dims, int *Nx, double *Lx, int *bcs,  int *Nv,
                double *v_sigma, int *discret, int *poissFlavor, double **m,
                double **Z, int *order, int *im_ex, double *dt, double *tfinal,
                int *numint, double **intervalLimits, double **ndens_int,
                double **velo_int, double **T_int, int *ecouple, int *ionFix,
                double *Te_start, double *Te_end, int *CL_type, int *ion_type,
                int *MT_or_TR, double **n, double **u, double **T,
                int *dataFreq, int *outputDist, double *RHS_tol, int *BGK_type,
                double *beta, int *hydro_flag, int *eqflag, double *eqrtol, int *input_file_data_flag,
                char *input_file_data_filename, char *inputFilename) {

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  char line[10000] = {"dummy"};
  char input_path[100] = {"./input/"};
  FILE *input_file;
  int i;
  int int_id = 0;

  /*Set input parameters to default values*/

  set_default_values(Nx, Lx, bcs, Nv, v_sigma, order, discret, im_ex, poissFlavor,
                     ecouple, ionFix, Te_start, Te_end, CL_type, ion_type,
                     MT_or_TR, dt, tfinal, BGK_type, beta, hydro_flag,
                     input_file_data_flag, dataFreq, outputDist, RHS_tol);

  strcat(input_path, inputFilename);
  printf("Opening input file %s\n", input_path);
  /*Open input file*/
  input_file = fopen(input_path, "r");

  if (input_file == NULL) {
    printf("Error - input file not found\n");
    exit(1);
  }

  *nspec = -1;

  /*Read input file*/
  while (strcmp(line, "Stop") != 0) {

    read_line(input_file, line);

    if (rank == 0)
      printf("%s\n", line);

    /*Mesh info*/
    if (strcmp(line, "nspec") == 0) {
      *nspec = read_int(input_file);
      if (rank == 0)
        printf("%d\n", *nspec);
    }

    if (strcmp(line, "dims") == 0) {
      *dims = read_int(input_file);
      if (rank == 0)
        printf("%d\n", *dims);
    }

    /*Number of physical space grid pts*/
    if (strcmp(line, "Nx") == 0) {
      *Nx = read_int(input_file);
      if (rank == 0)
        printf("%d\n", *Nx);
    }

    /*Length of physical domain in m*/
    if (strcmp(line, "Lx") == 0) {
      *Lx = read_double(input_file);
      if (rank == 0)
        printf("%g\n", *Lx);
    }

    if(strcmp(line, "Bcs") == 0){
      *bcs = read_int(input_file);
      if(rank == 0){
        printf("%d\n", *bcs);
      }
    }

    /*Number of velocity nodes in each dimension*/
    if (strcmp(line, "Nv") == 0) {
      *Nv = read_int(input_file);
      if (rank == 0)
        printf("%d\n", *Nv);
    }

    if (strcmp(line, "v_width") == 0) {
      *v_sigma = read_double(input_file);
      if (rank == 0)
        printf("%g\n", *v_sigma);
    }

    if (strcmp(line, "discret") == 0) {
      *discret = read_int(input_file);
      if (rank == 0)
        printf("%d\n", *discret);
    }

    /*Type of poisson solver*/
    if (strcmp(line, "poiss") == 0) {
      *poissFlavor = read_int(input_file);
      if (rank == 0)
        printf("%d\n", *poissFlavor);
    }

    /*Time-step in s*/
    if (strcmp(line, "Time_step") == 0) {
      *dt = read_double(input_file);
      if (rank == 0)
        printf("%g\n", *dt);
    }

    /*Number of time-steps in s*/
    if (strcmp(line, "Final_time") == 0) {
      *tfinal = read_double(input_file);
      if (rank == 0)
        printf("%g\n", *tfinal);
    }

    /*Order of accuracy of space/time discretization*/
    if (strcmp(line, "Space_order") == 0) {
      *order = read_int(input_file);
      if (rank == 0)
        printf("%d\n", *order);
    }

    /* Flag to use the hydro kinetic scheme*/
    if (strcmp(line, "Hydro_flag") == 0) {
      *hydro_flag = read_int(input_file);
      if (rank == 0)
        printf("%d\n", *hydro_flag);
    }

    //if eqflag >0, don't timestep, just run until the total entropy
    //change is below a certain tolerance, relative to the previous entropy.
    if (strcmp(line, "Equilibrate_flag") == 0) {
      *eqflag = read_int(input_file);
      if (rank == 0)
        printf("%d\n", *eqflag);
    }

    if (strcmp(line, "Equilibrate_rtol") == 0) {
      *eqrtol = read_double(input_file);
      if (rank == 0)
        printf("%g\n", *eqrtol);
    }

    /*Implicit solve (lagged) for the BGK operator*/
    if (strcmp(line, "Imp_exp") == 0) {
      *im_ex = read_int(input_file);
      if (rank == 0)
        printf("%d\n", *im_ex);
    }

    /*Species info*/

    /*species masses in kg*/
    if (strcmp(line, "mass") == 0) {
      if (*nspec == -1) {
        printf("ERROR: number of species not set\n");
        exit(37);
      }

      *m = malloc((*nspec) * sizeof(double));
      for (i = 0; i < (*nspec); i++) {
        (*m)[i] = read_double(input_file);
        if (rank == 0)
          printf("%g\n", (*m)[i]);
      }
    }

    /*species  ionizations*/
    if (strcmp(line, "Z") == 0) {
      if (*nspec == -1) {
        printf("ERROR: number of species not set\n");
        exit(37);
      }

      *Z = malloc((*nspec) * sizeof(double));
      for (i = 0; i < (*nspec); i++) {
        (*Z)[i] = read_double(input_file);
        if (rank == 0)
          printf("%g\n", (*Z)[i]);
      }
    }

    /*initial data info*/

    /*ICs for 0D problems*/

    // Number densities in cc
    if (strcmp(line, "n") == 0) {
      if (*nspec == -1) {
        printf("ERROR: number of species not set");
        exit(37);
      }

      *n = (double *)malloc((*nspec) * sizeof(double));
      for (i = 0; i < (*nspec); i++) {
        (*n)[i] = read_double(input_file);
        printf("%d %g\n", i, (*n)[i]);
      }
    }

    // x1 velocities in m/s
    if (strcmp(line, "v") == 0) {
      if (*nspec == -1) {
        printf("ERROR: number of species not set");
        exit(37);
      }

      *u = (double *)malloc((*nspec) * sizeof(double));
      for (i = 0; i < (*nspec); i++) {
        (*u)[i] = read_double(input_file);
        printf("%g\n", (*u)[i]);
      }
    }

    // Temperatures in eV
    if (strcmp(line, "T") == 0) {
      if (*nspec == -1) {
        printf("ERROR: number of species not set");
        exit(37);
      }

      *T = (double *)malloc((*nspec) * sizeof(double));
      for (i = 0; i < (*nspec); i++) {
        (*T)[i] = read_double(input_file);
        printf("%g\n", (*T)[i]);
      }
    }

    /*ICs for 1D problems*/
    if (strcmp(line, "Init_from_file") == 0) {
      *input_file_data_flag = 1;
      read_line(input_file, input_file_data_filename);
    }

    if (strcmp(line, "NumIntervals") == 0) {
      *numint = read_int(input_file);
      if (rank == 0)
        printf("%d\n", *numint);

      *intervalLimits = (double *)malloc((*numint) * sizeof(double));
      *ndens_int = (double *)malloc((*numint) * (*nspec) * sizeof(double));
      *velo_int = (double *)malloc((*numint) * (*nspec) * sizeof(double));
      *T_int = (double *)malloc((*numint) * (*nspec) * sizeof(double));

      // load data
      while (strcmp(line, "End_init") != 0) {
        read_line(input_file, line);
        if (rank == 0)
          printf("%s\n", line);

        if (strcmp(line, "Interval") == 0) {
          int_id = (int)read_int(input_file);
          int_id--;
          if (rank == 0)
            printf("%d\n", int_id);
        }

        if (strcmp(line, "x") == 0) {
          (*intervalLimits)[int_id] = read_double(input_file);
          if (rank == 0)
            printf("%g\n", (*intervalLimits)[int_id]);
        }

        if (strcmp(line, "n_i") == 0) {
          for (i = 0; i < (*nspec); i++) {
            (*ndens_int)[int_id * (*nspec) + i] = read_double(input_file);
            if (rank == 0)
              printf("%g\n", (*ndens_int)[int_id * (*nspec) + i]);
          }
        }

        if (strcmp(line, "v_i") == 0) {
          for (i = 0; i < (*nspec); i++) {
            (*velo_int)[int_id * (*nspec) + i] = read_double(input_file);
            if (rank == 0)
              printf("%g\n", (*velo_int)[int_id * (*nspec) + i]);
          }
        }

        if (strcmp(line, "T_i") == 0) {
          for (i = 0; i < (*nspec); i++) {
            (*T_int)[int_id * (*nspec) + i] = read_double(input_file);
            if (rank == 0)
              printf("%g\n", (*T_int)[int_id * (*nspec) + i]);
          }
        }
      }
    }

    /* collision information */

    // Multi-species BGK model to use
    // 0 - the Haack-Hauck-Murillo BGK
    // 1 - Greene BGK
    // 2 - NRL BGK
    if (strcmp(line, "BGKtype") == 0) {
      *BGK_type = read_int(input_file);
      if (rank == 0)
        printf("%d\n", *BGK_type);
    }

    // Free parameter for Greene model, if needed
    if (strcmp(line, "beta") == 0) {
      *beta = read_double(input_file);
    }

    // set to 1 for coupling with background electron field with specified (in
    // code, for now) eV
    // set to 2 to indicate that species 0 is electrons and to use B-Y cross
    // section
    // set to 3 for coupling with background electron field that matches mixture temp
    if (strcmp(line, "ecouple") == 0) {
      *ecouple = read_int(input_file);
      if (rank == 0)
        printf("%d\n", *ecouple);
    }

    if (strcmp(line, "ion_fix") == 0) {
      *ionFix = read_int(input_file);
      if (rank == 0)
        printf("%d\n", *ionFix);
    }

    if (strcmp(line, "Te_start") == 0) {
      *Te_start = read_double(input_file);
      if (rank == 0)
        printf("%g\n", *Te_start);
    }

    if (strcmp(line, "Te_end") == 0) {
      *Te_end = read_double(input_file);
      if (rank == 0)
        printf("%g\n", *Te_end);
    }

    // Identify coulomb log type, if needed
    // 0 is GMS
    // 1 is NRL
    if (strcmp(line, "Coulomb_type") == 0) {
      *CL_type = read_int(input_file);
      if (rank == 0)
        printf("%d\n", *CL_type);
    }

    // Chooses the source for ion-ion collision rates
    // 0 - Stanton-Murillo
    // 1 - NRL
    if (strcmp(line, "Ion_coll_type") == 0) {
      *ion_type = read_int(input_file);
      if (rank == 0)
        printf("%d\n", *ion_type);
    }

    // Chooses the flavor of coll rate for the given cross section
    // 0 - Momentum Transfer
    // 1 - Temperature relaxation
    if (strcmp(line, "Ion_coll_flavor") == 0) {
      *MT_or_TR = read_int(input_file);
      if (rank == 0)
        printf("%d\n", *MT_or_TR);
    }

    /*output files info*/

    /*Output file writing rate*/
    if (strcmp(line, "Data_writing_frequency") == 0) {
      *dataFreq = read_int(input_file);
      if (rank == 0)
        printf("%d\n", *dataFreq);
    }

    if (strcmp(line, "Dump_distro") == 0) {
      *outputDist = read_int(input_file);
      if (rank == 0)
        printf("%d\n", *outputDist);
    }

    if (strcmp(line, "RHS_tol") == 0)
      *RHS_tol = read_double(input_file);
  }

  if(rank == 0)
    printf("done with input file\n");
  fflush(stdout);

  fclose(input_file);
}

/*!
 *  This function sets the input parameters to their defualt values, and sets up
 * flags for a few if they are not set by the input file
 */
void set_default_values(int *Nx, double *Lx, int *bcs, int *Nv, double *v_sigma,
                        int *order, int *discret, int *im_ex, int *poissFlavor,
                        int *ecouple, int *ionFix, double *Te_start,
                        double *Te_end, int *CL_type, int *ion_type,
                        int *MT_or_TR, double *dt, double *tfinal,
                        int *BGK_type, double *beta, int *hydro_flag, int *eqflag, double *eqrtol,
                        int *input_file_data_flag, int *dataFreq, int *dumpDist,
                        double *RHS_tol) {

  /*Mesh info*/

  /*Number of physical space grid pts*/
  *Nx = 32;

  /*Length of physical domain in cm*/
  *Lx = 1000e-6; // 1000 nm;

  *bcs = 0;

  /*Number of velocity nodes in each dimension*/
  *Nv = 40;

  /* width of velocity domain in thermal velocities */
  *v_sigma = 6.0;

  /*Time-step in s*/
  *dt = 1e-16;

  /*Final time in s*/
  *tfinal = 1e-12;

  *Te_start = 100.0;

  *Te_end = 100.0;

  *discret = 1;

  *im_ex = 0;

  *poissFlavor = 22;

  *ecouple = 0;

  *ionFix = 0;

  *CL_type = 0;

  *ion_type = 0;

  *MT_or_TR = 0;

  *BGK_type = 0;

  *beta = 0;

  *hydro_flag = 0;
  *eqflag = 0;
  *eqrtol = 1e-8;


  *input_file_data_flag = 0;

  /*Order of accuracy of space/time discretization*/
  *order = 1;

  /*Output file writing rate*/
  *dataFreq = 1;

  *dumpDist = 0;

  *RHS_tol = 0.5;
}

/*!
 *  This function checks the return value of the fscanf function while reading
 * the input file
 */
void check_input(const int *flag) {
  if ((*flag) != 1) {
    printf("ERROR\nFlag: %d", (*flag));
    printf("\n%s\n", "ERROR: check_input");
    printf("\n%s\n", "Error while reading the input file!");
    exit(0);
  }
}

/*!
 * This function reads a string of characters from a file
 */
void read_line(FILE *file, char line[10000]) {

  fgets(line, 10000, file);
  if (line == NULL) {
    printf("\n%s\n", "ERROR: read_line");
    printf("\n%s\n", "Error while reading the input file!");
    exit(0);
  }
  size_t len = strlen(line);
  if (len > 0 && line[len - 1] == '\n')
    line[--len] = '\0';
}

/*!
 * This function reads an unsigned integer from a file
 */
size_t read_int(FILE *file) {
  int read;
  size_t int_value;

  read = fscanf(file, "%lu\n", &int_value);
  check_input(&read);

  return int_value;
}

/*!
 * This function reads a double from a file
 */
double read_double(FILE *file) {
  int read;
  double double_value;

  read = fscanf(file, "%lf\n", &double_value);
  check_input(&read);

  return double_value;
}

/*!
 * This function reads a string of characters from a file without moving on a
 * new line
 */
void read_line_no_adv(FILE *file, char line[80]) {
  int read;

  read = fscanf(file, "%s", line);
  check_input(&read);
}

/*!
 * This function reads an unsigned integer from a file without moving on a new
 * line
 */
size_t read_int_no_adv(FILE *file) {
  int read;
  size_t int_value;

  read = fscanf(file, "%lu", &int_value);
  check_input(&read);

  return int_value;
}

/*!
 * This function reads a double from a file without moving on a new line
 */
double read_double_no_adv(FILE *file) {
  int read;
  double double_value;

  read = fscanf(file, "%lf", &double_value);
  check_input(&read);

  return double_value;
}
