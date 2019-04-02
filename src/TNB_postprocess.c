/*

Post-processor for TNB. This assumes that you ran the main code with

mpirun -n <num_ranks> ./exec/MultiBGK_ <input_filename>

and Dump_Distro set to 1 for a 1D run

to build do

make postProc

to process this data run

mpirun -n <num_ranks> ./exec/postProc_ <input_filename> <nspec> <Nv> <step>

Where
num_ranks, input_filename, nspec, and Nv are the same as in the original run.
Select step to do the analysis for the specific timestep

Note: it will automatically find the number of x points on each rank, this is
stored by the code

*/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "TNB.h"
#include "io.h"

int main(int argc, char **argv) {

  MPI_Init(&argc, &argv);

  int rank, numRanks;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numRanks);

  char input_filename[100];

  strcpy(input_filename, argv[1]);
  int Nspec = atoi(argv[2]);
  int Nv = atoi(argv[3]);
  int step = atoi(argv[4]);

  int s, i;

  int Nx_rank = 0;
  double *Lv = malloc(Nspec * sizeof(double));
  double *m = malloc(Nspec * sizeof(double));

  io_init_inhomog(Nx_rank, Nv, Nspec, NULL, NULL);

  // printf("Io init the first time\n");

  load_grid_inhomog(Lv, &Nx_rank, m, input_filename, rank);

  // printf("Loaded inhomog grid\n");

  double ***f = malloc(Nx_rank * sizeof(double **));
  for (i = 0; i < Nx_rank; i++) {
    f[i] = malloc(Nspec * sizeof(double *));
    for (s = 0; s < Nspec; s++) {
      f[i][s] = malloc(Nv * Nv * Nv * sizeof(double));
    }
  }

  // Just setup for loading distributions.
  io_init_inhomog(Nx_rank, Nv, Nspec, NULL, NULL);

  load_distributions_inhomog(f, input_filename, step, rank);

  printf("Loaded distributions\n");

  // Now do the TNB...

  // TNB will need

  // Nv
  // Masses of species
  // velocity grids for each species
  // integration weights for each species
  // The distribution functions

  // So first we must build c and wts from the grid knowledge

  double **c = malloc(Nspec * sizeof(double *));
  double **wts = malloc(Nspec * sizeof(double *));

  for (int spec = 0; spec < Nspec; spec++) {
    c[spec] = malloc(Nv * sizeof(double));
    wts[spec] = malloc(Nv * sizeof(double));
    double dv = 2.0 * Lv[spec] / (Nv - 1.0);
    for (int velo = 0; velo < Nv; velo++) {
      c[spec][velo] = -Lv[spec] + dv * velo;
      wts[spec][velo] = dv;
    }
    wts[spec][0] *= 0.5;
    wts[spec][Nv - 1] *= 0.5;
  }

  initializeTNB(Nv, c, wts);

  char outputFileBuffer[50];
  sprintf(outputFileBuffer, "Data/TNB_DT_rank%d.dat", rank);

  FILE *F_dt = fopen(outputFileBuffer, "w");

  for (int xval = 0; xval < Nx_rank; xval++) {
    for (int spec1 = 0; spec1 < Nspec; spec1++) {
      for (int spec2 = 0; spec2 < Nspec; spec2++) {
        // printf("xval %d, spec1 %d, spec2 %d\n", xval, spec1, spec2);
        double mu = m[spec1] * m[spec2] / (m[spec1] + m[spec2]);
        double R_BGK_DT =
            GetReactivity_dt(mu, f[xval][spec1], f[xval][spec2], spec1, spec2);
        if (R_BGK_DT > 0)
          fprintf(F_dt, "%d %10.6e\n", xval, R_BGK_DT);
      }
    }
  }

  fclose(F_dt);

  MPI_Finalize();

  return 0;
}
