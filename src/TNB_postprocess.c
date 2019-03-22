/*

Post-processor for TNB. This assumes that you ran the main code with

mpirun -n <num_ranks> ./exec/MultiBGK_ <input_filename>

and Dump_Distro set to 1 for a 1D run

to build do

mpicc src/TNB_postprocess.c src/io.c src/TNB.c -o exec/postProc_

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

  io_init_inhomog(Nx_rank, Nv, Nspec, NULL);

  load_grid_inhomog(Lv, &Nx_rank, input_filename, rank);

  double ***f = malloc(Nx_rank * sizeof(double **));
  for (i = 0; i < Nx_rank; i++) {
    f[i] = malloc(Nspec * sizeof(double *));
    for (s = 0; s < Nspec; s++) {
      f[i][s] = malloc(Nv * Nv * Nv * sizeof(double));
    }
  }

  io_init_inhomog(Nx_rank, Nv, Nspec, NULL);

  load_distributions_inhomog(f, input_filename, step, rank);

  printf("Loaded distributions\n");

  // Now do the TNB...

  // TNB will need

  // Nv
  // Masses of species
  // velocity grids for each species
  // integration weights for each species
  // The distribution functions

  MPI_Finalize();

  return 0;
}
