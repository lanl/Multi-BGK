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

  load_distributions_inhomog(f, input_filename, 1, rank);

  printf("Loaded distributions\n");

  MPI_Finalize();

  return 0;
}
