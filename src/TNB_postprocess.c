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
#include "input.h"
#include "io.h"

int main(int argc, char **argv) {

  MPI_Init(&argc, &argv);

  int rank, numRanks;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
  MPI_Status status;

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

  char outputFileBuffer[500];
  sprintf(outputFileBuffer, "Data/TNB_DT_step%d.dat", step);

  double *valueBuffer = malloc((Nx_rank + 1) * sizeof(double));

  FILE *F_dt;

  if (rank == 0) {
    F_dt = fopen(outputFileBuffer, "w");
  }

  for (int xval = 0; xval < Nx_rank; xval++) {
    for (int spec1 = 0; spec1 < Nspec; spec1++) {
      for (int spec2 = 0; spec2 < Nspec; spec2++) {
        // printf("xval %d, spec1 %d, spec2 %d\n", xval, spec1, spec2);
        double mu = m[spec1] * m[spec2] / (m[spec1] + m[spec2]);

        if (mu < 2.e-24 && mu > 1.8e-24) { // Just get DT reaction

          double R_BGK_DT = GetReactivity_dt(mu, f[xval][spec1], f[xval][spec2],
                                             spec1, spec2);

          valueBuffer[xval] = R_BGK_DT;
        }
      }
    }
  }

  // Now print to file

  // Get the x values
  char xFileBuffer[500];
  sprintf(xFileBuffer, "./Data/%s_x", input_filename);
  printf("%s\n", xFileBuffer);
  FILE *x_file;
  x_file = fopen(xFileBuffer, "r");
  double x_point = 0;

  int *rankNx = malloc(sizeof(int));

  if (rank == 0) {
    for (int xval = 0; xval < Nx_rank; xval++) {
      x_point = read_double(x_file);
      fprintf(F_dt, "%10.6e %10.6e\n", x_point, valueBuffer[xval]);
    }

    // Gather from other ranks
    for (int recvrank = 1; recvrank < numRanks; recvrank++) {
      MPI_Recv(rankNx, 1, MPI_INT, recvrank, recvrank, MPI_COMM_WORLD, &status);
      MPI_Recv(valueBuffer, *rankNx, MPI_DOUBLE, recvrank, 100 + recvrank,
               MPI_COMM_WORLD, &status);
      for (int xval = 0; xval < *rankNx; xval++) {
        x_point = read_double(x_file);
        fprintf(F_dt, "%10.6e %10.6e\n", x_point, valueBuffer[xval]);
      }
    }

  } else {
    MPI_Send(&Nx_rank, 1, MPI_INT, 0, rank, MPI_COMM_WORLD);
    MPI_Send(valueBuffer, Nx_rank, MPI_DOUBLE, 0, 100 + rank, MPI_COMM_WORLD);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  fclose(F_dt);

  MPI_Finalize();

  return 0;
}
