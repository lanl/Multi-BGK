#include "io.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int Nx;
static int Nv;
static int nspec;
static double **c;

void io_init_homog(int numV, int numS, double **velos) {

  Nv = numV;
  nspec = numS;
  c = velos;
}

void io_init_inhomog(int numX, int numV, int numS, double **velos) {
  Nx = numX;
  Nv = numV;
  nspec = numS;
  c = velos;
}

void store_distributions_homog(double **f, double t, int step, char *fileName) {
  int s;
  char name_buffer[100], grid_buffer[100];
  FILE *fid_store, *fid_grids;

  sprintf(grid_buffer, "Data/%s_gridinfo.dat", fileName);
  fid_grids = fopen(grid_buffer, "w");

  for (s = 0; s < nspec; s++) {
    if (step >= 0)
      sprintf(name_buffer, "Data/%s_spec%d_step%d.dat", fileName, s, step);
    else
      sprintf(name_buffer, "Data/%s_spec%d.dat", fileName, s);

    fid_store = fopen(name_buffer, "w");

    fwrite(f[s], sizeof(double), Nv * Nv * Nv, fid_store);

    fprintf(fid_grids, "%d %d %g\n", s, Nv, -c[s][0]);

    fclose(fid_store);
  }

  fprintf(fid_grids, "%g\n", t);
  if (step < 0)
    fprintf(fid_grids, "%d\n", -1 * step);

  printf("Stored time %g\n", t);

  fclose(fid_grids);
}

void load_distributions_homog(double **f, char *fileName) {
  int s, readflag;
  char name_buffer[100];
  FILE *fid_load;

  for (s = 0; s < nspec; s++) {
    sprintf(name_buffer, "Data/%s_spec%d.dat", fileName, s);

    fid_load = fopen(name_buffer, "r");

    if (fid_load == NULL) {
      printf("Error: unable to open distribution function file %s\n",
             name_buffer);
      exit(1);
    }

    readflag = fread(f[s], sizeof(double), Nv * Nv * Nv, fid_load);
    if (readflag != Nv * Nv * Nv) {
      printf("Error reading distribution function data\n");
      exit(1);
    }

    fclose(fid_load);
  }
}

void store_grid(char *fileName) {
  int s;
  char grid_buffer[100];
  FILE *fid_grids;

  sprintf(grid_buffer, "Data/%s_gridinfo.dat", fileName);
  fid_grids = fopen(grid_buffer, "w");

  for (s = 0; s < nspec; s++) {
    fprintf(fid_grids, "%d %d %g\n", s, Nv, -c[s][0]);
  }
  fclose(fid_grids);
}

// Not looking to just dump binary data here - need something that I can load in
// via python.
void store_distributions_inhomog(double ***f, char *fileName, int t) {
  int s, xj;
  char name_buffer[100];
  FILE *fid_store;

  for (s = 0; s < nspec; s++) {
    sprintf(name_buffer, "Data/%s_spec%d.step%d.dat", fileName, s, t);

    fid_store = fopen(name_buffer, "w");

    for (xj = 0; xj < Nx; xj++)
      fwrite(f[xj][s], sizeof(double), Nv * Nv * Nv, fid_store);
    fclose(fid_store);
  }
}

void load_taus_homog(double **nu, char *filename) {
  int s, t, readflag;
  double tmp;
  char name_buffer[100];
  FILE *fid_load;

  sprintf(name_buffer, "Data/%s_tau.dat", filename);

  fid_load = fopen(name_buffer, "r");

  if (fid_load == NULL) {
    printf("Error: could not locate tau datafile %s\n", name_buffer);
    exit(1);
  }

  for (s = 0; s < nspec; s++) {
    printf("Loading the taus from %s for species %d\n", name_buffer, s);
    for (t = 0; t < nspec; t++) {
      readflag = fscanf(fid_load, "%lf\n", &tmp);
      nu[s][t] = 1.0 / tmp;
      printf("i: %d j: %d nu: %g\n", s, t, 1.0 / nu[s][t]);
    }
  }

  fclose(fid_load);
}

void load_grid_restart(double *Lv, double *t, int *nT, char *fileName) {
  int s, readflag;
  char name_buffer[100];
  char line[100];
  FILE *fid_load;

  int spec, num_v, nt_val;
  double Lv_val, t_val;

  sprintf(name_buffer, "Data/%s_gridinfo.dat", fileName);

  fid_load = fopen(name_buffer, "r");

  if (fid_load == NULL) {
    printf("Error: unable to open distribution function file %s\n",
           name_buffer);
    exit(1);
  }

  // printf("Opening %s to load grid information\n",name_buffer);

  for (s = 0; s < nspec; s++) {
    fgets(line, 100, fid_load);
    // printf("%s",line);

    readflag = sscanf(line, "%d %d %lf\n", &spec, &num_v, &Lv_val);
    Lv[s] = Lv_val;
    printf("Species %d velo limit is %g\n", s, Lv[s]);
  }
  if (fgets(line, 100, fid_load) != NULL) {
    readflag = sscanf(line, "%lf\n", &t_val);
  } else {
    printf("Error - specifiy initial time in second to last line of the "
           "gridinfo file\n");
    exit(1);
  }
  if (fgets(line, 100, fid_load) != NULL) {
    readflag = sscanf(line, "%d\n", &nt_val);
  } else {
    printf("Error - specifiy initial time step in last line of the gridinfo "
           "file\n");
    exit(1);
  }

  *t = t_val;
  *nT = nt_val;

  fclose(fid_load);
}
