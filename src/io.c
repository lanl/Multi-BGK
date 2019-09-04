#include "io.h"
#include <stdio.h>
#include <stdlib.h>

static int Nx;
static int Nv;
static int nspec;
static int order;
static double **c;
static double *m;

void io_init_homog(int numV, int numS, double **velos) {

  Nv = numV;
  nspec = numS;
  c = velos;
}

void io_init_inhomog(int numX, int numV, int numS, int o, double **velos,
                     double *mass) {
  Nx = numX;
  Nv = numV;
  nspec = numS;
  order = o;
  c = velos;
  m = mass;
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

void store_grid_inhomog(char *fileName, int rank) {
  int s;
  char grid_buffer[100];
  FILE *fid_grids;

  sprintf(grid_buffer, "Data/%s_gridinfo_rank%d.dat", fileName, rank);
  fid_grids = fopen(grid_buffer, "w");

  fprintf(fid_grids, "%d\n", Nx);

  for (s = 0; s < nspec; s++) {
    fprintf(fid_grids, "%d %d %lf %lg\n", s, Nv, -c[s][0], m[s]);
  }
  fclose(fid_grids);
}

void store_distributions_inhomog(double ***f, char *fileName, int nT, double t,
                                 int rank) {
  int s, xj;
  char name_buffer[10000];
  char time_buffer[10000];
  FILE *fid_store;
  FILE *fid_time;

  for (s = 0; s < nspec; s++) {
    sprintf(name_buffer, "Data/%s_spec%d.step%d_rank%d.dat", fileName, s, nT,
            rank);

    fid_store = fopen(name_buffer, "w");

    for (xj = 0; xj < Nx; xj++)
      fwrite(f[xj + order][s], sizeof(double), Nv * Nv * Nv, fid_store);
    fclose(fid_store);
  }

  if (rank == 0) {
    sprintf(time_buffer, "Data/%s_lasttimestep.dat", fileName);

    fid_time = fopen(time_buffer, "w");
    fprintf(fid_time, "%d %lg\n", nT, t);
    fclose(fid_time);
  }
}

void load_distributions_inhomog(double ***f, char *fileName, int nT, int rank) {
  int s, xj;
  char name_buffer[100];
  FILE *fid_load;
  int readflag;

  for (s = 0; s < nspec; s++) {
    sprintf(name_buffer, "Data/%s_spec%d.step%d_rank%d.dat", fileName, s, nT,
            rank);

    // printf("Loading %s\n",name_buffer);

    fid_load = fopen(name_buffer, "r");

    if (fid_load == NULL) {
      printf("Cannot find distribution file %s\n", name_buffer);
      exit(37);
    }

    for (xj = 0; xj < Nx; xj++) {
      readflag =
          fread(f[xj + order][s], sizeof(double), Nv * Nv * Nv, fid_load);
      if (readflag != Nv * Nv * Nv) {
        printf("Error reading distribution function data\n");
        exit(2);
      }
    }
    fclose(fid_load);

    // printf("Completed %s\n",name_buffer);
  }
}

void load_time_inhomog(char *fileName, int *nT, double *t) {

  FILE *fid_time;
  char filenamebuffer[1000];
  char line[100];
  int readflag;

  sprintf(filenamebuffer, "Data/%s_lasttimestep.dat", fileName);

  fid_time = fopen(filenamebuffer, "r");

  if (fid_time == NULL) {
    printf("Cannot find time file %s\n", filenamebuffer);
    exit(37);
  }

  fgets(line, 100, fid_time);

  readflag = sscanf(line, "%d %lg\n", nT, t);

  printf("Loading time from %s obtained %d %lg\n", filenamebuffer, *nT, *t);
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

// This is geared towards TNB postprocessing, not restarts

void load_grid_inhomog(double *Lv, int *Nx, double *mass, char *fileName,
                       int rank) {
  int s, readflag;
  char name_buffer[100];
  char line[100];
  FILE *fid_load;

  int spec, num_v;
  double Lv_val, m_val;

  sprintf(name_buffer, "Data/%s_gridinfo_rank%d.dat", fileName, rank);

  fid_load = fopen(name_buffer, "r");

  if (fid_load == NULL) {
    printf("Error: unable to open distribution function file %s\n",
           name_buffer);
    exit(1);
  }

  // printf("Opening %s to load grid information\n",name_buffer);

  fgets(line, 100, fid_load);
  readflag = sscanf(line, "%d\n", Nx);

  for (s = 0; s < nspec; s++) {
    fgets(line, 100, fid_load);
    // printf("%s",line);

    readflag = sscanf(line, "%d %d %lf %lg\n", &spec, &num_v, &Lv_val, &m_val);
    Lv[s] = Lv_val;
    mass[s] = m_val;
    // printf("Species %d velo limit is %g mass is %g\n", s, Lv_val, m_val);
  }

  fclose(fid_load);
}
