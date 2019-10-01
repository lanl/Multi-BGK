#include "io.h"
#include <stdio.h>
#include <stdlib.h>

#ifdef ALDR_ON
#include <alInterface.h>
#include <sqlite3.h>
#endif

static int Nx;
static int Nv;
static int nspec;
static double **c;

#ifdef ALDR_ON
static sqlite3 *db;
#endif

void io_init_homog(int numV, int numS, double **velos) {

  Nv = numV;
  nspec = numS;
  c = velos;

}

#ifdef ALDR_ON
void io_init_db(char *filename) {
  db = initDB(0, filename);  
}

void io_close_db() {
  closeDB(db);
}

#endif

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

void load_diffusion_homog(double **Dij, char *filename) {
  int s, s2, readflag;
  double tmp;
  char name_buffer[100];
  FILE *fid_load;

  sprintf(name_buffer, "Data/%s_Dij.dat", filename);

  fid_load = fopen(name_buffer, "r");

  if (fid_load == NULL) {
    printf("Error: could not locate diffusion coefficient datafile %s\n",
           name_buffer);
    exit(1);
  }

  for (s = 0; s < nspec; s++) {
    for (s2 = 0; s2 < nspec; s2++) {
      // printf("Loading the Dij from %s for species pair %d %d\n", name_buffer,
      // s, s2);
      readflag = fscanf(fid_load, "%lf\n", &tmp);
      Dij[s][s2] = tmp;
      // printf("i: %d j: %d Dij: %g\n", s, s2, Dij[s][s2]);
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

#ifdef ALDR_ON

void request_aldr(double *n, double *T, double *Z, char *tag, char *dbfile, double **D_ij) {

  bgk_request_t input;
  bgk_result_t output;

  //Initialize input request
  input.density[0] = 0.0;
  input.density[1] = 0.0;
  input.density[2] = 0.0;
  input.density[3] = 0.0;
  input.charges[0] = 0.0;
  input.charges[1] = 0.0;
  input.charges[2] = 0.0;
  input.charges[3] = 0.0;
  input.temperature = 0.0;

  output.viscosity = 0.0;
  output.thermalConductivity = 0.0;
  output.diffusionCoefficient[0] = 0.0;
  output.diffusionCoefficient[1] = 0.0;
  output.diffusionCoefficient[2] = 0.0;
  output.diffusionCoefficient[3] = 0.0;
  output.diffusionCoefficient[4] = 0.0;
  output.diffusionCoefficient[5] = 0.0;
  output.diffusionCoefficient[6] = 0.0;
  output.diffusionCoefficient[7] = 0.0;
  output.diffusionCoefficient[8] = 0.0;
  output.diffusionCoefficient[9] = 0.0;

  double Tmix = 0.0;
  double ntot = 0.0;;

  //Calculate mixture T
  for(int sp=0; sp < nspec; sp++) {
    Tmix += n[sp]*T[sp];
    ntot += n[sp];
  }
  Tmix /= ntot;
  
  input.temperature = Tmix;
  
  for(int sp=0; sp < nspec; sp++) {
    input.density[sp] = n[sp];
    input.charges[sp] = Z[sp];
  }

  output = bgk_req_single(input, 0, tag, db);
  
  D_ij[0][0] = output.diffusionCoefficient[0];
  D_ij[0][1] = output.diffusionCoefficient[1];
  D_ij[0][2] = output.diffusionCoefficient[2];
  D_ij[0][3] = output.diffusionCoefficient[3];
  D_ij[1][1] = output.diffusionCoefficient[4];
  D_ij[1][2] = output.diffusionCoefficient[5];
  D_ij[1][3] = output.diffusionCoefficient[6];
  D_ij[2][2] = output.diffusionCoefficient[7];
  D_ij[2][3] = output.diffusionCoefficient[8];
  D_ij[3][3] = output.diffusionCoefficient[9];

  // Symmetric components
  D_ij[1][0] = D_ij[0][1];
  D_ij[2][0] = D_ij[0][2];
  D_ij[3][0] = D_ij[0][3];
  D_ij[2][1] = D_ij[1][2];
  D_ij[3][1] = D_ij[1][3];
  D_ij[3][2] = D_ij[2][3];
}

void test_aldr() {

  bgk_request_t input;
  bgk_result_t output;

  char tag[100] = "dummy";
  char dbfile[100] = "testDB.db";

  sqlite3 *db_test;

  printf("testing ALDR linking\n");
  fflush(stdout);

  db_test = initDB(0, dbfile);

  input.temperature = 1.0;
  for (int i = 0; i < 4; i++) {
    input.density[i] = 1.0;
    input.charges[i] = 1.0;
  }

  output = bgk_req_single(input, 0, tag, db_test);

  printf("Output: viscosity %g, TC %g D11 %g\n", output.viscosity,
         output.thermalConductivity, output.diffusionCoefficient[0]);

  closeDB(db);
}

#endif
