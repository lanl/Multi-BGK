#include "io.h"
#include "units/unit_data.c"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef ALDR_ON
#include <alInterface.h>
#include <sqlite3.h>
#endif

#ifndef NDENS_TOL
#define NDENS_TOL 1e10
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
void io_init_db(char *filename) { db = initDB(0, filename); }

void io_close_db() { closeDB(db); }

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

void request_aldr_single(double *n, double *T, double *Z, char *tag,
                         char *dbfile, double **D_ij) {

  bgk_request_t input;
  bgk_result_t output;

  // Initialize input request
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
  double ntot = 0.0;

  input.temperature = Tmix;

  for (int sp = 0; sp < nspec; sp++) {
    if (n[sp] > NDENS_TOL)
      input.density[sp] = n[sp];
    else
      input.density[sp] = 0.0;

    input.charges[sp] = Z[sp];
  }

  // Calculate mixture T
  for (int sp = 0; sp < nspec; sp++) {
    if (n[sp] > NDENS_TOL) {
      Tmix += n[sp] * T[sp];
      ntot += n[sp];
    }
  }
  if (ntot == 0.0) {
    printf("Error - zero number density\n");
    exit(37);
  }

  Tmix /= ntot;

  output = bgk_req_single(input, 0, tag, db);

  D_ij[0][0] = output.diffusionCoefficient[0];
  D_ij[0][1] = output.diffusionCoefficient[1];
  D_ij[0][2] = output.diffusionCoefficient[2];
  D_ij[0][2] = output.diffusionCoefficient[3];
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

/*
Calculates the species coupling parameters at the given density, temperature,
and charge state Uses formulas rho_tot = \sum_j Z_j e n_j a_i = (3 Z_i e / (4
\pi \rho_tot)) ** (1/3) Gamma_i = (Z_i e)^2 / a_i T
*/
double get_max_coupling(double *n, double T, double *Z) {

  double rho_tot = 0.0;
  double Gamma_max = 0.0;
  double Gamma_i = 0.0;
  double a_i = 0.0;

  for (int sp = 0; sp < nspec; sp++)
    rho_tot += Z[sp] * n[sp];

  for (int sp = 0; sp < nspec; sp++) {
    a_i = pow(3.0 * Z[sp] / (4.0 * M_PI * rho_tot), 1.0 / 3.0);
    Gamma_i = Z[sp] * Z[sp] * E_02_CGS / a_i / T;
    Gamma_max = Gamma_i > Gamma_max ? Gamma_i : Gamma_max;
  }

  return Gamma_max;
}

// Sends a whole buncha results to the glue code
void request_aldr_batch(double **n, double **T, double **Z, char *tag,
                        char *dbfile, double ***D_ij) {

  bgk_request_t *input_list;
  bgk_result_t *output_list;
  // bgk_result_t *output_list = malloc(Nx * sizeof(bgk_result_t));

  unsigned strongly_coupled_cells_count = 0;
  unsigned strongly_coupled_cells[Nx];

  // Determine the number of strongly coupled cells
  for (unsigned x_node = 0; x_node < Nx; ++x_node) {

    double Tmix = 0.0;
    double ntot = 0.0;

    // Calculate mixture T
    for (int sp = 0; sp < nspec; sp++) {
      if (n[x_node][sp] > NDENS_TOL) {
        Tmix += n[x_node][sp] * T[x_node][sp];
        ntot += n[x_node][sp];
      }
    }

    Tmix /= ntot;

    double Gamma_cell_max = get_max_coupling(n[x_node], Tmix, Z[x_node]);
    if (Gamma_cell_max > 0.1) {
      strongly_coupled_cells_count++;
      strongly_coupled_cells[x_node] = 1;
    } else
      strongly_coupled_cells[x_node] = 0;
  }

  input_list = malloc(strongly_coupled_cells_count * sizeof(bgk_request_t));

  // Initialize input list to ensure zeros go in for unused species
  for (unsigned i = 0; i < strongly_coupled_cells_count; i++) {
    for (unsigned sp = 0; sp < 4; sp++) {
      input_list[i].density[sp] = 0.0;
      input_list[i].charges[sp] = 0.0;
    }
    input_list[i].temperature = 0.0;
  }

  unsigned input_list_count = 0;
  for (unsigned x_node = 0; x_node < Nx; ++x_node) {

    if (strongly_coupled_cells[x_node]) {
      double Tmix = 0.0;
      double ntot = 0.0;

      for (int sp = 0; sp < nspec; sp++) {
        if (n[x_node][sp] > NDENS_TOL) {
          input_list[input_list_count].density[sp] = n[x_node][sp];
        } else {
          input_list[input_list_count].density[sp] = 0.0;
        }
        input_list[input_list_count].charges[sp] = Z[x_node][sp];
      }

      // Calculate mixture T
      for (int sp = 0; sp < nspec; sp++) {
        if (n[x_node][sp] > NDENS_TOL) {
          Tmix += n[x_node][sp] * T[x_node][sp];
          ntot += n[x_node][sp];
        }
      }

      Tmix /= ntot;

      input_list[input_list_count].temperature = Tmix;
      input_list_count++;
    }
  }

  output_list =
      bgk_req_batch(input_list, strongly_coupled_cells_count, 0, tag, db);

  // Store the results
  int output_list_index = 0;
  for (unsigned x_node = 0; x_node < Nx; ++x_node) {
    if (strongly_coupled_cells[x_node]) {
      D_ij[x_node][0][0] =
          output_list[output_list_index].diffusionCoefficient[0];
      D_ij[x_node][0][1] =
          output_list[output_list_index].diffusionCoefficient[1];
      D_ij[x_node][0][2] =
          output_list[output_list_index].diffusionCoefficient[2];
      D_ij[x_node][0][3] =
          output_list[output_list_index].diffusionCoefficient[3];
      D_ij[x_node][1][1] =
          output_list[output_list_index].diffusionCoefficient[4];
      D_ij[x_node][1][2] =
          output_list[output_list_index].diffusionCoefficient[5];
      D_ij[x_node][1][3] =
          output_list[output_list_index].diffusionCoefficient[6];
      D_ij[x_node][2][2] =
          output_list[output_list_index].diffusionCoefficient[7];
      D_ij[x_node][2][3] =
          output_list[output_list_index].diffusionCoefficient[8];
      D_ij[x_node][3][3] =
          output_list[output_list_index].diffusionCoefficient[9];

      // Symmetric components
      D_ij[x_node][1][0] = D_ij[x_node][0][1];
      D_ij[x_node][2][0] = D_ij[x_node][0][2];
      D_ij[x_node][3][0] = D_ij[x_node][0][3];
      D_ij[x_node][2][1] = D_ij[x_node][1][2];
      D_ij[x_node][3][1] = D_ij[x_node][1][3];
      D_ij[x_node][3][2] = D_ij[x_node][2][3];

      output_list_index++;

      printf("l: %d ", x_node);
      for (int i = 0; i < 10; i++)
        printf("D[%d]: %g ", i, output_list[x_node].diffusionCoefficient[i]);
      printf("\n");
    } else { // USE SM, setting this to -1 is the flag
      D_ij[x_node][0][0] = -1;
      printf("l: %d ", x_node);
      printf("Is using SM values\n");
    }
  }

  free(input_list);
  free(output_list);
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
