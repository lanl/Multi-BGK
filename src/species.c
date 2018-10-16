#include "species.h"
#include <stdio.h>
#include <string.h>

// need to add line reading stuff

// For reference
// typedef struct {
//  size_t  id;              /*!< identifier (macroscopic)*/
//  size_t  num_levels;       /*!< number of internal energy levels*/
//  size_t *lev_id;          /*!< vector storing the global id of the internal
//  energy levels*/ double  Rgas;            /*!< gas constants [J/kg]*/ double
//  mass;            /*!< mass [kg]*/ double  mm;              /*!< molecular
//  mass [kg/mol]*/ double  d_ref;           /*!< reference diameter [m]*/
//  double  T_ref;           /*!< reference temperature [K]*/
//  double  mu_ref;          /*!< reference dynamic viscosity [N*s/m^2]*/
//  double  omega;           /*!< viscosity law temperature exponent*/
//  double  E0;              /*!< formation energy [J]*/
//  double *Ei;              /*!< vector storing internal energy levels [J]*/
//  double *gi;              /*!< vector storing internal energy level
//  degeneracies*/ char    name[80];        /*!< name*/
//} species;

// input function has read in the number of species and their names from the
// input file

void set_spec_default_values(species *mixture, int num_species) {
  int i;

  for (i = 0; i < num_species; i++) {
    mixture[i].id = i;
    mixture[i].num_levels = 0;
    mixture[i].Rgas = 1.;
    mixture[i].mass = 1.;
    mixture[i].mm = 1.;
    mixture[i].d_ref = 2.;
    mixture[i].T_ref = 1.;
    mixture[i].mu_ref = 1.;
    mixture[i].E0 = 0.;
    strcpy(mixture[i].name, "default");
  }
}

void load_and_allocate_spec(species **mixture, int num_species,
                            char **species_list) {
  size_t i;
  char line[80] = {"dummy"};
  char input_path[100] = {"./input/species/"};
  char name[80] = {"dummy"};
  FILE *input_file;

  int num_levels;
  double mass, E0, d_ref, T_ref, omega, mu_ref;

  *mixture = malloc(num_species * sizeof(species));

  set_spec_default_values(*mixture, num_species);

  printf("num species %d\n", num_species);
  for (i = 0; i < num_species; i++)
    printf("%s\n", species_list[i]);

  fflush(stdout);

  for (i = 0; i < num_species; i++) {
    if (strcmp(species_list[i], "default") == 0)
      continue;
    strcat(input_path, species_list[i]);
    printf("Opening species input file %s\n", input_path);
    /*Open input file*/
    input_file = fopen(input_path, "r");

    if (input_file == NULL) {
      printf("Error - input file not found\n");
      exit(1);
    }

    /*Read input file*/
    while (strcmp(line, "Stop") != 0) {
      read_line(input_file, line);

      if (strcmp(line, "name") == 0) {
        read_line(input_file, name);
        strcpy((*mixture)[i].name, name);
        printf("name: %s\n", (*mixture)[i].name);
      }

      /*Number of internal energy levels*/
      if (strcmp(line, "num_en_levels") == 0) {

        num_levels = read_int(input_file);
        (*mixture)[i].num_levels = num_levels;
        printf("num_levels:%d\n", num_levels);

        /*Read energy levels and degeneracies*/
        /*mixture[i].lev_id = malloc(max(1,num_levels) * sizeof(size_t));
        if ((num_levels > 0)) {
          mixture[i].Ei = malloc(num_levels * sizeof(double));
          mixture[i].gi = malloc(num_levels * sizeof(double));
          for (j = 0;j < num_levels;j++) {
            g = read_double_no_adv(input_file);
            E = read_double(input_file);
            mixture[i].gi[j] = g;
            mixture[i].Ei[j] = E;
          }
          }*/
      }

      /*Mass*/
      if (strcmp(line, "mass") == 0) {
        mass = read_double(input_file);
        // mm   = mass*NA;
        // Rgas = RU/mm;
        (*mixture)[i].mass = mass;
        printf("mass: %g\n", (*mixture)[i].mass);
        // mixture[s].mm   = mm;
        // mixture[s].Rgas = Rgas;
      }

      /*Formation energy*/
      if (strcmp(line, "En_for") == 0) {
        E0 = read_double(input_file);
        (*mixture)[i].E0 = E0;
        printf("Formation energy:%g\n", (*mixture)[i].E0);
      }

      /*Reference diameter*/
      if (strcmp(line, "d_ref") == 0) {
        d_ref = read_double(input_file);
        (*mixture)[i].d_ref = d_ref;
      }

      /*Reference temperature*/
      if (strcmp(line, "T_ref") == 0) {
        T_ref = read_double(input_file);
        (*mixture)[i].T_ref = T_ref;
      }

      /*Reference dynamic viscosity*/
      if (strcmp(line, "mu_ref") == 0) {
        mu_ref = read_double(input_file);
        (*mixture)[i].mu_ref = mu_ref;
      }

      /*Viscosity law exponent*/
      if (strcmp(line, "omega") == 0) {
        omega = read_double(input_file);
        (*mixture)[i].omega = omega;
      }
      fflush(stdout);
    }

    /*Close species file*/
    fclose(input_file);

    /*Re-initialize the file name and the file line being read*/
    strcpy(input_path, "./input/species/");
    strcpy(line, "dummy");
  }
}
