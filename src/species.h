#ifndef _SPECIES_H
#define _SPECIES_H

#include <stdlib.h>

/*!  \file  species.h
 *   \brief This file provides the definition of the (macroscopic) species type
 */

typedef struct {
  size_t id;         /*!< identifier (macroscopic)*/
  size_t num_levels; /*!< number of internal energy levels*/
  size_t
      *lev_id; /*!< vector storing the global id of the internal energy levels*/
  double Rgas; /*!< gas constants [J/kg]*/
  double mass; /*!< mass [kg]*/
  double mm;   /*!< molecular mass [kg/mol]*/
  double d_ref;  /*!< reference diameter [m]*/
  double T_ref;  /*!< reference temperature [K]*/
  double mu_ref; /*!< reference dynamic viscosity [N*s/m^2]*/
  double omega;  /*!< viscosity law temperature exponent*/
  double E0;     /*!< formation energy [J]*/
  double *Ei;    /*!< vector storing internal energy levels [J]*/
  double *gi;    /*!< vector storing internal energy level degeneracies*/
  char name[80]; /*!< name*/
} species;

void load_and_allocate_spec(species **mixture, int num_species,
                            char **species_list);

#endif
