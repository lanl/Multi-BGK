#include "BGK.h"
#include <stdio.h>
#include "units/unit_data.c"

static int Nx_rank;
static int nspec;
static double dt;

void initialize_fluid_electron(int Nx_rank_in, int nspec_in, double dt_in) {
  Nx_rank = Nx_rank_in;
  nspec = nspec_in;
  dt = dt_in;
}

//Note: one ion species only at hte moment
void fluid_electron_coll_only(double **n, double *Te, double **Ti) {

  /*
  // simply solving 
  // 
  // \frac32 n_e d_t(T_e) = \frac{1}{4 \eps} \nu_{ie} (T_i - T_e)
  //
  // will come back later to add advection, ohm heating, etc
  */

  //Simple first order

  double nu_ei;
  double nu_ie;
  double myn, myTe, myTi;
  double delta;
  double Z = 1.0;
  
  for(int i = 0; i < Nx_rank; ++i) {
    myn = n[i][0];
    myTe = Te[i];
    myTi = Ti[i][0];
        getColl(&myn, &myTi, myTe, &Z, &nu_ei, &nu_ie, -1, 0);
    delta = dt*0.5*nu_ie * (myTi - myTe);
    //printf("n: %g Te: %g Ti: %g nu_ie: %g dt: %g change: %g\n", myn, myTe, myTi, nu_ie, dt, delta);
    Te[i] = myTe + delta;
  }
    
}

//calculate eE = - d_x(n_e T_e) / n_e
//produces eE in units of erg / cm
//note - the ephi produced by the poiss routines are in erg
void fluid_electron_get_eE(int Nx, double *n_e, double dx, double *eE, double *T_e) {
  //calculate n_e

  //eE[0] = - (-3.0 * (n_e[0] * T_e[0]) + 4.0 * (n_e[1] * T_e[1]) - n_e[2] * T_e[2]) / (2.0 * dx) / n_e[0] / ERG_TO_EV_CGS;
  //eE[0] = - ((n_e[1] * T_e[1]) - 0.0) / (2.0 * dx) / n_e[0] / ERG_TO_EV_CGS;
  eE[0] = - ((n_e[1] * T_e[1]) - n_e[Nx-1] * T_e[Nx-1]) / (2.0 * dx) / n_e[0] / ERG_TO_EV_CGS;
  //eE[0] = 0.0;
  //printf("%e ", eE[0]  * ERG_TO_EV_CGS);
  for(int l=1; l<Nx; ++l) {
    eE[l] = - (n_e[l + 1] * T_e[l + 1] - n_e[l -1 ] * T_e[l - 1]) / (2.0 * dx) / n_e[l] / ERG_TO_EV_CGS;
    //printf("%e ", eE[l]  * ERG_TO_EV_CGS);
  }
  //eE[Nx-1] = - (3.0*(n_e[Nx - 1] * T_e[Nx - 1]) - 4.0 * (n_e[Nx - 2] * T_e[Nx - 2]) + n_e[Nx - 3] * T_e[Nx - 3]) / (2.0 * dx) / n_e[Nx - 1] / ERG_TO_EV_CGS;
  //eE[Nx-1] = - (0.0 - n_e[Nx - 2] * T_e[Nx - 2])  / (2.0 * dx) / n_e[Nx - 1] / ERG_TO_EV_CGS;
  eE[Nx-1] = - (n_e[0] * T_e[0] - n_e[Nx - 2] * T_e[Nx - 2])  / (2.0 * dx) / n_e[Nx - 1] / ERG_TO_EV_CGS;
  //eE[Nx-1] = 0.0;
  //printf("%e ", eE[Nx-1] * ERG_TO_EV_CGS);

  //eE[0] = 2*eE[1] - eE[2];
  //eE[Nx-1] = 2*eE[Nx-2] - eE[Nx-3];
  
  //printf("\n");

  //for(int l = 0; l < Nx; ++l)
  //  printf("l: %d n: %e T: %e E: %e\n", l, n_e[l], T_e[l], eE[l]*ERG_TO_EV_CGS);
  
  //fflush(stdout);
}


/* Calculates d T_e / dt + (2/3) v_e \cdot eE = 0  */
/* Factored out the n_e as in this operator split form, it should remain constant */
void fluid_electron_ohmic_term(double **n, double ***v, double *T_e, double *eE, double **Z) {
  
  double ne_approx = 0.0;
  double ve = 0.0;

  int debug = 0;
  
  double const two_third = 2.0/3.0;
  //Calculate v_e
  for(int l = 0; l < Nx_rank; ++l) {
    ne_approx = 0.0;
    ve = 0.0;

    if(debug) printf("l %d ", l);
    
    for(int s = 0; s < nspec; ++s) {
      if(debug) printf("s: %d Z: %e n: %e v: %e ", s, Z[l][s], n[l][s], v[l][s][0]);
      ne_approx += Z[l][s] * n[l][s];
      ve += Z[l][s] * n[l][s] * v[l][s][0];
    }
    ve = ve / ne_approx;

    if(debug) printf("n_e: %e ve: %e Te: %e ", ne_approx, ve, T_e[l]);

    //units  eV      s                erg / cm   cm/s  eV / erg
    T_e[l] = T_e[l] - dt * two_third * eE[l + 1] * ve * ERG_TO_EV_CGS;
    if(debug) printf("Te_new: %e \n", l, T_e[l]);
  }

  if(debug) {
    printf("\n");
    fflush(stdout);
  }
  
}
