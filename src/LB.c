#include "LB.h"
#include "momentRoutines.h"
#include <stdio.h>
#include <stdlib.h>
#include "BGK.h"
#include <math.h>
#include "units/unit_data.c"

/* *********************** 
Calculate the MS LB as

Q[fi, fj] = nu_ij d_v ( (c - v_ij) f_i + T_ij / m_i * d_v f_i)

*************************/


#define EPS_COLL 1e-6

static int Nv;
static int nspec;
static double *m;
static double **c;
static double **Q;
static double *n;
static double *rho;
static double **v;
static double *T;
static double **nu;

//effectively a ctor
void initialize_LB(int ns, int numV, double *mass, double **vels) {

  Nv = numV;
  nspec = ns;
  m = mass;
  c = vels;

  //Initialize local arrays
  n = malloc(nspec * sizeof(double));
  rho = malloc(nspec * sizeof(double));
  T = malloc(nspec * sizeof(double));
  
  nu = malloc(nspec * sizeof(double *));
  v = malloc(nspec * sizeof(double *));
  for (int i = 0; i < nspec; ++i) {
    v[i] = malloc(3 * sizeof(double));
    nu[i] = malloc(ns * sizeof(double));
  }
  
  //printf("Lb params: Nv: %d ns: %d mass: %g velo[Nv/2]: %g\n", Nv, nspec, m[0], c[0][Nv/2]);
  fflush(stdout);
}

void LB_ex(double **f, double **f_out, double *Z, double dt, double Te) {

  double ntot = 0;
  double rhotot = 0;

  // setup: initialize output array and accum array
  for (int i = 0; i < nspec; i++)
    for (int k = 0; k < Nv * Nv * Nv; k++) {
      f_out[i][k] = 0.0;
    }

  //First step - get the moments
  for(int i = 0; i < nspec; ++i) {
    n[i] = getDensity(f[i],i);
    rho[i] = m[i] * n[i];
    ntot += n[i];
    rhotot += rho[i];
    getBulkVel(f[i], v[i], n[i], i);
    T[i] = getTemp(m[i], n[i], v[i], f[i], i);
  }

  // check for blowup
  if (isnan(n[0])) {
    printf("ERROR ERROR  NaN detected in "
           "LB.c\n");
    for (int i = 0; i < nspec; i++) {
      printf("%d n: %g v: %g T: %g Z: %g Te: %g\n", i, n[i], v[i][0], T[i],
             Z[i], Te);
    }
    exit(1);
  }


  // get collision rates

  double nu11, nu12, nu21;
  //do ij and ji at the same time
  for(int i = 0; i < nspec; ++i) {
    for(int j = 0; j <  nspec; ++j) {
      nu11 = 1.0;
      nu12 = 1.0;
      nu21 = 1.0;
      
      //self collision case
      double friction_accumulatrix;
      double diffusion_accumulatrix;
      double L_z, U_z, L_y, U_y, L_x, U_x;
      double c_L_z, c_U_z, c_L_y, c_U_y, c_L_x, c_U_x;
      
      if( i == j ) {
	if ( n[j] > EPS_COLL ) {

	  getColl(n, T, Te, Z, &nu11, &nu11, i, j);
	  
	  double dv = c[i][1] - c[i][0];
	  for(int index_x = 0; index_x < Nv; index_x++) {
	    
	    for(int index_y = 0; index_y < Nv; index_y++) {
	      
	      for(int index_z = 0; index_z < Nv; index_z++) {
		
		L_x = (index_x == 0) ? 0.0 : f[i][index_z + Nv * (index_y + Nv * (index_x - 1))];
		L_y = (index_y == 0) ? 0.0 : f[i][index_z + Nv * (index_y - 1 + Nv * (index_x))];
		L_z = (index_z == 0) ? 0.0 : f[i][index_z - 1 + Nv * (index_y + Nv * (index_x))];

		c_L_x = (index_x == 0) ? c[i][0] - dv : c[i][index_x - 1];
		c_L_y = (index_y == 0) ? c[i][0] - dv : c[i][index_y - 1];
		c_L_z = (index_z == 0) ? c[i][0] - dv : c[i][index_z - 1];

		U_x = (index_x == (Nv - 1)) ? 0.0 : f[i][index_z + Nv * (index_y + Nv * (index_x + 1))];
		U_y = (index_y == (Nv - 1)) ? 0.0 : f[i][index_z + Nv * (index_y + 1 + Nv * (index_x))];
		U_z = (index_z == (Nv - 1)) ? 0.0 : f[i][index_z + 1 + Nv * (index_y + Nv * (index_x))];

		c_U_x = (index_x == (Nv - 1)) ? c[i][Nv - 1] + dv : c[i][index_x + 1];
		c_U_y = (index_y == (Nv - 1)) ? c[i][Nv - 1] + dv : c[i][index_y + 1];
		c_U_z = (index_z == (Nv - 1)) ? c[i][Nv - 1] + dv : c[i][index_z + 1];


		int k = index_z + Nv * (index_y + Nv * (index_x));
		friction_accumulatrix = 0.0;
		diffusion_accumulatrix = 0.0;

		//Do I need to upwind this?
		/*
		if ((c[i][index_x] - v[i][0]) > 0) {
		  friction_accumulatrix += ((c[i][index_x] - v[i][0]) * f[i][k] - (c_L_x - v[i][0]) * L_x) / dv;
		}
		else {
		  friction_accumulatrix += ((c_U_x - v[i][0]) * U_x - (c[i][index_x] - v[i][0]) * f[i][k]) / dv;
		}
		
		if ((c[i][index_y] - v[i][1]) > 0) {
		  friction_accumulatrix += ((c[i][index_y] - v[i][1]) * f[i][k] - (c_L_y - v[i][1]) * L_y) / dv;
		}
		else {
		  friction_accumulatrix += ((c_U_y - v[i][1]) * U_x - (c[i][index_y] - v[i][1]) * f[i][k]) / dv;
		}
		
		if ((c[i][index_z] - v[i][2]) > 0) {
		  friction_accumulatrix += ((c[i][index_z] - v[i][2]) * f[i][k] - (c_L_z - v[i][2]) * L_z) / dv;
		}
		else {
		  friction_accumulatrix += ((c_U_z - v[i][2]) * U_x - (c[i][index_z] - v[i][2]) * f[i][k]) / dv;
		}
*/
		
		friction_accumulatrix += ((c_U_x - v[i][0]) * U_x - (c_L_x - v[i][0]) * L_x) / dv * 0.5;
		friction_accumulatrix += ((c_U_y - v[i][1]) * U_y - (c_L_y - v[i][1]) * L_y) / dv * 0.5;
		friction_accumulatrix += ((c_U_z - v[i][2]) * U_z - (c_L_z - v[i][2]) * L_z) / dv * 0.5;
		
		/*
		friction_accumulatrix += ((c[i][index_x + 1] - v[i][0]) * f[i][index_z + Nv * (index_y + Nv * (index_x))] -
					  (c[i][index_x - 1] - v[i][0]) * f[i][index_z + Nv * (index_y + Nv * (index_x))]) / dv;
		friction_accumulatrix += ((c[i][index_y + 1] - v[i][1]) * f[i][index_z + Nv * (index_y + 1 + Nv * (index_x))] -
					  (c[i][index_y - 1] - v[i][1]) * f[i][index_z + Nv * (index_y - 1 + Nv * (index_x))]) / dv;
		friction_accumulatrix += ((c[i][index_z + 1] - v[i][2]) * f[i][index_z + 1 + Nv * (index_y + Nv * (index_x))] -
					  (c[i][index_z - 1] -v[i][2]) * f[i][index_z - 1+ Nv * (index_y + Nv * (index_x))]) / dv;
		*/

		diffusion_accumulatrix += (U_x + L_x - 2*f[i][k]);
		diffusion_accumulatrix += (U_y + L_y - 2*f[i][k]);
		diffusion_accumulatrix += (U_z + L_z - 2*f[i][k]);

		/*
		diffusion_accumulatrix += (f[i][index_z + 1 + Nv * (index_y + Nv * (index_x))] + f[i][index_z - 1 + Nv * (index_y + Nv * (index_x))] - 2.0 * f[i][k]);
		diffusion_accumulatrix += (f[i][index_z + Nv * (index_y + 1 + Nv * (index_x))] + f[i][index_z + Nv * (index_y - 1 + Nv * (index_x))] - 2.0 * f[i][k]);
		diffusion_accumulatrix += (f[i][index_z + Nv * (index_y + Nv * (index_x + 1))] + f[i][index_z + Nv * (index_y + Nv * (index_x - 1))] - 2.0 * f[i][k]);
		*/
		
		diffusion_accumulatrix *= (T[i] / ERG_TO_EV_CGS) / m[i] / dv / dv;
		
		f_out[i][k] = nu11 * (friction_accumulatrix + diffusion_accumulatrix);
	      }
	    }
	  }
	}
      }
      else {
	printf("LB not implemented for multispecies yet\n");
	exit(1);	  
      }
      

      
    }    
  }
    

  
}
