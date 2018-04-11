#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "initialize_sol.h"
#include "units/unit_data.c"

void initialize_sol_inhom(double ***f, int nint, double *int_loc, double *ndens_in, double *v_in, double *T_in, int Nx, double *x, int nspec, int Nv, double **c, double *m, double **n_oned, double ***v_oned, double **T_oned) {

  //accessing multid arrays
  //f - Nx nspec Nv*Nv*Nv
  //int_loc - nint
  //rho/v/T_in - nint nspec
  //x - nv
  //m - nspec
  //c - nspec nv

  int i,j,k,l,s;
  int int_id;
  int vIndex,inputIndex;

  for(l=0;l<Nx;l++) {

    //Figure out the interval
    //If no interval is found it defaults to zero
    int_id = 0; 
    for(j=0;j<nint;j++)
      if(x[l] > int_loc[j])
	int_id = j;
    
    for(s=0;s<nspec;s++) {
      inputIndex = int_id*nspec + s;

      //set moment data
      if(ndens_in[inputIndex] != 0) {
	n_oned[l][s] = ndens_in[inputIndex];
	v_oned[l][s][0] = v_in[inputIndex];
	v_oned[l][s][1] = 0.0;
	v_oned[l][s][2] = 0.0;
	T_oned[l][s] = T_in[inputIndex];	     
      }
      else {
	n_oned[l][s] = 0.0;
	v_oned[l][s][0] = 0.0;
	v_oned[l][s][1] = 0.0;
	v_oned[l][s][2] = 0.0;
	T_oned[l][s] = 0.0;
      }

      //set distribution data
      for(i=0;i<Nv;i++) 
	for(j=0;j<Nv;j++)
	  for(k=0;k<Nv;k++) {
	    vIndex = k + Nv*(j + Nv*i);	    

 	    if(ndens_in[inputIndex] != 0) {
	      f[l][s][vIndex] =  ndens_in[inputIndex]*pow(m[s]/(2.0*M_PI*T_in[inputIndex]/ERG_TO_EV_CGS),1.5)*
		exp(-(0.5*m[s]/(T_in[inputIndex]/ERG_TO_EV_CGS))*
		    ( (c[s][i]-v_in[inputIndex])*(c[s][i]-v_in[inputIndex]) 
	            + (c[s][j]*c[s][j])
		    + (c[s][k]*c[s][k]) ));
	    }
	    else {
	      f[l][s][vIndex] = 0.0;
	    }
	  }	
    }
  }


  
}


