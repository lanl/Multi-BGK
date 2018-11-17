#include <stdlib.h>
#include <stdio.h>
#include "../src/implicit.h"
#include <gsl/gsl_matrix.h>

int main() {

  int nspec = 5;
  double dt = 1.0e-12;
  int i,j, dim;

  int conserveFlag = 1;

  double **nu = malloc(nspec * sizeof(double *));
  for(i=0;i<nspec; i++) {
    nu[i] = malloc(nspec * sizeof(double));
  }

  double **v     = malloc(nspec * sizeof(double *));
  double **vnew  = malloc(nspec * sizeof(double *));
  double ***vmix = malloc(nspec * sizeof(double **));
  for(i=0;i< nspec;i++) {
    v[i]     = malloc(3*sizeof(double));
    vnew[i]  = malloc(3*sizeof(double));
    vmix[i]  = malloc(nspec * sizeof(double *));
    for(j=0;j<nspec;j++)
      vmix[i][j] = malloc(3*sizeof(double));
  }

  double m[nspec];
  double n[nspec];

  double T[nspec];
  double Tnew[nspec];
  double **Tmix = malloc(nspec * sizeof(double *));
  for(i=0;i<nspec;i++)
    Tmix[i] = malloc(nspec*sizeof(double));
  

  //Put in some dummy values
  for(i=0;i<nspec;i++) {
    m[i] = 1.0;
    n[i] = 2.0;

    v[i][0] =  0.0;
    v[i][1] =  0.0;
    v[i][2] =  0.0;

    T[i] = 1.0 + i;

    for(j=0;j<nspec;j++)
      nu[i][j] = 1.0e10;
  }

  double amu = 1.6605e-24;
  m[0] = 12*amu;
  m[1] = amu;
  m[2] = 2.0*amu;
  m[3] = 3.0*amu;
  m[4] = amu / 1876;

  n[0] = 3.6e21;
  n[1] = 3.6e21;
  n[2] = 1.0e10;
  n[3] = 1.0e10;

  T[0] = 100.0;
  T[1] = 100.0;
  T[2] = 100.0;
  T[3] = 100.0;
  T[4] = 100.0;

  for(j=0;j<nspec;j++)
    nu[4][j] = 0.0;

  implicitGetVelocitiesTemperaturesLinear(n, v, T, nu, m, dt, nspec, conserveFlag, vnew, vmix, Tnew, Tmix);

  printf("New velocity\n");
  for (i = 0; i < nspec; i++)
    printf("%5e ", vnew[i][0]);
  printf("\n");

  /*
  printf("New mix velocity\n");
  printf("----------------\n");
  for (i = 0; i < nspec; i++)
    for (j = 0; j < nspec; j++)
      printf("%5e ", vmix[i][j][0]);
  printf("\n");
  printf("----------------\n\n");
  */

  printf("New Temperature\n");
  for (i = 0; i < nspec; i++)
    printf("%5e ", Tnew[i]);
  printf("\n");

  /*
  printf("New mix Temperatures\n");
  printf("----------------\n");
  for (i = 0; i < nspec; i++)
    for (j = 0; j < nspec; j++)
      printf("%5e ", Tmix[i][j]);
  printf("\n");
  printf("----------------\n\n");
  */

  /*
  printf("\n\n\n\n nonlinear solve\n\n\n------------------\n");

  implicitGetVelocitiesTemperaturesNonlinear(n, v, T, m, dt, nspec, conserveFlag, vnew, vmix, Tnew, Tmix);

  printf("New velocity\n");
  printf("----------------\n");
  for (i = 0; i < nspec; i++)
    printf("%5e ", vnew[i][0]);
  printf("\n");
  printf("----------------\n\n");

  printf("New mix velocity\n");
  printf("----------------\n");
  for (i = 0; i < nspec; i++)
    for (j = 0; j < nspec; j++)
      printf("%5e ", vmix[i][j][0]);
  printf("\n");
  printf("----------------\n\n");


  printf("New Temperature\n");
  printf("----------------\n");
  for (i = 0; i < nspec; i++)
    printf("%5e ", Tnew[i]);
  printf("\n");
  printf("----------------\n\n");

  printf("New mix Temperatures\n");
  printf("----------------\n");
  for (i = 0; i < nspec; i++)
    for (j = 0; j < nspec; j++)
      printf("%5e ", Tmix[i][j]);
  printf("\n");
  printf("----------------\n\n");
  */

  return 0;
}
