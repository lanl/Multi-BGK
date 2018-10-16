#include "momentRoutines.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "units/unit_data.c"

static int N;
static int nspec;
static double **c;
static double **wts;

//static double KB;
//static const double KB_true = 1.380658e-23; //Boltzmann constant


/*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/

void initialize_moments(int nodes, int ns, double **vel, double **weights) {

  N = nodes;
  nspec = ns;
  c = vel;
  wts = weights;
}


double getDensity(double *in, int sp)
{
  double result = 0.0;
  int i, j, k;
  
  for(i=0;i<N;i++) 
    for(j=0;j<N;j++)
      for(k=0;k<N;k++)	
	result += wts[sp][i]*wts[sp][j]*wts[sp][k]*in[k + N*(j + N*i)];		  

  return result;
}

/*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/

void getBulkVel(double *in, double *out, double n, int sp)
{
  double temp1, temp2, temp3;
  int i, j, k;
  
  out[0] = 0.0;
  out[1] = 0.0;
  out[2] = 0.0;
  
  for(i=0;i<N;i++)
    for(j=0;j<N;j++)
      for(k=0;k<N;k++)
	{
	  temp1 = c[sp][i]*wts[sp][i]*wts[sp][j]*wts[sp][k];
	  temp2 = c[sp][j]*wts[sp][i]*wts[sp][j]*wts[sp][k];
	  temp3 = c[sp][k]*wts[sp][i]*wts[sp][j]*wts[sp][k];
	  
	  out[0] += temp1*in[k + N*(j + N*i)];
	  out[1] += temp2*in[k + N*(j + N*i)];
	  out[2] += temp3*in[k + N*(j + N*i)];
	}
  
  if(n != 0) {
    out[0] = out[0]/n;
    out[1] = out[1]/n;
    out[2] = out[2]/n;
  }
}

/*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/

double getTemp(double m, double n, double *u, double *in, int sp)
{
  double result = 0.0;

  int i, j, k;
  for(i=0;i<N;i++)
    for(j=0;j<N;j++)
      for(k=0;k<N;k++) {
	result += wts[sp][i]*wts[sp][j]*wts[sp][k]*in[k + N*(j + N*i)]*((c[sp][i]-u[0])*(c[sp][i]-u[0]) + (c[sp][j]-u[1])*(c[sp][j]-u[1]) + (c[sp][k]-u[2])*(c[sp][k]-u[2]));
      }

  //Int(|v-u|^2 f) = 3* nT/m    

  if(n != 0)
    return result*m*ERG_TO_EV_CGS/(3*n);
  else
    return 0.0;
  
}

double getH(double n, double *in, int sp) {
  double result = 0.0;

  int i, j, k;
  for(i=0;i<N;i++)
    for(j=0;j<N;j++)
      for(k=0;k<N;k++) {
	result += wts[sp][i]*wts[sp][j]*wts[sp][k]*in[k + N*(j + N*i)]*log(in[k + N*(j + N*i)]);
      }

  if(n != 0)
    return result;
  else
    return 0.0;
  
}
