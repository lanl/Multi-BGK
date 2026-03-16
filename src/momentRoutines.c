#include "momentRoutines.h"
#include "units/unit_data.c"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int N;
static int nspec;
static double **c;
static double **wts;

// static double KB;
// static const double KB_true = 1.380658e-23; //Boltzmann constant

/*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/

void initialize_moments(int nodes, int ns, double **vel, double **weights) {

  N = nodes;
  nspec = ns;
  c = vel;
  wts = weights;
}

double getDensity(double *in, int sp) {
  double result = 0.0;
  int i, j, k;

  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      for (k = 0; k < N; k++)
        result +=
            wts[sp][i] * wts[sp][j] * wts[sp][k] * in[k + N * (j + N * i)];

  return result;
}

/*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/

void getBulkVel(double *in, double *out, double n, int sp) {
  double temp1, temp2, temp3;
  int i, j, k;

  out[0] = 0.0;
  out[1] = 0.0;
  out[2] = 0.0;

  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      for (k = 0; k < N; k++) {
        temp1 = c[sp][i] * wts[sp][i] * wts[sp][j] * wts[sp][k];
        temp2 = c[sp][j] * wts[sp][i] * wts[sp][j] * wts[sp][k];
        temp3 = c[sp][k] * wts[sp][i] * wts[sp][j] * wts[sp][k];

        out[0] += temp1 * in[k + N * (j + N * i)];
        out[1] += temp2 * in[k + N * (j + N * i)];
        out[2] += temp3 * in[k + N * (j + N * i)];
      }

  if (n != 0) {
    out[0] = out[0] / n;
    out[1] = out[1] / n;
    out[2] = out[2] / n;
  }
}

/*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/

double getTemp(double m, double n, double *u, double *in, int sp) {
  double result = 0.0;

  int i, j, k;
  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      for (k = 0; k < N; k++) {
        result += wts[sp][i] * wts[sp][j] * wts[sp][k] *
                  in[k + N * (j + N * i)] *
                  ((c[sp][i] - u[0]) * (c[sp][i] - u[0]) +
                   (c[sp][j] - u[1]) * (c[sp][j] - u[1]) +
                   (c[sp][k] - u[2]) * (c[sp][k] - u[2]));
      }

  // Int(|v-u|^2 f) = 3* nT/m

  if (n != 0)
    return result * m * ERG_TO_EV_CGS / (3 * n);
  else
    return 0.0;
}

//gets directional temperature
double getDirTemp(double m, double n, double *u, double *in, int sp, int dir) {
  double result = 0.0;
  
  if ((dir < 0) || (dir > 2)) {
    printf("Error calculating directional temperature, please set dir to 0 (x), 1 (y), or 2 (z)\n");
    fflush(stdout);
    exit(1);
  }
		    		   
  int i, j, k, mydir;
  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      for (k = 0; k < N; k++) {
	if(dir == 0)
	  mydir = i;
	else if (dir == 1)
	  mydir = j;
	else
	  mydir = k;
	
        result += wts[sp][i] * wts[sp][j] * wts[sp][k] *
	          in[k + N * (j + N * i)] *
      	          ((c[sp][mydir] - u[dir]) * (c[sp][mydir] - u[dir]));
	  }
  
  // Int(|v-u|^2 f) = 3* nT/m

  if (n != 0)
    return result * m * ERG_TO_EV_CGS / n;
  else
    return 0.0;
}

//computes n Theta, the pressure/stress tensor
void getStressTensor(double m, double n, double *u, double *in, double *out, int sp) {

  for(int i = 0; i < 3; ++i) {
    for(int j = 0; j < 3; ++j) {
      double accum = 0;
      int index1, index2;
      for(int vx = 0; vx < N; ++vx)
	for(int vy = 0; vy < N; ++vy)
	  for(int vz = 0; vz < N; ++vz) {
	    if(i == 0)
	      index1 = vx;
	    else if (i == 1)
	      index1 = vy;
	    else
	      index1 = vz;
		
	    if(j == 0)
	      index2 = vx;
	    else if (j == 1)
	      index2 = vy;
	    else
	      index2 = vz;

	    //units are 1/cm^3 * cm^2/s^2 = 1/cm/s^2
	    accum += wts[sp][vz] * wts[sp][vy] * wts[sp][vx] * in[vz + N * (vy + N*vx)] *
	      ((c[sp][index1] - u[i]) * (c[sp][index2] - u[j]));
	  }
      out[j + 3*i] = m * accum; //units are g / cm / s^2 
    }
  }
      
  
}

double getH(double n, double *in, int sp) {
  double result = 0.0;

  int i, j, k;
  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      for (k = 0; k < N; k++) {
        result += wts[sp][i] * wts[sp][j] * wts[sp][k] *
                  in[k + N * (j + N * i)] * log(in[k + N * (j + N * i)]);
      }

  if (n != 0)
    return result;
  else
    return 0.0;
}

// Calculates the marginal f(vx) = \int f(vx, vy, vz) dvy dvz
void getMarginal(double *in, double *out, int sp) {

  int i, j, k;
  double result;

  for (i = 0; i < N; ++i) {
    result = 0.0;
    for (j = 0; j < N; ++j) {
      for (k = 0; k < N; ++k) {
        result += wts[sp][j] * wts[sp][k] * in[k + N * (j + N * i)];
      }
    }
    out[i] = result;
  }
}
