#include "units/unit_data.c"
#include <gsl/gsl_linalg.h>
#include <math.h>
#include <stdlib.h>

gsl_matrix *computeAlpha(double *rho, double **nu, int nspec) {

  int i, j;

  gsl_matrix *alpha = gsl_matrix_calloc(nspec, nspec);

  for (i = 0; i < nspec; i++) {
    for (j = 0; j < nspec; j++) {
      if ((nu[i][j] == 0.0) && (nu[j][i] == 0.0))
        gsl_matrix_set(alpha, i, j, 0.0);
      else
        gsl_matrix_set(alpha, i, j,
                       rho[i] * nu[i][j] /
                           (rho[i] * nu[i][j] + rho[j] * nu[j][i]));
    }
  }

  return alpha;
}

gsl_matrix *computeBeta(double *n, double **nu, int nspec) {

  int i, j;

  gsl_matrix *beta = gsl_matrix_calloc(nspec, nspec);

  for (i = 0; i < nspec; i++) {
    for (j = 0; j < nspec; j++) {
      if ((nu[i][j] == 0.0) && (nu[j][i] == 0.0))
        gsl_matrix_set(beta, i, j, 0.0);
      else
        gsl_matrix_set(beta, i, j,
                       n[i] * nu[i][j] / (n[i] * nu[i][j] + n[j] * nu[j][i]));
    }
  }

  return beta;
}

gsl_matrix *computeGamma(double *m, gsl_matrix *beta, double *s, double **smix,
                         int nspec) {

  int i, j;
  double value;

  gsl_matrix *gamma = gsl_matrix_calloc(nspec, nspec);

  for (i = 0; i < nspec; i++) {
    gsl_matrix_set(gamma, i, i, 0.0);
    for (j = 0; j < i; j++) {
      value = (1.0 / 3.0) *
              (m[i] * gsl_matrix_get(beta, i, j) * (s[i] - smix[i][j]) +
               m[j] * gsl_matrix_get(beta, j, i) * (s[j] - smix[j][i]));
      gsl_matrix_set(gamma, i, j, value);
      gsl_matrix_set(gamma, j, i, value);
    }
  }

  return gamma;
}

void computeMixtureVelocities(double **v_new, gsl_matrix *alpha, int nspec,
                              double ***v_mix) {
  int i, j, dim;

  for (i = 0; i < nspec; i++)
    for (j = 0; j < nspec; j++)
      for (dim = 0; dim < 3; dim++)
        v_mix[i][j][dim] = gsl_matrix_get(alpha, i, j) * v_new[i][dim] +
                           gsl_matrix_get(alpha, j, i) * v_new[j][dim];
}

void computeMixtureTemperatures(double *T_new, gsl_matrix *beta,
                                gsl_matrix *gamma, int nspec, double **T_mix) {
  int i, j;

  for (i = 0; i < nspec; i++)
    for (j = 0; j < nspec; j++)
      T_mix[i][j] = gsl_matrix_get(beta, i, j) * T_new[i] +
                    gsl_matrix_get(beta, j, i) * T_new[j] +
                    gsl_matrix_get(gamma, i, j);
}

void implicitVelocityUpdate(double **v_old, double *rho, double **nu,
                            gsl_matrix *alpha, double dt, int nspec,
                            int conserveFlag, double **vnew) {

  int i, j;
  double rowsum, value;

  gsl_matrix *A = gsl_matrix_calloc(nspec, nspec);
  gsl_matrix *D = gsl_matrix_calloc(nspec, nspec);
  gsl_matrix *M = gsl_matrix_calloc(nspec, nspec);

  gsl_vector *vnew_vec = gsl_vector_calloc(nspec);

  gsl_vector *b = gsl_vector_calloc(nspec);

  gsl_matrix_set_identity(M);

  // Construct A and D
  for (i = 0; i < nspec; i++) {
    rowsum = 0;
    for (j = 0; j < nspec; j++) {
      if (conserveFlag == 1) { // electron background case
        value = nu[i][j];
      } else {
        value = nu[i][j] * gsl_matrix_get(alpha, j, i);
      }
      gsl_matrix_set(A, i, j, value);
      rowsum += value;
    }
    gsl_matrix_set(D, i, i, rowsum);
  }

  // Compute M = I + dt*(D - A)
  gsl_matrix_sub(D, A);
  gsl_matrix_scale(D, dt);
  gsl_matrix_add(M, D);

  int status;
  gsl_permutation *P = gsl_permutation_alloc(nspec);
  gsl_linalg_LU_decomp(M, P, &status);

  // Solve M vnew = v_old for each velocity component

  int dim;
  for (dim = 0; dim < 3; dim++) {

    for (i = 0; i < nspec; i++)
      gsl_vector_set(b, i, v_old[i][dim]);

    gsl_linalg_LU_solve(M, P, b, vnew_vec);

    for (i = 0; i < nspec; i++)
      vnew[i][dim] = gsl_vector_get(vnew_vec, i);
  }

  // Clean up
  gsl_permutation_free(P);
  gsl_matrix_free(A);
  gsl_matrix_free(D);
  gsl_matrix_free(M);
  gsl_vector_free(vnew_vec);
  gsl_vector_free(b);
}

void implicitTemperatureUpdate(double *Told, double *m, double *n,
                               double *v2old, double *v2new, double **v2mix,
                               double **nu, gsl_matrix *beta, double dt,
                               int nspec, int conserveFlag, double *Tnew) {

  int i, j;
  double valueB, valueG, rowvalB, rowvalG;

  gsl_matrix *B = gsl_matrix_calloc(nspec, nspec);
  gsl_matrix *G = gsl_matrix_calloc(nspec, nspec);
  gsl_matrix *F = gsl_matrix_calloc(nspec, nspec);
  gsl_matrix *M = gsl_matrix_calloc(nspec, nspec);

  gsl_vector *S = gsl_vector_calloc(nspec);
  gsl_vector *Tnew_vec = gsl_vector_calloc(nspec);

  for (i = 0; i < nspec; i++) {
    rowvalB = 0.0;
    rowvalG = 0.0;
    for (j = 0; j < nspec; j++) {
      if (conserveFlag == 1) {
        valueB = nu[i][j];
      } else {
        valueB = nu[i][j] * gsl_matrix_get(beta, j, i);
      }
      valueG = valueB * (-m[i] * (v2new[i] - v2mix[i][j]) +
                         m[j] * (v2new[j] - v2mix[i][j]));

      rowvalB += valueB;
      rowvalG += valueG;

      gsl_matrix_set(B, i, j, valueB); // 1 / s
      gsl_matrix_set(G, i, j, valueG); // erg / s
    }
    gsl_vector_set(S, i,
                   (m[i] * (v2old[i] - v2new[i]) + dt * rowvalG / 3.0) *
                           ERG_TO_EV_CGS +
                       Told[i]);
    //                          erg                      erg             eV/erg
    //                          eV
    gsl_matrix_set(F, i, i, rowvalB); // 1 / s
  }

  // Form M = I + dt(F-B)
  gsl_matrix_set_identity(M);

  gsl_matrix_sub(F, B);
  gsl_matrix_scale(F, dt);
  gsl_matrix_add(M, F);

  int status;
  gsl_permutation *P = gsl_permutation_alloc(nspec);
  gsl_linalg_LU_decomp(M, P, &status);

  // Solve MT = S
  gsl_linalg_LU_solve(M, P, S, Tnew_vec);

  for (i = 0; i < nspec; i++)
    Tnew[i] = gsl_vector_get(Tnew_vec, i);

  // Clean up
  gsl_matrix_free(B);
  gsl_matrix_free(G);
  gsl_matrix_free(F);
  gsl_matrix_free(M);
  gsl_vector_free(S);
  gsl_vector_free(Tnew_vec);
  gsl_permutation_free(P);
}

void implicitGetVelocitiesTemperaturesLinear(double *n, double **v, double *T,
                                             double **nu, double *m, double dt,
                                             int nspec, int conserveFlag,
                                             double **vnew, double ***vmix_new,
                                             double *Tnew, double **Tmix_new) {
  int i, j, dim;

  double rho[nspec];

  double **v2mix = malloc(nspec * sizeof(double *));
  double v2[nspec];
  double v2new[nspec];
  for (i = 0; i < nspec; i++) {
    v2mix[i] = malloc(nspec * sizeof(double));
  }

  gsl_matrix *alpha;
  gsl_matrix *beta;
  gsl_matrix *gamma;

  for (i = 0; i < nspec; i++)
    rho[i] = m[i] * n[i];

  alpha = computeAlpha(rho, nu, nspec);
  beta = computeBeta(n, nu, nspec);

  implicitVelocityUpdate(v, rho, nu, alpha, dt, nspec, conserveFlag, vnew);

  // Compute mixture velocities
  computeMixtureVelocities(vnew, alpha, nspec, vmix_new);

  // Define the v2 variables
  for (i = 0; i < nspec; i++) {
    v2[i] = 0.0;
    v2new[i] = 0.0;
    for (dim = 0; dim < 3; dim++) {
      v2[i] += v[i][dim] * v[i][dim];
      v2new[i] += vnew[i][dim] * vnew[i][dim];
    }

    for (j = 0; j < nspec; j++) {
      v2mix[i][j] = 0.0;
      for (dim = 0; dim < 3; dim++)
        v2mix[i][j] += vmix_new[i][j][dim];
    }
  }

  // Solve for temperature
  beta = computeBeta(n, nu, nspec);
  gamma = computeGamma(m, beta, v2, v2mix, nspec);

  implicitTemperatureUpdate(T, m, n, v2, v2new, v2mix, nu, beta, dt, nspec,
                            conserveFlag, Tnew);

  computeMixtureTemperatures(Tnew, beta, gamma, nspec, Tmix_new);

  // Clean up
  gsl_matrix_free(alpha);
  gsl_matrix_free(beta);
  gsl_matrix_free(gamma);

  for (i = 0; i < nspec; i++) {
    free(v2mix[i]);
  }
  free(v2mix);
}

void implicitGetVelocitiesTemperaturesNonlinear(
    double *n, double **v, double *T, double *m, double dt, int nspec,
    int conserveFlag, double **vnew, double ***vmix_new, double *Tnew,
    double **Tmix_new) {
  int i, j, dim;

  double rho[nspec];
  double Told[nspec];
  double **vold = malloc(nspec * sizeof(double));

  double **nu = malloc(nspec * sizeof(double *));
  double **v2mix = malloc(nspec * sizeof(double *));
  double v2[nspec];
  double v2new[nspec];
  for (i = 0; i < nspec; i++) {
    vold[i] = malloc(3 * sizeof(double));
    v2mix[i] = malloc(nspec * sizeof(double));
    nu[i] = malloc(nspec * sizeof(double));
  }

  gsl_matrix *alpha;
  gsl_matrix *beta;
  gsl_matrix *gamma;

  double relErr = 1.0;
  double absErr = 1.0;
  int iter = 0;
  double Tnorm;

  double relTol = 1e-8;
  double absTol = 1e-8;
  int maxIter = 100;

  // some setup stuff

  for (i = 0; i < nspec; i++) {
    rho[i] = m[i] * n[i];
    vold[i][0] = v[i][0];
    vold[i][1] = v[i][1];
    vold[i][2] = v[i][2];
    Told[i] = T[i];
  }

  while ((relErr > relTol) && (absErr > absTol) && (iter < maxIter)) {

    // Compute new values for nu
    // Temporary value...
    for (i = 0; i < nspec; i++) {
      for (j = 0; j < nspec; j++) {
        double mu = m[i] * m[j] / (m[i] + m[j]);
        nu[i][j] = 64 * M_PI * mu * mu * n[i] * n[j] *
                   pow(m[i] * Told[j] + m[j] * Told[i], 0.5) /
                   (3.0 * pow(2.0 * M_PI * m[i] * m[j], 1.5));
      }
    }

    alpha = computeAlpha(rho, nu, nspec);
    beta = computeBeta(n, nu, nspec);

    implicitVelocityUpdate(v, rho, nu, alpha, dt, nspec, conserveFlag, vnew);

    // Compute mixture velocities
    computeMixtureVelocities(vnew, alpha, nspec, vmix_new);

    // Define the v2 variables
    for (i = 0; i < nspec; i++) {
      v2[i] = 0.0;
      v2new[i] = 0.0;
      for (dim = 0; dim < 3; dim++) {
        v2[i] += vold[i][dim] * vold[i][dim];
        v2new[i] += vnew[i][dim] * vnew[i][dim];
      }

      for (j = 0; j < nspec; j++) {
        v2mix[i][j] = 0.0;
        for (dim = 0; dim < 3; dim++)
          v2mix[i][j] += vmix_new[i][j][dim];
      }
    }

    // Solve for temperature
    beta = computeBeta(n, nu, nspec);
    gamma = computeGamma(m, beta, v2, v2mix, nspec);

    implicitTemperatureUpdate(T, m, n, v2, v2new, v2mix, nu, beta, dt, nspec,
                              conserveFlag, Tnew);

    computeMixtureTemperatures(Tnew, beta, gamma, nspec, Tmix_new);

    // Test convergence of temperatures
    absErr = 0.0;
    Tnorm = 0.0;
    for (i = 0; i < nspec; i++) {
      absErr += (Tnew[i] - Told[i]) * (Tnew[i] - Told[i]);
      Tnorm += Told[i] * Told[i];
      vold[i][0] = vnew[i][0];
      vold[i][1] = vnew[i][1];
      vold[i][2] = vnew[i][2];
      Told[i] = Tnew[i];
    }
    relErr = absErr / Tnorm;
    iter += 1;
  }

  // printf("Iters %d, absErr %g, relErr %g\n\n\n", iter, absErr, relErr);

  // Clean up
  gsl_matrix_free(alpha);
  gsl_matrix_free(beta);
  gsl_matrix_free(gamma);

  for (i = 0; i < nspec; i++) {
    free(vold[i]);
    free(v2mix[i]);
    free(nu[i]);
  }
  free(vold);
  free(v2mix);
  free(nu);
}

void implicitUpdateDistributionsLinear(double **f, double *n, double **v,
                                       double *T, double **nu, double *m,
                                       double dt, int nspec, int conserveFlag,
                                       int Nv, double **fnew) {}

void implicitUpdateDistributionsNonlinear(double **f, double *n, double **v,
                                          double *T, double **nu, double *m,
                                          double dt, int nspec,
                                          int conserveFlag, int Nv,
                                          double **fnew) {}
