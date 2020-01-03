#include "units/unit_data.c"
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef EPS_COLL
#define EPS_COLL 1e6
#endif

static int Nv;
static double **c;
static double **wts;

static int first_DD = 1;
static int first_DT = 1;

void initializeTNB(int Nv_in, double **c_in, double **wts_in) {
  printf("Initializing TNB\n");
  Nv = Nv_in;
  c = c_in;
  wts = wts_in;
}

// Added by Abdou
// This function calculates the rate of change of the distribution of species i
// due to its reaction with j.
// \int f_i(v_i) f_j(v_j) |v_i-v_j| cross_section dv_j

double GetTNB_dt(double mu, double *in, double *c1, int sp, int sp2) {
  double result = 0.0;
  double f2;
  double g1, g2, g3, g;
  double E_COM;
  double Nuclear_factor;
  double Nuclear_factor_num;
  double Nuclear_factor_dem;
  double cross_section;
  // TNB constants
  double a1 = 6.927e4;
  double a2 = 7.454e8;
  double a3 = 2.050e6;
  double a4 = 5.2002e4;
  double a5 = 0;
  double b1 = 6.38e1;
  double b2 = -9.95e-1;
  double b3 = 6.981e-5;
  double b4 = 1.728e-4;
  double B_G = 34.33827;

  if (mu > 2.e-24 || mu < 1.8e-24) {
    result += 0.;

    return result;
  }

  int i, j, k;
#pragma omp parallel for reduction(+ : result) private(                        \
    i, j, k, g1, g2, g3, g, E_COM, Nuclear_factor_num, Nuclear_factor_dem,     \
    Nuclear_factor, cross_section, f2)
  for (i = 0; i < Nv; i++)
    for (j = 0; j < Nv; j++)
      for (k = 0; k < Nv; k++) {

        g1 = (c1[0] - c[sp2][i]);
        g2 = (c1[1] - c[sp2][j]);
        g3 = (c1[2] - c[sp2][k]);

        g = sqrt(pow(g1, 2) + pow(g2, 2) + pow(g3, 2));
        //                double mu=m[i]*m2[j]/(m[i]+m2[j]);
        E_COM = 0.5 * mu * pow(g, 2) * ERG_TO_EV_CGS *
                1e-3; // Energy of the center-of-mass in keV

        Nuclear_factor_num =
            a1 + E_COM * (a2 + E_COM * (a3 + E_COM * (a4 + E_COM * a5)));
        Nuclear_factor_dem =
            1.0 + E_COM * (b1 + E_COM * (b2 + E_COM * (b3 + E_COM * b4)));
        Nuclear_factor = Nuclear_factor_num / Nuclear_factor_dem;

        // Calculate sigma_DT here:

        if (E_COM != 0)
          cross_section = Nuclear_factor * exp(-B_G / sqrt(E_COM)) / E_COM;
        else
          cross_section = 0;

        cross_section *= 1e-27; // convert from mbar to cm^2

        // Calculate the reactivity here:
        f2 = wts[sp2][i] * wts[sp2][j] * wts[sp2][k];

        result += g * cross_section * f2 * in[k + Nv * (j + Nv * i)];
      }

  return result;
}

//

// This function calculates the total reactivity of f1 and f2.
// \int f_i(v_i) f_j(v_j) |v_i-v_j| cross_section dv_j
double GetReactivity_dt(double mu, double *in, double *in2, int sp, int sp2) {
  double result = 0.0;
  double g1, g2, g3, g; // cm / s
  double f1, f2;        // 1 / cm^3 / (cm^3 / s^3)
  double E_COM;         // keV
  double Nuclear_factor;
  double Nuclear_factor_num;
  double Nuclear_factor_dem;
  double cross_section;

  // TNB constants
  double a1 = 6.927e4;
  double a2 = 7.454e8;
  double a3 = 2.050e6;
  double a4 = 5.2002e4;
  double a5 = 0;
  double b1 = 6.38e1;
  double b2 = -9.95e-1;
  double b3 = 6.981e-5;
  double b4 = 1.728e-4;
  double B_G = 34.33827;

  // check to see that we have DT
  if (mu > 2.e-24 || mu < 1.8e-24) {
    result += 0.;

    return result;
  }

  int i, j, k, l, m, n;

#pragma omp parallel for reduction(+ : result) private(                        \
    i, j, k, l, m, n, g1, g2, g3, g, E_COM, Nuclear_factor_num,                \
    Nuclear_factor_dem, Nuclear_factor, cross_section, f1, f2)
  for (i = 0; i < Nv; i++)
    for (j = 0; j < Nv; j++)
      for (k = 0; k < Nv; k++)
        for (l = 0; l < Nv; l++)
          for (m = 0; m < Nv; m++)
            for (n = 0; n < Nv; n++) {

              g1 = (c[sp][i] - c[sp2][l]);
              g2 = (c[sp][j] - c[sp2][m]);
              g3 = (c[sp][k] - c[sp2][n]);

              g = sqrt(pow(g1, 2) + pow(g2, 2) + pow(g3, 2));
              //                double mu=m[i]*m2[j]/(m[i]+m2[j]);

              E_COM = 0.5 * mu * pow(g, 2) * ERG_TO_EV_CGS *
                      1e-3; // Energy of the center-of-mass in keV

              Nuclear_factor_num =
                  a1 + E_COM * (a2 + E_COM * (a3 + E_COM * (a4 + E_COM * a5)));
              Nuclear_factor_dem =
                  1.0 + E_COM * (b1 + E_COM * (b2 + E_COM * (b3 + E_COM * b4)));
              Nuclear_factor = Nuclear_factor_num / Nuclear_factor_dem;

              // Calculate sigma_DT here:

              if (E_COM != 0)
                cross_section =
                    Nuclear_factor * exp(-B_G / sqrt(E_COM)) / E_COM;
              else
                cross_section = 0;

              cross_section *= 1e-27; // convert from mbar to cm^2

              // Calculate the reactivity here:

              f1 = wts[sp][i] * wts[sp][j] * wts[sp][k];
              f2 = wts[sp2][l] * wts[sp2][m] * wts[sp2][n];
              result += g * cross_section * f1 * in[k + Nv * (j + Nv * i)] *
                        f2 * in2[n + Nv * (m + Nv * l)];
            }

  return result;
}

double GetTNB_dd_He(double mu, double *in, double *c1, int sp, int sp2) {
  // TNB constants
  double a1 = 5.3701e4;
  double a2 = 3.3027e2;
  double a3 = -1.2706e-1;
  double a4 = 2.9327e-5;
  double a5 = -2.5151e-9;
  double b1 = 0.0;
  double b2 = 0.0;
  double b3 = 0.0;
  double b4 = 0.0;
  double B_G = 31.3970;

  double result = 0.0;
  double g1, g2, g3, g;
  double f2;
  double E_COM;
  double Nuclear_factor;
  double Nuclear_factor_num;
  double Nuclear_factor_dem;
  double cross_section;

  // check to see that we have DD
  if (mu > 1.7e-24 || mu < 1.6e-24) {
    result += 0.;

    return result;
  }

  int i, j, k;
#pragma omp parallel for reduction(+ : result) private(                        \
    i, j, k, g1, g2, g3, g, E_COM, Nuclear_factor_num, Nuclear_factor_dem,     \
    Nuclear_factor, cross_section, f2)

  for (i = 0; i < Nv; i++)
    for (j = 0; j < Nv; j++)
      for (k = 0; k < Nv; k++) {

        g1 = (c1[0] - c[sp2][i]);
        g2 = (c1[1] - c[sp2][j]);
        g3 = (c1[2] - c[sp2][k]);

        g = sqrt(pow(g1, 2) + pow(g2, 2) + pow(g3, 2));
        //                double mu=m[i]*m2[j]/(m[i]+m2[j]);
        E_COM = 0.5 * mu * pow(g, 2) * ERG_TO_EV_CGS *
                1e-3; // Energy of the center-of-mass in keV

        Nuclear_factor_num =
            a1 + E_COM * (a2 + E_COM * (a3 + E_COM * (a4 + E_COM * a5)));
        Nuclear_factor_dem =
            1.0 + E_COM * (b1 + E_COM * (b2 + E_COM * (b3 + E_COM * b4)));
        Nuclear_factor = Nuclear_factor_num / Nuclear_factor_dem;

        // Calculate sigma_DT here:

        if (E_COM != 0)
          cross_section = Nuclear_factor * exp(-B_G / sqrt(E_COM)) / E_COM;
        else
          cross_section = 0;

        cross_section *= 1e-27; // convert from mbar to cm^2

        // Calculate the reactivity here:
        f2 = wts[sp2][i] * wts[sp2][j] * wts[sp2][k];

        result += g * cross_section * f2 * in[k + Nv * (j + Nv * i)];
      }

  return result;
}

double GetTNB_dd_T(double mu, double *in, double *c1, int sp, int sp2) {
  // TNB constants
  double a1 = 5.5576e4;
  double a2 = 2.1054e2;
  double a3 = -3.2638e-2;
  double a4 = 1.4987e-6;
  double a5 = 1.8181e-10;
  double b1 = 0.0;
  double b2 = 0.0;
  double b3 = 0.0;
  double b4 = 0.0;
  double B_G = 31.3970;

  double result = 0.0;
  double g1, g2, g3, g;
  double f2;
  double E_COM;
  double Nuclear_factor;
  double Nuclear_factor_num;
  double Nuclear_factor_dem;
  double cross_section;

  // check to see that we have DD
  if (mu > 1.7e-24 || mu < 1.6e-24) {
    result += 0.;

    return result;
  }

  int i, j, k;
#pragma omp parallel for reduction(+ : result) private(                        \
    i, j, k, g1, g2, g3, g, E_COM, Nuclear_factor_num, Nuclear_factor_dem,     \
    Nuclear_factor, cross_section, f2)

  for (i = 0; i < Nv; i++)
    for (j = 0; j < Nv; j++)
      for (k = 0; k < Nv; k++) {

        g1 = (c1[0] - c[sp2][i]);
        g2 = (c1[1] - c[sp2][j]);
        g3 = (c1[2] - c[sp2][k]);

        g = sqrt(pow(g1, 2) + pow(g2, 2) + pow(g3, 2));
        //                double mu=m[i]*m2[j]/(m[i]+m2[j]);

        E_COM = 0.5 * mu * pow(g, 2) * ERG_TO_EV_CGS *
                1e-3; // Energy of the center-of-mass in keV

        Nuclear_factor_num =
            a1 + E_COM * (a2 + E_COM * (a3 + E_COM * (a4 + E_COM * a5)));
        Nuclear_factor_dem =
            1.0 + E_COM * (b1 + E_COM * (b2 + E_COM * (b3 + E_COM * b4)));
        Nuclear_factor = Nuclear_factor_num / Nuclear_factor_dem;

        // Calculate sigma_DT here:

        if (E_COM != 0)
          cross_section = Nuclear_factor * exp(-B_G / sqrt(E_COM)) / E_COM;
        else
          cross_section = 0;

        cross_section *= 1e-27; // convert from mbar to cm^2

        // Calculate the reactivity here:
        f2 = wts[sp2][i] * wts[sp2][j] * wts[sp2][k];

        result += g * cross_section * f2 * in[k + Nv * (j + Nv * i)];
      }

  return result;
}

double GetReactivity_dd_He(double mu, double *in, double *in2, int sp,
                           int sp2) {
  // TNB constants
  double a1 = 5.3701e4;
  double a2 = 3.3027e2;
  double a3 = -1.2706e-1;
  double a4 = 2.9327e-5;
  double a5 = -2.5151e-9;
  double b1 = 0.0;
  double b2 = 0.0;
  double b3 = 0.0;
  double b4 = 0.0;
  double B_G = 31.3970;

  double result = 0.0;
  double g1, g2, g3, g;
  double f1, f2;
  double E_COM;
  double Nuclear_factor;
  double Nuclear_factor_num;
  double Nuclear_factor_dem;
  double cross_section;

  // check to see that we have DD
  if (mu > 1.7e-24 || mu < 1.6e-24) {
    result += 0.;

    return result;
  }

  int i, j, k, l, m, n;

#pragma omp parallel for reduction(+ : result) private(                        \
    i, j, k, l, m, n, g1, g2, g3, g, E_COM, Nuclear_factor_num,                \
    Nuclear_factor_dem, Nuclear_factor, cross_section, f1, f2)
  for (i = 0; i < Nv; i++)
    for (j = 0; j < Nv; j++)
      for (k = 0; k < Nv; k++)
        for (l = 0; l < Nv; l++)
          for (m = 0; m < Nv; m++)
            for (n = 0; n < Nv; n++) {

              g1 = (c[sp][i] - c[sp2][l]);
              g2 = (c[sp][j] - c[sp2][m]);
              g3 = (c[sp][k] - c[sp2][n]);

              g = sqrt(pow(g1, 2) + pow(g2, 2) + pow(g3, 2));
              //                double mu=m[i]*m2[j]/(m[i]+m2[j]);

              E_COM = 0.5 * mu * pow(g, 2) * ERG_TO_EV_CGS *
                      1e-3; // Energy of the center-of-mass in keV

              Nuclear_factor_num =
                  a1 + E_COM * (a2 + E_COM * (a3 + E_COM * (a4 + E_COM * a5)));
              Nuclear_factor_dem =
                  1.0 + E_COM * (b1 + E_COM * (b2 + E_COM * (b3 + E_COM * b4)));
              Nuclear_factor = Nuclear_factor_num / Nuclear_factor_dem;

              // Calculate sigma_DT here:

              if (E_COM != 0)
                cross_section = Nuclear_factor * exp(-B_G / sqrt(E_COM)) /
                                E_COM; // millibarn
              else
                cross_section = 0;

              cross_section *= 1e-27; // converted from mbarn to cm^2

              // Calculate the reactivity

              f1 = wts[sp][i] * wts[sp][j] * wts[sp][k];
              f2 = wts[sp2][l] * wts[sp2][m] * wts[sp2][n];

              result += g * cross_section * f1 * in[k + Nv * (j + Nv * i)] *
                        f2 * in2[n + Nv * (m + Nv * l)];
            }

  return result;
}

double GetReactivity_dd_T(double mu, double *in, double *in2, int sp, int sp2) {
  // TNB constants
  double a1 = 5.5576e4;
  double a2 = 2.1054e2;
  double a3 = -3.2638e-2;
  double a4 = 1.4987e-6;
  double a5 = 1.8181e-10;
  double b1 = 0.0;
  double b2 = 0.0;
  double b3 = 0.0;
  double b4 = 0.0;
  double B_G = 31.3970;

  double result = 0.0;
  double g1, g2, g3, g;
  double f1, f2;
  double E_COM;
  double Nuclear_factor;
  double Nuclear_factor_num;
  double Nuclear_factor_dem;
  double cross_section;

  // check to see that we have DD
  if (mu > 1.7e-24 || mu < 1.6e-24) {
    result += 0.;

    return result;
  }

  int i, j, k, l, m, n;
#pragma omp parallel for reduction(+ : result) private(                        \
    i, j, k, l, m, n, g1, g2, g3, g, E_COM, Nuclear_factor_num,                \
    Nuclear_factor_dem, Nuclear_factor, cross_section, f1, f2)
  for (i = 0; i < Nv; i++)
    for (j = 0; j < Nv; j++)
      for (k = 0; k < Nv; k++)
        for (l = 0; l < Nv; l++)
          for (m = 0; m < Nv; m++)
            for (n = 0; n < Nv; n++) {

              g1 = (c[sp][i] - c[sp2][l]);
              g2 = (c[sp][j] - c[sp2][m]);
              g3 = (c[sp][k] - c[sp2][n]);

              g = sqrt(pow(g1, 2) + pow(g2, 2) + pow(g3, 2));
              //                double mu=m[i]*m2[j]/(m[i]+m2[j]);

              E_COM = 0.5 * mu * pow(g, 2) * ERG_TO_EV_CGS *
                      1e-3; // Energy of the center-of-mass in keV

              Nuclear_factor_num =
                  a1 + E_COM * (a2 + E_COM * (a3 + E_COM * (a4 + E_COM * a5)));
              Nuclear_factor_dem =
                  1.0 + E_COM * (b1 + E_COM * (b2 + E_COM * (b3 + E_COM * b4)));
              Nuclear_factor = Nuclear_factor_num / Nuclear_factor_dem;

              // Calculate sigma_DT here:

              if (E_COM != 0)
                cross_section = Nuclear_factor * exp(-B_G / sqrt(E_COM)) /
                                E_COM; // millibarn
              else
                cross_section = 0;

              cross_section *= 1e-27; // converted from mbarn to cm^2

              // Calculate the reactivity here:

              f1 = wts[sp][i] * wts[sp][j] * wts[sp][k];
              f2 = wts[sp2][l] * wts[sp2][m] * wts[sp2][n];

              result += g * cross_section * f1 * in[k + Nv * (j + Nv * i)] *
                        f2 * in2[n + Nv * (m + Nv * l)];
            }

  return result;
}

double GetTNB_tt(double mu, double *in, double *c1, int sp, int sp2) {
  // TNB constants
  double A1 = 38.39;
  double A2 = 448.;
  double A3 = 1.02e-3;
  double A4 = 2.09;
  double A5 = 0.;

  double result = 0.0;
  double g1, g2, g3, g;
  double f2;
  double E_COM;
  double Nuclear_factor;
  double Nuclear_factor_num;
  double Nuclear_factor_dem;
  double cross_section;

  // check to see that we have TT
  if (mu > 2.6e-24 || mu < 2.3e-24) {
    result += 0.;

    return result;
  }

  int i, j, k;
#pragma omp parallel for reduction(+ : result) private(                        \
    i, j, k, g1, g2, g3, g, E_COM, Nuclear_factor_num, Nuclear_factor_dem,     \
    Nuclear_factor, cross_section, f2)

  for (i = 0; i < Nv; i++)
    for (j = 0; j < Nv; j++)
      for (k = 0; k < Nv; k++) {

        g1 = (c1[0] - c[sp2][i]);
        g2 = (c1[1] - c[sp2][j]);
        g3 = (c1[2] - c[sp2][k]);

        g = sqrt(pow(g1, 2) + pow(g2, 2) + pow(g3, 2));
        //                double mu=m[i]*m2[j]/(m[i]+m2[j]);
        E_COM = 0.5 * mu * pow(g, 2) * ERG_TO_EV_CGS *
                1e-3; // Energy of the center-of-mass in keV

        Nuclear_factor_num = A5 + A2 / (pow(A4 - A3 * E_COM, 2) + 1.0);
        Nuclear_factor_dem = E_COM * (exp(-A1 / sqrt(E_COM)) - 1.0);

        // Calculate sigma_TT here:

        cross_section = 1e-24 * Nuclear_factor_num /
                        Nuclear_factor_dem; // convert from bar to cm^2

        // Calculate the reactivity here:
        f2 = wts[sp2][i] * wts[sp2][j] * wts[sp2][k];

        result += g * cross_section * f2 * in[k + Nv * (j + Nv * i)];
      }

  return result;
}

double GetReactivity_tt(double mu, double *in, double *in2, int sp, int sp2) {
  // TNB constants
  double A1 = 38.39;
  double A2 = 448.;
  double A3 = 1.02e-3;
  double A4 = 2.09;
  double A5 = 0.;

  double result = 0.0;
  double g1, g2, g3, g;
  double f1, f2;
  double E_COM;
  double Nuclear_factor;
  double Nuclear_factor_num;
  double Nuclear_factor_dem;
  double cross_section;

  // check to see that we have TT
  if (mu > 2.6e-24 || mu < 2.3e-24) {
    result += 0.;

    return result;
  }

  int i, j, k, l, m, n;
#pragma omp parallel for reduction(+ : result) private(                        \
    i, j, k, l, m, n, g1, g2, g3, g, E_COM, Nuclear_factor_num,                \
    Nuclear_factor_dem, Nuclear_factor, cross_section, f2)
  for (i = 0; i < Nv; i++)
    for (j = 0; j < Nv; j++)
      for (k = 0; k < Nv; k++)
        for (l = 0; l < Nv; l++)
          for (m = 0; m < Nv; m++)
            for (n = 0; n < Nv; n++) {

              g1 = (c[sp][i] - c[sp2][l]);
              g2 = (c[sp][j] - c[sp2][m]);
              g3 = (c[sp][k] - c[sp2][n]);

              g = sqrt(pow(g1, 2) + pow(g2, 2) + pow(g3, 2));
              //                double mu=m[i]*m2[j]/(m[i]+m2[j]);

              E_COM = 0.5 * mu * pow(g, 2) * ERG_TO_EV_CGS *
                      1e-3; // Energy of the center-of-mass in keV

              Nuclear_factor_num = A5 + A2 / (pow(A4 - A3 * E_COM, 2) + 1.0);
              Nuclear_factor_dem = E_COM * (exp(-A1 / sqrt(E_COM)) - 1.0);

              // Calculate sigma_TT here:

              cross_section = 1e-24 * Nuclear_factor_num /
                              Nuclear_factor_dem; // convert from bar to cm^2

              // Calculate the reactivity here:

              f1 = wts[sp][i] * wts[sp][j] * wts[sp][k];
              f2 = wts[sp2][l] * wts[sp2][m] * wts[sp2][n];

              result += g * cross_section * f1 * in[k + Nv * (j + Nv * i)] *
                        f2 * in2[n + Nv * (m + Nv * l)];
            }

  return result;
}

void TNB_DD(double **f, double **f_out, int sp, int rank, int TNBFlag,
            double mu, double *n, double **v, double *T) {

  double R_BGK_DD_HE, R_BGK_DD_T, R_BGK_DD;
  char buffer[50];
  double c1_TNB[3];
  FILE *fpii;

  if (n[sp] > EPS_COLL) {
    R_BGK_DD_HE = GetReactivity_dd_He(mu, f[sp], f[sp], sp, sp);
    R_BGK_DD_T = GetReactivity_dd_T(mu, f[sp], f[sp], sp, sp);
    R_BGK_DD = R_BGK_DD_HE + R_BGK_DD_T;
    printf("DD Reactivity: %g - n %g v %g, T %g \n", R_BGK_DD, n[sp], v[sp][0],
           T[sp]);

    if (TNBFlag == 2) {
      for (int vx = 0; vx < Nv; vx++)
        for (int vy = 0; vy < Nv; vy++)
          for (int vz = 0; vz < Nv; vz++) {
            int index = vz + Nv * (vy + Nv * vx);
            c1_TNB[0] = c[sp][vx];
            c1_TNB[0] = c[sp][vy];
            c1_TNB[0] = c[sp][vz];

            // tail depletion
            f_out[sp][index] -=
                f[sp][index] * GetTNB_dd_He(mu, f[sp], c1_TNB, sp, sp);
            f_out[sp][index] -=
                f[sp][index] * GetTNB_dd_T(mu, f[sp], c1_TNB, sp, sp);
          }
    }

    sprintf(buffer, "Data/TNB_DD_%d.dat", rank);
    if (first_DD) {
      fpii = fopen(buffer, "w");
      first_DD = 0;
    } else
      fpii = fopen(buffer, "a");
    fprintf(fpii, "%5.2e %5.2e %10.6e %10.6e\n", T[sp], T[sp], R_BGK_DD_HE,
            R_BGK_DD_T);
    fclose(fpii);
  } else {
    if (first_DD) {
      fpii = fopen(buffer, "w");
      first_DD = 0;
    } else
      fpii = fopen(buffer, "a");

    fprintf(fpii, "%5.2e %5.2e %10.6e %10.6e\n", T[sp], T[sp], 0.0, 0.0);
    fclose(fpii);
  }
}

void TNB_DT(double **f, double **f_out, int sp, int sp2, int rank, int TNBFlag,
            double mu, double *n, double **v, double *T) {

  char buffer[50];
  double c1_TNB[3];
  double deplete = 0;

  FILE *fpij;

  if ((n[sp] > EPS_COLL) && (n[sp2] > EPS_COLL)) {
    double R_BGK_DT = GetReactivity_dt(mu, f[sp], f[sp2], sp, sp2);

    printf("DT Reactivity: %g\n", R_BGK_DT);

    if (TNBFlag == 2) {
      for (int vx = 0; vx < Nv; vx++)
        for (int vy = 0; vy < Nv; vy++)
          for (int vz = 0; vz < Nv; vz++) {
            int index = vz + Nv * (vy + Nv * vx);
            c1_TNB[0] = c[sp][vx];
            c1_TNB[1] = c[sp][vy];
            c1_TNB[2] = c[sp][vz];

            // tail depletion of D and T
            deplete = GetTNB_dt(mu, f[sp2], c1_TNB, sp, sp2);
            f_out[sp][index] -= f[sp][index] * deplete;
            f_out[sp2][index] -= f[sp2][index] * deplete;
          }
    }

    sprintf(buffer, "Data/TNB_DT_%d.dat", rank);
    if (first_DT) {
      fpij = fopen(buffer, "w");
      first_DT = 0;
    } else
      fpij = fopen(buffer, "a");

    fprintf(fpij, "%5.2e %5.2e %10.6e \n", T[sp], T[sp2], R_BGK_DT);
    fclose(fpij);
  } else {
    if (first_DT) {
      fpij = fopen(buffer, "w");
      first_DT = 0;
    } else
      fpij = fopen(buffer, "a");
    fprintf(fpij, "%5.2e %5.2e %10.6e %10.6e\n", T[sp], T[sp2], 0.0, 0.0);
    fclose(fpij);
  }
}
