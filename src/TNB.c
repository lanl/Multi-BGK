#include "units/unit_data.c"
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

static int Nv;
static double **c;
static double **wts;

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
  double g1, g2, g3, g;
  double f1, f2;
  double E_COM;
  double Nuclear_factor;
  double Nuclear_factor_num;
  double Nuclear_factor_dem;
  double cross_section;
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

        g1 = (c1[0] - c[sp2][i]); //*wts[sp][i]*wts[sp][j]*wts[sp][k];  //
                                  // v_D-v_T=relative velocity
        g2 = (c1[1] - c[sp2][j]); //*wts[sp][i]*wts[sp][j]*wts[sp][k];
        g3 = (c1[2] - c[sp2][k]); //*wts[sp][i]*wts[sp][j]*wts[sp][k];

        g = sqrt(
            pow(g1, 2) + pow(g2, 2) +
            pow(g3, 2)); // Calculate the norm of the relative velocity here.
                         //                double mu=m[i]*m2[j]/(m[i]+m2[j]);
        E_COM = 0.5 * mu * pow(g, 2) * ERG_TO_EV_CGS *
                1e-3; // Energy of the center-of-mass in keV

        Nuclear_factor_num =
            a1 + E_COM * (a2 + E_COM * (a3 + E_COM * (a4 + E_COM * a5)));
        Nuclear_factor_dem =
            1.0 + E_COM * (b1 + E_COM * (b2 + E_COM * (b3 + E_COM * b4)));
        Nuclear_factor = Nuclear_factor_num / Nuclear_factor_dem;

        // Calculate sigma_DT here:

        cross_section = Nuclear_factor * exp(-B_G / sqrt(E_COM)) / E_COM;
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

              g1 = (c[sp][i] - c[sp2][l]); //*wts[sp][i]*wts[sp][j]*wts[sp][k];
                                           //// v_D-v_T=relative velocity
              g2 = (c[sp][j] - c[sp2][m]); //*wts[sp][i]*wts[sp][j]*wts[sp][k];
              g3 = (c[sp][k] - c[sp2][n]); //*wts[sp][i]*wts[sp][j]*wts[sp][k];

              g = sqrt(
                  pow(g1, 2) + pow(g2, 2) +
                  pow(g3,
                      2)); // Calculate the norm of the relative velocity here.
              //                double mu=m[i]*m2[j]/(m[i]+m2[j]);

              E_COM = 0.5 * mu * pow(g, 2) * ERG_TO_EV_CGS *
                      1e-3; // Energy of the center-of-mass in keV

              Nuclear_factor_num =
                  a1 + E_COM * (a2 + E_COM * (a3 + E_COM * (a4 + E_COM * a5)));
              Nuclear_factor_dem =
                  1.0 + E_COM * (b1 + E_COM * (b2 + E_COM * (b3 + E_COM * b4)));
              Nuclear_factor = Nuclear_factor_num / Nuclear_factor_dem;

              // Calculate sigma_DT here:

              cross_section =
                  Nuclear_factor * exp(-B_G / sqrt(E_COM)) / E_COM; // millibarn
              cross_section *= 1e-27; // converted from mbarn to cm^2

              // Calculate the reactivity here:

              f1 = wts[sp][i] * wts[sp][j] * wts[sp][k];    // s^3 / cm^3
              f2 = wts[sp2][l] * wts[sp2][m] * wts[sp2][n]; // s^3 / cm^3
              result += g * cross_section * f1 * in[k + Nv * (j + Nv * i)] *
                        f2 * in2[n + Nv * (m + Nv * l)];
              //       cm/s cm^2           /cm^3                    /cm^3
              //     -> 1 / s cm^3
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

  int i, j, k;
#pragma omp parallel for reduction(+ : result) private(                        \
    i, j, k, g1, g2, g3, g, E_COM, Nuclear_factor_num, Nuclear_factor_dem,     \
    Nuclear_factor, cross_section, f2)

  for (i = 0; i < Nv; i++)
    for (j = 0; j < Nv; j++)
      for (k = 0; k < Nv; k++) {

        g1 = (c1[0] - 0. * c[sp2][i]); //*wts[sp][i]*wts[sp][j]*wts[sp][k];  //
                                       // v_D-v_T=relative velocity
        g2 = (c1[1] - 0. * c[sp2][j]); //*wts[sp][i]*wts[sp][j]*wts[sp][k];
        g3 = (c1[2] - 0. * c[sp2][k]); //*wts[sp][i]*wts[sp][j]*wts[sp][k];

        g = sqrt(
            pow(g1, 2) + pow(g2, 2) +
            pow(g3, 2)); // Calculate the norm of the relative velocity here.
        //                double mu=m[i]*m2[j]/(m[i]+m2[j]);
        E_COM = 0.5 * mu * pow(g, 2) * ERG_TO_EV_CGS *
                1e-3; // Energy of the center-of-mass in keV

        Nuclear_factor_num =
            a1 + E_COM * (a2 + E_COM * (a3 + E_COM * (a4 + E_COM * a5)));
        Nuclear_factor_dem =
            1.0 + E_COM * (b1 + E_COM * (b2 + E_COM * (b3 + E_COM * b4)));
        Nuclear_factor = Nuclear_factor_num / Nuclear_factor_dem;

        // Calculate sigma_DT here:

        cross_section = Nuclear_factor * exp(-B_G / sqrt(E_COM)) / E_COM;
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

  int i, j, k;
#pragma omp parallel for reduction(+ : result) private(                        \
    i, j, k, g1, g2, g3, g, E_COM, Nuclear_factor_num, Nuclear_factor_dem,     \
    Nuclear_factor, cross_section, f2)

  for (i = 0; i < Nv; i++)
    for (j = 0; j < Nv; j++)
      for (k = 0; k < Nv; k++) {

        g1 = (c1[0] - 0. * c[sp2][i]); //*wts[sp][i]*wts[sp][j]*wts[sp][k];  //
                                       // v_D-v_T=relative velocity
        g2 = (c1[1] - 0. * c[sp2][j]); //*wts[sp][i]*wts[sp][j]*wts[sp][k];
        g3 = (c1[2] - 0. * c[sp2][k]); //*wts[sp][i]*wts[sp][j]*wts[sp][k];

        g = sqrt(
            pow(g1, 2) + pow(g2, 2) +
            pow(g3, 2)); // Calculate the norm of the relative velocity here.
        //                double mu=m[i]*m2[j]/(m[i]+m2[j]);
        E_COM = 0.5 * mu * pow(g, 2) * ERG_TO_EV_CGS *
                1e-3; // Energy of the center-of-mass in keV

        Nuclear_factor_num =
            a1 + E_COM * (a2 + E_COM * (a3 + E_COM * (a4 + E_COM * a5)));
        Nuclear_factor_dem =
            1.0 + E_COM * (b1 + E_COM * (b2 + E_COM * (b3 + E_COM * b4)));
        Nuclear_factor = Nuclear_factor_num / Nuclear_factor_dem;

        // Calculate sigma_DT here:

        cross_section = Nuclear_factor * exp(-B_G / sqrt(E_COM)) / E_COM;
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
    i, j, k, g1, g2, g3, g, E_COM, Nuclear_factor_num, Nuclear_factor_dem,     \
    Nuclear_factor, cross_section, f2)

  for (i = 0; i < Nv; i++)
    for (j = 0; j < Nv; j++)
      for (k = 0; k < Nv; k++)
        for (l = 0; l < Nv; l++)
          for (m = 0; m < Nv; m++)
            for (n = 0; n < Nv; n++) {

              g1 = (c[sp][i] -
                    0. * c[sp2][l]); //*wts[sp][i]*wts[sp][j]*wts[sp][k];  //
                                     // v_D-v_T=relative velocity
              g2 = (c[sp][j] -
                    0. * c[sp2][m]); //*wts[sp][i]*wts[sp][j]*wts[sp][k];
              g3 = (c[sp][k] -
                    0. * c[sp2][n]); //*wts[sp][i]*wts[sp][j]*wts[sp][k];

              g = sqrt(
                  pow(g1, 2) + pow(g2, 2) +
                  pow(g3,
                      2)); // Calculate the norm of the relative velocity here.
              //                double mu=m[i]*m2[j]/(m[i]+m2[j]);

              E_COM = 0.5 * mu * pow(g, 2) * ERG_TO_EV_CGS *
                      1e-3; // Energy of the center-of-mass in keV

              Nuclear_factor_num =
                  a1 + E_COM * (a2 + E_COM * (a3 + E_COM * (a4 + E_COM * a5)));
              Nuclear_factor_dem =
                  1.0 + E_COM * (b1 + E_COM * (b2 + E_COM * (b3 + E_COM * b4)));
              Nuclear_factor = Nuclear_factor_num / Nuclear_factor_dem;

              // Calculate sigma_DT here:

              cross_section =
                  Nuclear_factor * exp(-B_G / sqrt(E_COM)) / E_COM; // millibarn
              cross_section *= 1e-27; // converted from mbarn to cm^2

              // Calculate the reactivity here:

              f1 = wts[sp][i] * wts[sp][j] * wts[sp][k];    // s^3 / cm^3
              f2 = wts[sp2][l] * wts[sp2][m] * wts[sp2][n]; // s^3 / cm^3

              result += g * cross_section * f1 * in[k + Nv * (j + Nv * i)] *
                        f2 * in2[n + Nv * (m + Nv * l)];
              //       cm/s cm^2           /cm^3                    /cm^3
              //     -> 1 / s cm^3
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
    i, j, k, g1, g2, g3, g, E_COM, Nuclear_factor_num, Nuclear_factor_dem,     \
    Nuclear_factor, cross_section, f2)

  for (i = 0; i < Nv; i++)
    for (j = 0; j < Nv; j++)
      for (k = 0; k < Nv; k++)
        for (l = 0; l < Nv; l++)
          for (m = 0; m < Nv; m++)
            for (n = 0; n < Nv; n++) {

              g1 = (c[sp][i] -
                    0. * c[sp2][l]); //*wts[sp][i]*wts[sp][j]*wts[sp][k];  //
                                     // v_D-v_T=relative velocity
              g2 = (c[sp][j] -
                    0. * c[sp2][m]); //*wts[sp][i]*wts[sp][j]*wts[sp][k];
              g3 = (c[sp][k] -
                    0. * c[sp2][n]); //*wts[sp][i]*wts[sp][j]*wts[sp][k];

              g = sqrt(
                  pow(g1, 2) + pow(g2, 2) +
                  pow(g3,
                      2)); // Calculate the norm of the relative velocity here.
              //                double mu=m[i]*m2[j]/(m[i]+m2[j]);

              E_COM = 0.5 * mu * pow(g, 2) * ERG_TO_EV_CGS *
                      1e-3; // Energy of the center-of-mass in keV

              Nuclear_factor_num =
                  a1 + E_COM * (a2 + E_COM * (a3 + E_COM * (a4 + E_COM * a5)));
              Nuclear_factor_dem =
                  1.0 + E_COM * (b1 + E_COM * (b2 + E_COM * (b3 + E_COM * b4)));
              Nuclear_factor = Nuclear_factor_num / Nuclear_factor_dem;

              // Calculate sigma_DT here:

              cross_section =
                  Nuclear_factor * exp(-B_G / sqrt(E_COM)) / E_COM; // millibarn
              cross_section *= 1e-27; // converted from mbarn to cm^2

              // Calculate the reactivity here:

              f1 = wts[sp][i] * wts[sp][j] * wts[sp][k];    // s^3 / cm^3
              f2 = wts[sp2][l] * wts[sp2][m] * wts[sp2][n]; // s^3 / cm^3

              result += g * cross_section * f1 * in[k + Nv * (j + Nv * i)] *
                        f2 * in2[n + Nv * (m + Nv * l)];
              //       cm/s cm^2           /cm^3                    /cm^3
              //     -> 1 / s cm^3
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

  int i, j, k;
#pragma omp parallel for reduction(+ : result) private(                        \
    i, j, k, g1, g2, g3, g, E_COM, Nuclear_factor_num, Nuclear_factor_dem,     \
    Nuclear_factor, cross_section, f2)

  for (i = 0; i < Nv; i++)
    for (j = 0; j < Nv; j++)
      for (k = 0; k < Nv; k++) {

        g1 = (c1[0] - 0. * c[sp2][i]); //*wts[sp][i]*wts[sp][j]*wts[sp][k];  //
                                       // v_D-v_T=relative velocity
        g2 = (c1[1] - 0. * c[sp2][j]); //*wts[sp][i]*wts[sp][j]*wts[sp][k];
        g3 = (c1[2] - 0. * c[sp2][k]); //*wts[sp][i]*wts[sp][j]*wts[sp][k];

        g = sqrt(
            pow(g1, 2) + pow(g2, 2) +
            pow(g3, 2)); // Calculate the norm of the relative velocity here.
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
    i, j, k, g1, g2, g3, g, E_COM, Nuclear_factor_num, Nuclear_factor_dem,     \
    Nuclear_factor, cross_section, f2)

  for (i = 0; i < Nv; i++)
    for (j = 0; j < Nv; j++)
      for (k = 0; k < Nv; k++)
        for (l = 0; l < Nv; l++)
          for (m = 0; m < Nv; m++)
            for (n = 0; n < Nv; n++) {

              g1 = (c[sp][i] -
                    0. * c[sp2][l]); //*wts[sp][i]*wts[sp][j]*wts[sp][k];  //
                                     // v_D-v_T=relative velocity
              g2 = (c[sp][j] -
                    0. * c[sp2][m]); //*wts[sp][i]*wts[sp][j]*wts[sp][k];
              g3 = (c[sp][k] -
                    0. * c[sp2][n]); //*wts[sp][i]*wts[sp][j]*wts[sp][k];

              g = sqrt(
                  pow(g1, 2) + pow(g2, 2) +
                  pow(g3,
                      2)); // Calculate the norm of the relative velocity here.
              //                double mu=m[i]*m2[j]/(m[i]+m2[j]);

              E_COM = 0.5 * mu * pow(g, 2) * ERG_TO_EV_CGS *
                      1e-3; // Energy of the center-of-mass in keV

              Nuclear_factor_num = A5 + A2 / (pow(A4 - A3 * E_COM, 2) + 1.0);
              Nuclear_factor_dem = E_COM * (exp(-A1 / sqrt(E_COM)) - 1.0);

              // Calculate sigma_TT here:

              cross_section = 1e-24 * Nuclear_factor_num /
                              Nuclear_factor_dem; // convert from bar to cm^2

              // Calculate the reactivity here:

              f1 = wts[sp][i] * wts[sp][j] * wts[sp][k];    // s^3 / cm^3
              f2 = wts[sp2][l] * wts[sp2][m] * wts[sp2][n]; // s^3 / cm^3

              result += g * cross_section * f1 * in[k + Nv * (j + Nv * i)] *
                        f2 * in2[n + Nv * (m + Nv * l)];
              //       cm/s cm^2           /cm^3                    /cm^3
              //     -> 1 / s cm^3
            }

  return result;
}
