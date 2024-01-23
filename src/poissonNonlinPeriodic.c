#include "poissonNonlinPeriodic.h"
#include "units/unit_data.c"
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_fermi_dirac.h>
#include <gsl/gsl_vector.h>
#include <math.h>

// Adapted from Cory's matlab code

static double absTol = 1e-8;
static double relTol = 1e-8;
static int phiflag = 0;

/***********************************************/

// Declare nonlinear subroutines here for definiton later
void electronSource(gsl_vector *phi, double *g, double *gPrime, double ne0,
                    double *Te);
void electronSource_TF(gsl_vector *phi, double *g, double *gPrime, double mu,
                       double *Te);
double chemPot_TF(double *source, int N, double *Te, double mu0);

// subroutines for checking charge conservation
double chargeTallyYukLin(gsl_vector *phi, double ne0, double Te);
double chargeTallyYukNonlin(gsl_vector *phi, double ne0, double Te);
double chargeTallyFermiLin(gsl_vector *phi, double mu, double Te);
double chargeTallyFermiNonin(gsl_vector *phi, double mu, double Te);

/***********************************************/

/**************************************************************/
/* PoissNonlinPeriodic1D                                      */
/*                                                            */
/* DESCRIPTION                                                */
/* Computes the nonlinear poisson equation                    */
/* based on the Yukawa formulation                            */
/*                                                            */
/* INPUTS                                                     */
/* N: number of cells in domain                               */
/* source: RHS vector, length N                               */
/* dx: size of each cell (fixed)                              */
/* Te: electron temperature (eV)                              */
/*                                                            */
/* OUTPUT                                                     */
/* phi: solution to the equation                              */
/*      (this is actually e phi, converted to ergs)           */
/**************************************************************/

/**************************************************************
The Poisson equation solved by this routine, for reference:

(e phi)_xx = 4 pi e^2 (sum_i Z_i n_i - n_e0 exp(e phi / T_e))

units:
(e phi) : eV
e^2     : eV-cm
n_i     : 1/cm^3
T_e     : eV

**************************************************************/

void PoissNonlinPeriodic1D(int N, double *source, double dx, double Lx,
                           double *phi, double *Te) {

  int i, loop;

  // Set up vectors
  gsl_vector *phiVec = gsl_vector_calloc(N);
  gsl_vector *phiNext = gsl_vector_calloc(N);
  gsl_vector *dphi = gsl_vector_calloc(N);
  gsl_vector *RHS = gsl_vector_calloc(N);
  gsl_vector *phixx = gsl_vector_calloc(N);
  double *g = malloc(N * sizeof(double));
  double *gPrime = malloc(N * sizeof(double));

  // initialize
  if (phiflag == 0) {
    phiflag = 1;
    for (i = 0; i < N; i++)
      gsl_vector_set(phiVec, i, 1e-10);
  } else // use previous field as guess
    for (i = 0; i < N; i++)
      gsl_vector_set(phiVec, i, phi[i]);

  // Set up matrices
  gsl_matrix *A = gsl_matrix_calloc(N, N); //-d_xx operator
  gsl_matrix *B = gsl_matrix_calloc(N, N); // Jacobian matrix

  // initialize the matrix A
  gsl_matrix_set(A, 0, 0, 2.0 / dx / dx);
  gsl_matrix_set(A, 0, 1, -1.0 / dx / dx);
  gsl_matrix_set(A, 0, N - 1, -1.0 / dx / dx);
  for (i = 1; i < (N - 1); i++) {
    gsl_matrix_set(A, i, i, 2.0 / dx / dx);
    gsl_matrix_set(A, i, i - 1, -1.0 / dx / dx);
    gsl_matrix_set(A, i, i + 1, -1.0 / dx / dx);
  }
  gsl_matrix_set(A, N - 1, N - 1, 2.0 / dx / dx);
  gsl_matrix_set(A, N - 1, N - 2, -1.0 / dx / dx);
  gsl_matrix_set(A, N - 1, 0, -1.0 / dx / dx);

  gsl_matrix_memcpy(B, A);

  // set up tolerance stuff
  double relErr = relTol + 1.0;
  double absErr = absTol + 1.0;

  // Done with setup, start the iteration
  gsl_permutation *Permut = gsl_permutation_alloc(N);
  int signum;

  // Find total charge - for setting n_e0
  double Ci = 0;

  for (i = 0; i < N; i++)
    Ci += source[i];

  // distribute charge evenly over each cell
  double ne0 = (Ci * dx / Lx);

  loop = 0;

  /*
  while ( (( relErr > relTol ) || ( absErr > absTol )) && (loop < 50) ) {

    //Calculates exp(e phi/T_e) and its derivative
    electronSource(phiVec,g,gPrime,ne0,Te);

    //Multiplies A*(e phi) and stores in phixx
    gsl_blas_dgemv(CblasNoTrans,1.0,A,phiVec,0.0,phixx);

    //finishes setting up Jacobian matrix with the nonlinear part
    //Also sets up RHS of the linear solve
    for(i=0;i<N;i++) {
      gsl_matrix_set(B,i,i,gsl_matrix_get(A,i,i) + 4.0*M_PI*E_02_CGS*gPrime[i]);
      gsl_vector_set(RHS,i, -gsl_vector_get(phixx,i) +
  4.0*M_PI*E_02_CGS*(source[i] - g[i]));
    }

    //solves B*dphi = RHS
    gsl_linalg_LU_decomp(B,Permut,&signum);
    gsl_linalg_LU_solve(B,Permut,RHS,dphi);

    //check to see if we have converged or broke something
    relErr = gsl_blas_dnrm2(dphi)/gsl_blas_dnrm2(phiVec);
    absErr = gsl_blas_dnrm2(dphi);

    if(isnan(absErr) != 0) {
      printf("Error: nan detected in PoissonNonlin\n");
      printf("g              gp         phi \n");
      for(i=0;i<N;i++)
        printf("%le  %le  %le \n",g[i],gPrime[i],gsl_vector_get(phiVec,i));
      exit(1);
    }


    //now we know the update, set phinext

    gsl_blas_daxpy(1.0,dphi,phiVec);
    //translating from BLAS jargon, this sets phi^(n+1) = dphi + phi^n
    loop++;
  }
  */

  /* Faster version suggested by Cory */

  while (((relErr > relTol) || (absErr > absTol)) && (loop < 50)) {
    electronSource(phiVec, g, gPrime, ne0, Te);

    //get a fresh B to get off diagonal terms and remove changes done by LU decomp
    gsl_matrix_memcpy(B, A);    
    
    for (i = 0; i < N; i++) {
      gsl_matrix_set(
          B, i, i, gsl_matrix_get(A, i, i) + 4.0 * M_PI * E_02_CGS * gPrime[i]);
      gsl_vector_set(
          RHS, i,
          4.0 * M_PI * E_02_CGS *
              (source[i] - g[i] + gPrime[i] * gsl_vector_get(phiVec, i)));
    }
    // solves B*phiNext = RHS
    gsl_linalg_LU_decomp(B, Permut, &signum);
    gsl_linalg_LU_solve(B, Permut, RHS, phiNext);

    // check error tol
    gsl_vector_memcpy(dphi, phiNext);
    gsl_blas_daxpy(-1.0, phiVec, dphi);
    relErr = gsl_blas_dnrm2(dphi) / gsl_blas_dnrm2(phiVec);
    absErr = gsl_blas_dnrm2(dphi);
    if (isnan(absErr) != 0) {
      printf("Error: nan detected in PoissonNonlin\n");
      printf("g              gp         phi \n");
      for (i = 0; i < N; i++)
        printf("%le  %le  %le \n", g[i], gPrime[i], gsl_vector_get(phiVec, i));
      exit(1);
    }

    gsl_vector_memcpy(phiVec, phiNext);
    loop++;
  }

  if (loop == 50)
    printf("Newton method failed to converge\n");

  double ion_tot = 0.0;
  double e_tot = 0.0;

  for (i = 0; i < N; i++) {
    ion_tot += source[i] * dx;
    e_tot += ne0 * exp(gsl_vector_get(phiVec, i) / Te[i]) * dx;
  }
  // printf("Charges:| ion: %g  electron %g abs diff %g reldiff %g
  // \n",ion_tot,e_tot,fabs(ion_tot-e_tot),fabs(ion_tot-e_tot)/ion_tot);

  // put results in phi, converted to ergs
  for (i = 0; i < N; i++)
    phi[i] = gsl_vector_get(phiVec, i) / ERG_TO_EV_CGS;

  gsl_vector_free(phiVec);
  gsl_vector_free(phiNext);
  gsl_vector_free(RHS);
  gsl_vector_free(dphi);
  gsl_vector_free(phixx);
  gsl_matrix_free(A);
  gsl_matrix_free(B);
  gsl_permutation_free(Permut);
  free(g);
  free(gPrime);
}

void electronSource(gsl_vector *phi, double *g, double *gPrime, double ne0,
                    double *Te) {

  int N = (*phi).size;
  int i;

  for (i = 0; i < N; i++) {
    g[i] = ne0 * exp(gsl_vector_get(phi, i) / Te[i]);
    gPrime[i] = (1.0 / Te[i]) * ne0 * exp(gsl_vector_get(phi, i) / Te[i]);
  }
}

/**************************************************************/
/* PoissLinPeriodic1D                                         */
/*                                                            */
/* DESCRIPTION                                                */
/* Computes the linear poisson equation                       */
/* based on the Yukawa formulation                            */
/*                                                            */
/* INPUTS                                                     */
/* N: number of cells in domain                               */
/* source: RHS vector, length N                               */
/* Lx: length of domain                                       */
/* Te: Electron emperature in eV                              */
/*                                                            */
/* OUTPUT                                                     */
/* e phi: solution to the equation (in ergs)                  */
/**************************************************************/

/**************************************************************
The Poisson equation solved by this routine, for reference:

-(e phi)_xx + (4pi e^2 / Te) n_e0 (e phi) = 4pi e^2   (source - n_e0)
eV / cm^2         eV-cm  eV  1/cc (eV   )       eV-cm  1/cc     1/cc

**************************************************************/

void PoissLinPeriodic1D(int N, double *source, double dx, double Lx,
                        double *phi, double *Te) {

  int i;

  // Set up vectors
  gsl_vector *phiVec = gsl_vector_calloc(N);
  gsl_vector *RHS = gsl_vector_calloc(N);

  // initialize phi
  for (i = 0; i < N; i++)
    gsl_vector_set(phiVec, i, 1.0);

  // Set up matrices
  gsl_matrix *A = gsl_matrix_calloc(N, N);

  // Find total charge
  double Ci = 0;

  for (i = 0; i < N; i++)
    Ci += source[i];

  // distribute charge evenly over each cell
  double ne0 = (Ci * dx / Lx);

  double ne0_diag = (4.0 * M_PI * E_02_CGS) * ne0;

  // initialize the matrix A
  gsl_matrix_set(A, 0, 0, 2.0 / dx / dx + ne0_diag / Te[0]);
  gsl_matrix_set(A, 0, 1, -1.0 / dx / dx);
  gsl_matrix_set(A, 0, N - 1, -1.0 / dx / dx);
  for (i = 1; i < (N - 1); i++) {
    gsl_matrix_set(A, i, i, 2.0 / dx / dx + ne0_diag / Te[i]);
    gsl_matrix_set(A, i, i - 1, -1.0 / dx / dx);
    gsl_matrix_set(A, i, i + 1, -1.0 / dx / dx);
  }
  gsl_matrix_set(A, N - 1, N - 1, 2.0 / dx / dx + ne0_diag / Te[N - 1]);
  gsl_matrix_set(A, N - 1, N - 2, -1.0 / dx / dx);
  gsl_matrix_set(A, N - 1, 0, -1.0 / dx / dx);

  // set up RHS for solve
  for (i = 0; i < N; i++)
    gsl_vector_set(RHS, i, 4.0 * M_PI * E_02_CGS * (source[i] - ne0));

  // Set up the LA solve
  gsl_permutation *Permut = gsl_permutation_alloc(N);
  int signum;

  // solves A*phiVec = RHS
  gsl_linalg_LU_decomp(A, Permut, &signum);
  gsl_linalg_LU_solve(A, Permut, RHS, phiVec);

  double phitot = 0.0;

  // check charge conservation
  double ion_tot = 0.0;
  double e_tot = 0.0;

  for (i = 0; i < N; i++) {
    ion_tot += source[i] * dx;
    e_tot += ne0 * (1 + gsl_vector_get(phiVec, i) / Te[i]) * dx;
  }
  // printf("Charges:| ion: %g  electron %g abs diff %g reldiff %g
  // \n",ion_tot,e_tot,fabs(ion_tot-e_tot),fabs(ion_tot-e_tot)/ion_tot);

  // put results in phi and convert to ergs
  for (i = 0; i < N; i++) {
    phi[i] = gsl_vector_get(phiVec, i) / ERG_TO_EV_CGS;
    phitot += phi[i] / Te[i];
  }
  // printf("phitot: %g\n",phitot);

  gsl_vector_free(phiVec);
  gsl_vector_free(RHS);
  gsl_matrix_free(A);
}

/**************************************************************/
/* PoissLinPeriodic1D_TF                                      */
/*                                                            */
/* DESCRIPTION                                                */
/* Computes the linear poisson equation                       */
/* based on the Thomas-Fermi formulation                      */
/*                                                            */
/* INPUTS                                                     */
/* N: number of cells in domain                               */
/* source: RHS vector, length N                               */
/* Lx: length of domain                                       */
/* Te: Electron emperature in eV, length N                    */
/*                                                            */
/* OUTPUT                                                     */
/* e phi: solution to the equation (in ergs)                  */
/**************************************************************/

/**************************************************************
The Poisson equation solved by this routine, for reference:

-(e phi)_xx + [ ( 8pi e^2 (2pim_e)^3/2 Te^1/2 / (2pi hbar)^3) F_(-1/2)(beta mu)
] (e phi)
eV / cm^2     [       eV-cm   g        eV            eV-s              1/eV eV
]  eV

= 4pi e^2   (source - 2(2pi m_e T_e)^(3/2) / (2 pi hbar)^3 F_(1/2)(beta mu) )
      ev-cm  1/cm^3         g   eV                 eV-s            1/eV eV

**************************************************************/

void PoissLinPeriodic1D_TF(int N, double *source, double dx, double Lx,
                           double *phi, double *Te) {
  int i;

  // Set up vectors
  gsl_vector *phiVec = gsl_vector_calloc(N);
  gsl_vector *RHS = gsl_vector_calloc(N);

  // initialize phi
  for (i = 0; i < N; i++)
    gsl_vector_set(phiVec, i, 1.0);

  // Set up matrices
  gsl_matrix *A = gsl_matrix_calloc(N, N);

  // Find average chemical potential
  double mu = chemPot_TF(source, N, Te, 1.0); // eV

  // Fermi integrals
  double F_mhalf;
  double F_half;

  double diag = 8.0 * M_PI * E_02_CGS *
                pow(2.0 * M_PI * M_ELEC_CGS * ERG_TO_EV_CGS, 1.5) /
                pow(2.0 * M_PI * HBAR_CGS, 3);
  // double diag =
  // F_mhalf*8.0*M_PI*E_02_CGS*pow(2.0*M_PI*M_ELEC_CGS*Te*ERG_TO_EV_CGS,1.5)/pow(2.0*M_PI*HBAR_CGS,3)/Te;
  //                             eV cm                g^3/2       eV^3/2 s^3 /
  //                             g^3/2 cm^3         eV^3 s^3

  double fdfac = sqrt(M_PI) * 0.5;

  // initialize the matrix A
  F_mhalf = fdfac * gsl_sf_fermi_dirac_mhalf(mu / Te[0]);
  gsl_matrix_set(A, 0, 0,
                 2.0 / dx / dx + diag * sqrt(Te[0]) * F_mhalf); // 1/cm^2
  gsl_matrix_set(A, 0, 1, -1.0 / dx / dx);
  gsl_matrix_set(A, 0, N - 1, -1.0 / dx / dx);
  for (i = 1; i < (N - 1); i++) {
    F_mhalf = fdfac * gsl_sf_fermi_dirac_mhalf(mu / Te[i]);
    gsl_matrix_set(A, i, i,
                   2.0 / dx / dx + diag * sqrt(Te[i]) * F_mhalf); // 1/cm^2
    gsl_matrix_set(A, i, i - 1, -1.0 / dx / dx);
    gsl_matrix_set(A, i, i + 1, -1.0 / dx / dx);
  }
  F_mhalf = fdfac * gsl_sf_fermi_dirac_mhalf(mu / Te[N - 1]);
  gsl_matrix_set(A, N - 1, N - 1,
                 2.0 / dx / dx + diag * sqrt(Te[N - 1]) * F_mhalf); // 1/cm^2
  gsl_matrix_set(A, N - 1, N - 2, -1.0 / dx / dx);
  gsl_matrix_set(A, N - 1, 0, -1.0 / dx / dx);

  // set up RHS for solve
  for (i = 0; i < N; i++) {
    F_half = fdfac * gsl_sf_fermi_dirac_half(mu / Te[i]);
    gsl_vector_set(
        RHS, i,
        4.0 * M_PI * E_02_CGS *
            (source[i] -
             F_half * 2.0 *
                 pow(2.0 * M_PI * M_ELEC_CGS * Te[i] * ERG_TO_EV_CGS, 1.5) /
                 pow(2.0 * M_PI * HBAR_CGS, 3))); // eV /cm^2
  }
  // Set up the LA solve
  gsl_permutation *Permut = gsl_permutation_alloc(N);
  int signum;

  // solves A*phiVec = RHS
  gsl_linalg_LU_decomp(A, Permut, &signum);
  gsl_linalg_LU_solve(A, Permut, RHS, phiVec);

  double phitot = 0.0;

  double ion_tot = 0.0;
  double e_tot = 0.0;

  for (i = 0; i < N; i++) {
    ion_tot += source[i] * dx;
    e_tot += (F_half + gsl_vector_get(phiVec, i) * F_mhalf / Te[i]) *
             (2.0 * pow(2.0 * M_PI * M_ELEC_CGS * Te[i] * ERG_TO_EV_CGS, 1.5) /
              pow(2.0 * M_PI * HBAR_CGS, 3)) *
             dx;
  }
  // printf("Charges:| ion: %g  electron %g abs diff %g reldiff %g
  // \n",ion_tot,e_tot,fabs(ion_tot-e_tot),fabs(ion_tot-e_tot)/ion_tot);

  // put results in phi and convert to ergs
  for (i = 0; i < N; i++) {
    phi[i] = gsl_vector_get(phiVec, i) / ERG_TO_EV_CGS;
    phitot += phi[i] / Te[i];
  }

  gsl_vector_free(phiVec);
  gsl_vector_free(RHS);
  gsl_matrix_free(A);
}

/**************************************************************/
/* PoissNonLinPeriodic1D_TF                                   */
/*                                                            */
/* DESCRIPTION                                                */
/* Computes the nonlinear poisson equation                    */
/* based on the Thomas-Fermi formulation                      */
/*                                                            */
/* INPUTS                                                     */
/* N: number of cells in domain                               */
/* source: RHS vector, length N                               */
/* Lx: length of domain                                       */
/* Te: Electron emperature in eV                              */
/*                                                            */
/* OUTPUT                                                     */
/* e phi: solution to the equation (in ergs)                  */
/**************************************************************/

/**************************************************************
The Poisson equation solved by this routine, for reference:

-(e phi)_xx = 4pi e^2   (source - 2(2pi m_e T_e)^(3/2) / (2 pi hbar)^3
F_(1/2)(beta (mu + e phi)) )
  eV / cm^2       ev-cm  1/cm^3         g   eV                 eV-s
1/eV eV    eV


**************************************************************/

void PoissNonlinPeriodic1D_TF(int N, double *source, double dx, double Lx,
                              double *phi, double *Te) {

  int i, loop;

  // Set up vectors
  gsl_vector *phiVec = gsl_vector_calloc(N);
  gsl_vector *phiNext = gsl_vector_calloc(N);
  gsl_vector *dphi = gsl_vector_calloc(N);
  gsl_vector *RHS = gsl_vector_calloc(N);
  gsl_vector *phixx = gsl_vector_calloc(N);
  double *g = malloc(N * sizeof(double));
  double *gPrime = malloc(N * sizeof(double));

  double fdfac = 0.5 * sqrt(M_PI);

  // initialize
  if (phiflag == 0) {
    printf("Initializing phi\n");
    phiflag = 1;
    for (i = 0; i < N; i++)
      gsl_vector_set(phiVec, i, 1e1);
  } else // use previous field as guess
    for (i = 0; i < N; i++)
      gsl_vector_set(phiVec, i, phi[i]);

  // Set up matrices
  gsl_matrix *A = gsl_matrix_calloc(N, N); //-d_xx operator
  gsl_matrix *B = gsl_matrix_calloc(N, N); // Jacobian matrix

  // initialize the matrix A
  gsl_matrix_set(A, 0, 0, 2.0 / dx / dx);
  gsl_matrix_set(A, 0, 1, -1.0 / dx / dx);
  gsl_matrix_set(A, 0, N - 1, -1.0 / dx / dx);
  for (i = 1; i < (N - 1); i++) {
    gsl_matrix_set(A, i, i, 2.0 / dx / dx);
    gsl_matrix_set(A, i, i - 1, -1.0 / dx / dx);
    gsl_matrix_set(A, i, i + 1, -1.0 / dx / dx);
  }
  gsl_matrix_set(A, N - 1, N - 1, 2.0 / dx / dx);
  gsl_matrix_set(A, N - 1, N - 2, -1.0 / dx / dx);
  gsl_matrix_set(A, N - 1, 0, -1.0 / dx / dx);

  gsl_matrix_memcpy(B, A);

  // set up tolerance stuff
  double relErr = relTol + 1.0;
  double absErr = absTol + 1.0;

  // Done with setup, start the iteration
  gsl_permutation *Permut = gsl_permutation_alloc(N);
  int signum;

  // Get chemical potential approximation
  double mu = chemPot_TF(source, N, Te, 1.0);

  loop = 0;

  /*
  while ( (( relErr > relTol ) || ( absErr > absTol )) && (loop < 50) ) {

    //Calculates exp(e phi/T_e) and its derivative
    electronSource_TF(phiVec,g,gPrime,mu,Te);

    //Multiplies A*(e phi) and stores in phixx
    gsl_blas_dgemv(CblasNoTrans,1.0,A,phiVec,0.0,phixx);

    //finishes setting up Jacobian matrix with the nonlinear part
    //Also sets up RHS of the linear solve
    for(i=0;i<N;i++) {
      gsl_matrix_set(B,i,i,gsl_matrix_get(A,i,i) + 4.0*M_PI*E_02_CGS*gPrime[i]);
      gsl_vector_set(RHS,i, -gsl_vector_get(phixx,i) +
  4.0*M_PI*E_02_CGS*(source[i] - g[i]));
    }

    //solves B*dphi = RHS
    gsl_linalg_LU_decomp(B,Permut,&signum);
    gsl_linalg_LU_solve(B,Permut,RHS,dphi);

    //check to see if we have converged or broke something
    relErr = gsl_blas_dnrm2(dphi)/gsl_blas_dnrm2(phiVec);
    absErr = gsl_blas_dnrm2(dphi);

    if(isnan(absErr) != 0) {
      printf("Error: nan detected in PoissonNonlin_TF\n");
      printf("g              gp         phi \n");
      for(i=0;i<N;i++)
        printf("%le  %le  %le \n",g[i],gPrime[i],gsl_vector_get(phiVec,i));
      exit(1);
    }


    //now we know the update, set phinext

    gsl_blas_daxpy(1.0,dphi,phiVec);
    //translating from BLAS jargon, this sets phi^(n+1) = dphi + phi^n
    loop++;
  }
  */

  /* Faster version suggested by Cory */

  while (((relErr > relTol) || (absErr > absTol)) && (loop < 50)) {
    electronSource_TF(phiVec, g, gPrime, mu, Te);

    //get a fresh B to get off diagonal terms and remove changes done by LU decomp
    gsl_matrix_memcpy(B, A);    

    for (i = 0; i < N; i++) {
      gsl_matrix_set(
          B, i, i, gsl_matrix_get(A, i, i) + 4.0 * M_PI * E_02_CGS * gPrime[i]);
      gsl_vector_set(
          RHS, i,
          4.0 * M_PI * E_02_CGS *
              (source[i] - g[i] + gPrime[i] * gsl_vector_get(phiVec, i)));
    }
    // solves B*phiNext = RHS
    gsl_linalg_LU_decomp(B, Permut, &signum);
    gsl_linalg_LU_solve(B, Permut, RHS, phiNext);

    // check error tol
    gsl_vector_memcpy(dphi, phiNext);
    gsl_blas_daxpy(-1.0, phiVec, dphi);
    relErr = gsl_blas_dnrm2(dphi) / gsl_blas_dnrm2(phiVec);
    absErr = gsl_blas_dnrm2(dphi);
    if (isnan(absErr) != 0) {
      printf("Error: nan detected in PoissonNonlin\n");
      printf("loop: %d \n", loop);
      printf("g              gp         phi \n");
      for (i = 0; i < N; i++)
        printf("%le  %le  %le \n", g[i], gPrime[i], gsl_vector_get(phiVec, i));
      exit(1);
    }

    gsl_vector_memcpy(phiVec, phiNext);
    loop++;
  }

  if (loop == 50)
    printf("Newton method failed to converge\n");

  double ion_tot = 0.0;
  double e_tot = 0.0;

  for (i = 0; i < N; i++) {
    ion_tot += source[i] * dx;
    e_tot += fdfac *
             gsl_sf_fermi_dirac_half((mu + gsl_vector_get(phiVec, i)) / Te[i]) *
             (2.0 * pow(2.0 * M_PI * M_ELEC_CGS * Te[i] * ERG_TO_EV_CGS, 1.5) /
              pow(2.0 * M_PI * HBAR_CGS, 3)) *
             dx;
  }
  // printf("Charges:| ion: %g  electron %g abs diff %g reldiff %g
  // \n",ion_tot,e_tot,fabs(ion_tot-e_tot),fabs(ion_tot-e_tot)/ion_tot);

  // put results in phi, converted to ergs
  for (i = 0; i < N; i++)
    phi[i] = gsl_vector_get(phiVec, i) / ERG_TO_EV_CGS;

  gsl_vector_free(phiVec);
  gsl_vector_free(phiNext);
  gsl_vector_free(RHS);
  gsl_vector_free(dphi);
  gsl_vector_free(phixx);
  gsl_matrix_free(A);
  gsl_matrix_free(B);
  gsl_permutation_free(Permut);
  free(g);
  free(gPrime);
}

void simplePoisson(int N, double *source, double dx, double Lx, double *phi) {
  int i;

  // Set up vectors
  gsl_vector *phiVec = gsl_vector_calloc(N);
  gsl_vector *RHS = gsl_vector_calloc(N);

  // initialize phi
  for (i = 0; i < N; i++)
    gsl_vector_set(phiVec, i, 1.0);

  // Set up matrices
  gsl_matrix *A = gsl_matrix_calloc(N, N);

  // initialize the matrix A
  gsl_matrix_set(A, 0, 0, 2.0 / dx / dx);
  gsl_matrix_set(A, 0, 1, -1.0 / dx / dx);
  gsl_matrix_set(A, 0, N - 1, -1.0 / dx / dx);
  for (i = 1; i < (N - 1); i++) {
    gsl_matrix_set(A, i, i, 2.0 / dx / dx);
    gsl_matrix_set(A, i, i - 1, -1.0 / dx / dx);
    gsl_matrix_set(A, i, i + 1, -1.0 / dx / dx);
  }
  for (i = 0; i < N; i++) // set last row to 1 to avoid constant solution issue
    gsl_matrix_set(A, N - 1, i, 1.0);

  // set up RHS for solve
  for (i = 0; i < N; i++)
    gsl_vector_set(RHS, i, 4.0 * M_PI * E_02_CGS * source[i]);
  gsl_vector_set(RHS, N - 1, 0);

  // Set up the LA solve
  gsl_permutation *Permut = gsl_permutation_alloc(N);
  int signum;

  // solves A*phiVec = RHS
  gsl_linalg_LU_decomp(A, Permut, &signum);
  gsl_linalg_LU_solve(A, Permut, RHS, phiVec);

  // put results in phi and convert to ergs
  for (i = 0; i < N; i++) {
    phi[i] = gsl_vector_get(phiVec, i) / ERG_TO_EV_CGS;
  }

  gsl_vector_free(phiVec);
  gsl_vector_free(RHS);
  gsl_matrix_free(A);
}

// Calculates the average value of the chemical potential for all cells
// by inverting the Fermi integral

double chemPot_TF(double *source, int N, double *Te, double mu0) {
  double RHS;
  double RHS_prefac = 0.5 * pow(2.0 * M_PI * HBAR_CGS, 3) /
                      pow(2.0 * M_PI * M_ELEC_CGS * ERG_TO_EV_CGS, 1.5);
  //                                   ev^3 s^3                 g^3/2 -> eV^3/2
  //                                   * cm^3

  double x0;
  double xn, xnp1;
  double absErr, relErr;

  double musum = 0.0;

  int i, loop;

  double fdfac = 0.5 * sqrt(M_PI);

  for (i = 0; i < N; i++) {
    x0 = mu0 / Te[i]; // eV / eV = unitless
    absErr = 1.0 + absTol;
    relErr = 1.0 + relTol;

    RHS = RHS_prefac * source[i] / pow(Te[i], 1.5); //[source] * cm^3 = unitless

    loop = 0;
    xn = x0;

    while (((absErr > absTol) || (relErr > relTol)) && (loop < 50)) {
      xnp1 = xn - (fdfac * gsl_sf_fermi_dirac_half(xn) - RHS) /
                      (fdfac * gsl_sf_fermi_dirac_mhalf(xn));
      absErr = fabs(xnp1 - xn);
      relErr = fabs((xnp1 - xn)) / fabs(xn);
      xn = xnp1;
      loop++;
    }

    // printf("i %d loop %d mu %g abs %g rel %g\n",i,loop,xn, absErr, relErr);
    musum += xn * Te[i];
  }

  return musum / (double)N;
}

void electronSource_TF(gsl_vector *phi, double *g, double *gPrime, double mu,
                       double *Te) {

  int N = (*phi).size;
  int i;

  double fdfac = 0.5 * sqrt(M_PI);

  double prefac = fdfac * 2.0 *
                  pow(2.0 * M_PI * M_ELEC_CGS * ERG_TO_EV_CGS, 1.5) /
                  pow(2.0 * M_PI * HBAR_CGS, 3);

  for (i = 0; i < N; i++) {
    g[i] = prefac * pow(Te[i], 1.5) *
           gsl_sf_fermi_dirac_half((mu + gsl_vector_get(phi, i)) / Te[i]);
    gPrime[i] = (1.0 / Te[i]) * prefac * pow(Te[i], 1.5) *
                gsl_sf_fermi_dirac_mhalf((mu + gsl_vector_get(phi, i)) / Te[i]);
  }
}

// returns a uniform electron temperature on the grid
void get_uniform_Te(double *Te, int Nx, double T0) {
  int i;
  for (i = 0; i < Nx; i++)
    Te[i] = T0;
}

// linearly ramps electron temperature from T_start at time zero to T_end at
// time tfinal
void get_ramp_Te(double *Te, int Nx, double T_start, double T_end, double t,
                 double tfinal) {
  int l;
  double rampTemp;

  // calculate cutoff time
  if (t > tfinal)
    rampTemp = T_end;
  else
    rampTemp = T_start + t * (T_end - T_start) / tfinal;

  for (l = 0; l < Nx; l++) {
    Te[l] = rampTemp;
  }
}

void get_ramp_Te_cubic(double *Te, int Nx, double alpha, double shift,
                       double t) {

  for (int l = 0; l < Nx; ++l) {
    Te[l] = alpha * pow((t + shift) / 1e-9, 3);
  }
}
