#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <gmp.h>
#include <mps/mps.h>
#include <omp.h>

/*Function to calculate the Polynomial Coefficients*/
void coefi_pol(int l, mpf_t al_real[], mpf_t al_imag[], mpf_t coef_real[],
               mpf_t coef_imag[], int LMAX) {
  mpz_t int_binom;
  mpf_t float_binom, raiz_binom, coef_negativ;
  int index;

  // Initialize GMP variables
  mpz_init(int_binom);
  mpf_inits(float_binom, raiz_binom, coef_negativ, NULL);

  // Initialize coefficient arrays
  mpf_inits(coef_real[l], coef_imag[l], NULL);

  // Calculate the binomial term and its square root
  mpz_bin_uiui(int_binom, 2 * l, l);
  mpf_set_z(float_binom, int_binom);
  mpf_sqrt(raiz_binom, float_binom);

  // Multiply the multipole moments by the binomial root
  mpf_mul(coef_real[l], raiz_binom, al_real[l]);
  mpf_mul(coef_imag[l], raiz_binom, al_imag[l]);

  for (int m = 1; m <= l; m++) {
    index = ((m * ((2 * LMAX) + 1 - m)) / 2) + l;

    // Initialize coefficient arrays for current m
    mpf_inits(coef_real[l - m], coef_imag[l - m], coef_real[l + m],
              coef_imag[l + m], NULL);

    // Calculate the binomial term and its square root
    mpz_bin_uiui(int_binom, 2 * l, m + l);
    mpf_set_z(float_binom, int_binom);
    mpf_sqrt(raiz_binom, float_binom);

    // Calculate real part for l - m
    mpf_mul(coef_real[l - m], raiz_binom, al_real[index]);
    mpf_set_si(coef_negativ, pow(-1, m));
    mpf_mul(coef_real[l - m], coef_real[l - m], coef_negativ);

    // Calculate imaginary part for l - m
    mpf_mul(coef_imag[l - m], raiz_binom, al_imag[index]);
    mpf_set_si(coef_negativ, pow(-1, m) * (-1));
    mpf_mul(coef_imag[l - m], coef_imag[l - m], coef_negativ);

    // Calculate real and imaginary parts for l + m
    mpf_mul(coef_real[l + m], raiz_binom, al_real[index]);
    mpf_mul(coef_imag[l + m], raiz_binom, al_imag[index]);
  }

  // Clear GMP variables
  mpz_clear(int_binom);
  mpf_clears(float_binom, raiz_binom, coef_negativ, NULL);
}

/*Function to extract the roots*/
void raizes_pol(int l, mpf_t coef_real[], mpf_t coef_imag[], double *raiz_real,
                double *raiz_imag) {
  // Define rational coefficients arrays
  mpq_t rat_coef_real[(2 * l) + 1], rat_coef_imag[(2 * l) + 1];

  // Convert floating-point coefficients to rational coefficients
  for (int i = 0; i < (2 * l) + 1; i++) {
    mpq_inits(rat_coef_real[i], rat_coef_imag[i], NULL);
    mpq_set_f(rat_coef_real[i], coef_real[i]);
    mpq_set_f(rat_coef_imag[i], coef_imag[i]);
  }

  // Initialize MPSolve context and polynomial
  mps_context *s = mps_context_new();
  if (s == NULL) {
    // Handle error
    return;
  }

  mps_monomial_poly *p = mps_monomial_poly_new(s, 2 * l);
  if (p == NULL) {
    // Handle error
    mps_context_free(s);
    return;
  }
  mps_context_select_algorithm(s, MPS_ALGORITHM_SECULAR_GA);

  // Set polynomial coefficients
  for (int i = 0; i < (2 * l) + 1; i++) {
    mps_monomial_poly_set_coefficient_q(s, p, i, rat_coef_real[i],
                                        rat_coef_imag[i]);
  }

  // Set input polynomial in the context
  mps_context_set_input_poly(s, MPS_POLYNOMIAL(p));

  // Allocate memory for results
  cplx_t *results = cplx_valloc(2 * l);
  if (results == NULL) {
    // Handle error
    mps_polynomial_free(s, MPS_POLYNOMIAL(p));
    mps_context_free(s);
    return;
  }

  // Solve the polynomial
  mps_mpsolve(s);
  mps_context_get_roots_d(s, &results, NULL);

  // Extract real and imaginary parts of the roots
  for (int i = 0; i < 2 * l; i++) {
    raiz_real[i] = cplx_Re(results[i]);
    raiz_imag[i] = cplx_Im(results[i]);
  }

  // Clear rational coefficients
  for (int i = 0; i < (2 * l) + 1; i++) {
    mpq_clears(rat_coef_real[i], rat_coef_imag[i], NULL);
  }

  // Free allocated memory
  mps_polynomial_free(s, MPS_POLYNOMIAL(p));
  mps_context_free(s);
  cplx_vfree(results);
}

/*Function to find the coordinates of multipole vectors*/
void coord_pol(int l, double raiz_real[], double raiz_imag[], double *theta,
               double *phi) {
  // Check for NULL pointers
  if (raiz_real == NULL || raiz_imag == NULL || theta == NULL || phi == NULL) {
    return;
  }

  // Define constants
  const int num_roots = 2 * l;

  // Allocate arrays for magnitudes and complex numbers
  double R[num_roots];
  double complex z[num_roots];

  // Calculate polar coordinates
  for (int i = 0; i < num_roots; i++) {
    z[i] = raiz_real[i] + (raiz_imag[i] * I);
    R[i] = cabs(z[i]);
    phi[i] = carg(z[i]);           // Phi domain is (0, 2*pi)
    theta[i] = 2 * atan(1 / R[i]); // Theta domain is (0, pi)
  }
}

/*Function to calculate the multipole vectors*/
void multipol_vec(double *al_real, double *al_imag, double *coef_theta,
                  double *coef_phi, int LMAX) {

  // Define constants
  //    const int MVS_NUMERO = (LMAX * (LMAX + 1)) - 2;
  int ALMS_NUMERO_LINHAS = ((LMAX + 1) * (LMAX + 2)) / 2;

  // Allocate memory for mpf_t arrays
  mpf_t *al_real_mpf = (mpf_t *)malloc(ALMS_NUMERO_LINHAS * sizeof(mpf_t));
  mpf_t *al_imag_mpf = (mpf_t *)malloc(ALMS_NUMERO_LINHAS * sizeof(mpf_t));
  if (al_real_mpf == NULL || al_imag_mpf == NULL) {
    fprintf(stderr, "Memory allocation failed\n");
    free(al_real_mpf);
    free(al_imag_mpf);
    return;
  }

  // Initialize mpf_t arrays and set values from double pointers
  for (int i = 0; i < ALMS_NUMERO_LINHAS; i++) {
    mpf_init(al_real_mpf[i]);
    mpf_set_d(al_real_mpf[i], al_real[i]);
    mpf_init(al_imag_mpf[i]);
    mpf_set_d(al_imag_mpf[i], al_imag[i]);
  }

  // Parallel loop for calculations
  #pragma omp parallel for ordered
  for (int l = 2; l <= LMAX; l++) {
      mpf_t coef_real[(2 * l) + 1], coef_imag[(2 * l) + 1];
      double raiz_real[2 * l], raiz_imag[2 * l];
      double theta[2 * l], phi[2 * l];

      // Calculate polynomial coefficients
      coefi_pol(l, al_real_mpf, al_imag_mpf, coef_real, coef_imag, LMAX);

      // Find polynomial roots
      raizes_pol(l, coef_real, coef_imag, raiz_real, raiz_imag);

      // Clear polynomial coefficients
      for (int i = 0; i < ((2 * l) + 1); i++) {
      mpf_clears(coef_real[i], coef_imag[i], NULL);
      }

      // Convert roots to polar coordinates
      coord_pol(l, raiz_real, raiz_imag, theta, phi);

      // Store theta and phi values
      int k = 0;
      for (int j = 0; j < 2 * l; j++) {
      k = (((l - 1) * (l - 1)) + (l - 1) - 2) + j;
      coef_theta[k] = theta[j];
      coef_phi[k] = phi[j];
      }
  }

  // Free allocated memory
  free(al_real_mpf);
  free(al_imag_mpf);
}
