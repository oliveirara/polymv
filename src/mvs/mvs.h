#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <gmp.h>
#include <mps/mps.h>
#include <omp.h>

void coefi_pol(int l, mpf_t al_real[], mpf_t al_imag[], mpf_t coef_real[],
               mpf_t coef_imag[], int LMAX);
void raizes_pol(int l, mpf_t coef_real[], mpf_t coef_imag[], double *raiz_real,
                double *raiz_imag);
void coord_pol(int l, double raiz_real[], double raiz_imag[], double *theta,
               double *phi);
void multipol_vec(double *al_real, double *al_imag, double *coef_theta,
                  double *coef_phi, int LMAX);
