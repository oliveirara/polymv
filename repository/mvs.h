#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <gmp.h>
#include <mps/mps.h>
#include <omp.h>

#define NSIDE 64
#define BUFFER_SIZE 8192
#define FILENAME_SIZE 400
#define ACOS(ab) acos (((ab) > 1.0) ? 1.0 : (((ab) < -1.0) ? -1.0 : (ab)))
#define USAGE_MESSAGE "Usage: %s <mc> <LMAX> <filename>\n"


void coefi_pol (int l, mpf_t al_real[], mpf_t al_imag[], mpf_t coef_real[], mpf_t coef_imag[], int LMAX);
void raizes_pol(int l, mpf_t coef_real[], mpf_t coef_imag[], double *raiz_real, double *raiz_imag);
void coord_pol(int l, double raiz_real[], double raiz_imag[], double *theta, double *phi);

void multipol_vec(double* al_real, double* al_imag, double* coef_theta, double* coef_phi, int LMAX);

