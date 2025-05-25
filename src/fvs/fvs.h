#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <complex.h>
#include <chealpix.h>
#include <nlopt.h>
#include <omp.h>

typedef struct _FrechetData {
  int l;
  double *restrict x;
  double *restrict y;
  double *restrict z;
} FrechetData;

void guess(int ell, const double *x, const double *y, const double *z,
           double *s);
double frechet_pol_min(unsigned n, const double *x, double *grad,
                       void *my_func_data);
void frechet_pol(int l, double *restrict theta, double *restrict phi,
                 double *frechet_vec_theta, double *frechet_vec_phi);
void frechet_vec(double *coef_theta, double *coef_phi, double *frechet_theta,
                 double *frechet_phi, int LMAX);
