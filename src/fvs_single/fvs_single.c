#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <nlopt.h>
#include <chealpix.h>
#include <omp.h>

#define M_PI 3.14159265358979323846
#define NSIDE 64
#define ACOS(ab) acos(((ab) > 1.0) ? 1.0 : (((ab) < -1.0) ? -1.0 : (ab)))

//This function might be overdesigned, I was trying to fix a problem where the code would instalock on the poles, but I think it might be fine now. I'm leaving it as is for now, but we can come back and simplify if we want to.
// The solution for the instalocking seems to be adding the antipodes BEFORE running this function, which is a negative side of having a generic function, it cannot predict if the input is all in the north pole

typedef struct _SingleData {
  int n;
  double *restrict x;
  double *restrict y;
  double *restrict z;
} SingleData;

/*This function calculates the psi function of MVs for each ipix and returns the
 * initial guess for minimize, copied straight from fvs.c */
void guess_single(int ell, const double *x, const double *y, const double *z,
           double *s) {
  // Check for NULL pointers
  if (x == NULL || y == NULL || z == NULL || s == NULL) {
    return;
  }
  // We'll do 2 passes just like the original python code did, this is to avoid having the code instalock on the poles
  // Define constants
  const int npix_coarse =
      (12 * NSIDE * NSIDE ) ; // Full sky search with NSIDE=64, no antipode pairing, so we use all pixels

  // Use heap allocation for large arrays
  double (*pixel_coords)[3] = malloc(sizeof(*pixel_coords) * npix_coarse);
  if (pixel_coords == NULL) {
    fprintf(stderr, "Memory allocation failed for pixel_coords\n");
    return;
  }

  // Initialize variables for the minimum psi and the guess coordinates
  double psi_min = 1.0e300;
  double guess[3];

  // Calculate pixel coordinates
  for (int ipix = 0; ipix < npix_coarse; ipix++) {
    double vec[3];
    pix2vec_ring(NSIDE, ipix, vec);
    pixel_coords[ipix][0] = vec[0];
    pixel_coords[ipix][1] = vec[1];
    pixel_coords[ipix][2] = vec[2];
  }

  // Iterate over each pixel to find the minimum psi
  for (int ipix = 0; ipix < npix_coarse; ipix++) {
    double sum_arccos_squared = 0.0;

    // Calculate the sum of squared arccosines
    for (int j = 0; j < ell; j++) {
      double dot_product = (pixel_coords[ipix][0] * x[j]) +
                           (pixel_coords[ipix][1] * y[j]) +
                           (pixel_coords[ipix][2] * z[j]);

      // Clamp to avoid domain error in acos
      double dp_clamped = fmax(-1.0, fmin(1.0, dot_product));
      double acos_val = acos(dp_clamped);

      sum_arccos_squared += (acos_val * acos_val);
    }

    // Update the minimum psi and guess coordinates if a new minimum is found
    if (sum_arccos_squared < psi_min) {
      psi_min = sum_arccos_squared;
      guess[0] = pixel_coords[ipix][0];
      guess[1] = pixel_coords[ipix][1];
      guess[2] = pixel_coords[ipix][2];
    }
  }
  // Second pass, finer grain search around the first guess
  const long NSIDE_FINE = 1024;
  const long npix_fine = 12 * NSIDE_FINE * NSIDE_FINE ; // Full sky search with NSIDE=1024, no antipode pairing, so we use all pixels
  const double search_radius = cos(2 * M_PI / 180.0); // 2 degree search radius
  psi_min = 1.0e300; // reset psi_min for the finer search
  for (long ipix = 0; ipix < npix_fine; ipix++) {
    double vec[3];
    pix2vec_ring(NSIDE_FINE, ipix, vec);

    // Check if the pixel is within the search radius of the coarse guess
    double dot_product = (vec[0] * guess[0]) + (vec[1] * guess[1]) + (vec[2] * guess[2]);
    if (dot_product > search_radius) {
      double sum_arccos_squared = 0.0;
      for (int j = 0; j < ell; j++) {
        double dp_clamped = fmax(-1.0, fmin(1.0, (vec[0] * x[j]) + (vec[1] * y[j]) + (vec[2] * z[j])));
        double acos_val = acos(dp_clamped);
        sum_arccos_squared += (acos_val * acos_val);
      }
      if (sum_arccos_squared < psi_min) {
        psi_min = sum_arccos_squared;
        guess[0] = vec[0];
        guess[1] = vec[1];
        guess[2] = vec[2];
      }
    }
  }
  // Convert the guess coordinates to spherical coordinates
  double theta_frechet[1], phi_frechet[1];
  vec2ang(guess, theta_frechet, phi_frechet);

  // Store the result in the output array
  s[0] = theta_frechet[0];
  s[1] = phi_frechet[0];

  free(pixel_coords); // I'm leaving this at the end of the whole function in case we come back and change the second iteration to include this variable. If we do so, this will free the memory anyways.
}


double frechet_pol_min_single(unsigned dim, const double *x, double *grad,
                       void *data) {
  // Check for NULL pointers
  if (x == NULL || data == NULL) {
    return HUGE_VAL; // Return an error value
  }

  SingleData *sd = (SingleData *)data;
  // Extract theta and phi from input array
  const double theta = x[0];
  const double phi = x[1];
  double frechet_mean = 0.0;

  // Calculate Cartesian coordinates
  const double x_v = sin(theta) * cos(phi);
  const double y_v = sin(theta) * sin(phi);
  const double z_v = cos(theta);

  // Calculate Frechet mean
  for (int i = 0; i < sd->n; i++) {
    const double ab = (x_v * sd->x[i]) + (y_v * sd->y[i]) + (z_v * sd->z[i]);
    const double ab_clamped = fmax(-1.0, fmin(1.0, ab));
    const double acos_ab = acos(ab_clamped);
    frechet_mean += acos_ab * acos_ab;
  }

  // Calculate gradient if required
  if (grad != NULL) {
    const double dxdtheta = cos(theta) * cos(phi);
    const double dydtheta = cos(theta) * sin(phi);
    const double dzdtheta = -sin(theta);

    const double dxdphi = -sin(theta) * sin(phi);
    const double dydphi = sin(theta) * cos(phi);

    grad[0] = 0.0;
    grad[1] = 0.0;

    for (int i = 0; i < sd->n; i++) {
      const double ab = (x_v * sd->x[i]) + (y_v * sd->y[i]) + (z_v * sd->z[i]);
      const double ab_clamped = fmax(-1.0, fmin(1.0, ab));
      const double acos_ab = acos(ab_clamped);
      const double sqrt_1mab = sqrt(1.0 - ab_clamped * ab_clamped);

      if (fabs(ab_clamped) == 1.0) {
      // Gradient is ill-defined here; handle carefully or skip
        continue;
      }

      grad[0] -=
          2.0 * acos_ab *
          (dxdtheta * sd->x[i] + dydtheta * sd->y[i] + dzdtheta * sd->z[i]) /
          sqrt_1mab;
      grad[1] -=
          2.0 * acos_ab * (dxdphi * sd->x[i] + dydphi * sd->y[i]) / sqrt_1mab;
    }
  }

  return frechet_mean;
}

/*
 * fvs_single - Computes one Fréchet vector for an arbitrary set of
 * polar unit vectors on S2. No antipode pairing is assumed.
 *
 * Parameters:
 *   theta    - input polar angles, length n_vecs  (radians, [0, pi])
 *   phi      - input azimuthal angles, length n_vecs (radians, [0, 2pi])
 *   n_vecs   - number of input vectors (any positive integer)
 *   out_theta - output polar angle of the Fréchet mean
 *   out_phi   - output azimuthal angle of the Fréchet mean
 */
void fvs_single(double *theta, double *phi, int n_vecs, double *out_theta, double *out_phi) {

  if (theta == NULL || phi == NULL || out_theta == NULL || out_phi == NULL) {
    return;
  }
  *out_theta = NAN;
  *out_phi = NAN;
  if (n_vecs <= 0) {
    fprintf(stderr, "fvs_single: n_vecs must be positive (got %d)\n", n_vecs);
    return;
  }

  double *x = malloc(sizeof(double) * n_vecs);
  double *y = malloc(sizeof(double) * n_vecs);
  double *z = malloc(sizeof(double) * n_vecs);
  if (x == NULL || y == NULL || z == NULL) {
    free(x); free(y); free(z);
    fprintf(stderr, "fvs_single: memory allocation failed\n");
    return;
  }
  for (int i = 0; i < n_vecs; i++) {
    x[i] = sin(theta[i]) * cos(phi[i]);
    y[i] = sin(theta[i]) * sin(phi[i]);
    z[i] = cos(theta[i]);
  }

  SingleData sd = {n_vecs, x, y, z};

  double s[2] = {0.0, 0.0};  // will be overwritten by guess_single
  double lb[2] = {0.0, 0.0};
  double ub[2] = {M_PI, 2.0 * M_PI};
  double f;

  // STEP 1: get initial guess from pixel scan (must happen first)
  guess_single(n_vecs, x, y, z, s);
  // fprintf(stderr, "DEBUG guess: theta=%.4f (%.1f°) phi=%.4f (%.1f°)\n", s[0], s[0]*180/M_PI, s[1], s[1]*180/M_PI);

  // STEP 2: coarse Nelder-Mead from the guess
  nlopt_opt opt1 = nlopt_create(NLOPT_LN_NELDERMEAD, 2);
  if (opt1 == NULL ||
    nlopt_set_lower_bounds(opt1, lb) < 0 ||
    nlopt_set_upper_bounds(opt1, ub) < 0 ||
    nlopt_set_min_objective(opt1, frechet_pol_min_single, &sd) < 0 ||
    nlopt_set_xtol_rel(opt1, 1e-4) < 0 ||
    nlopt_optimize(opt1, s, &f) < 0) {
    if (opt1 != NULL) {
      nlopt_destroy(opt1);
    }
    fprintf(stderr, "fvs_single: coarse optimization failed\n");
    free(x); free(y); free(z);
    return;
  }
  nlopt_destroy(opt1);

  // Try guess and its antipode, keep the better result
  double s_antipode[2] = {
    M_PI - s[0],                                    // flip theta
    s[1] > M_PI ? s[1] - M_PI : s[1] + M_PI        // flip phi
  };

  double f1 = 1e300, f2 = 1e300;
  double s1[2] = {s[0], s[1]};
  double s2[2] = {s_antipode[0], s_antipode[1]};

  // Optimize from original guess
  nlopt_opt opt2a = nlopt_create(NLOPT_LD_LBFGS, 2);
  nlopt_set_lower_bounds(opt2a, lb);
  nlopt_set_upper_bounds(opt2a, ub);
  nlopt_set_min_objective(opt2a, frechet_pol_min_single, &sd);
  nlopt_set_xtol_rel(opt2a, 1e-7);
  nlopt_optimize(opt2a, s1, &f1);
  nlopt_destroy(opt2a);

  // Optimize from antipodal guess
  nlopt_opt opt2b = nlopt_create(NLOPT_LD_LBFGS, 2);
  nlopt_set_lower_bounds(opt2b, lb);
  nlopt_set_upper_bounds(opt2b, ub);
  nlopt_set_min_objective(opt2b, frechet_pol_min_single, &sd);
  nlopt_set_xtol_rel(opt2b, 1e-7);
  nlopt_optimize(opt2b, s2, &f2);
  nlopt_destroy(opt2b);

  // Pick the better result
  if (f1 <= f2) {
    *out_theta = s1[0];
    *out_phi   = s1[1];
  } else {
    *out_theta = s2[0];
    *out_phi   = s2[1];
  }

  free(x); free(y); free(z);
}
