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

  // Define constants
  const int npix =
      (12 * NSIDE * NSIDE) / 2; // Total number of pixels divided by 2

  // Use heap allocation for large arrays
  double (*pixel_coords)[3] = malloc(sizeof(*pixel_coords) * npix);
  if (pixel_coords == NULL) {
    fprintf(stderr, "Memory allocation failed for pixel_coords\n");
    return;
  }

  // Initialize variables for the minimum psi and the guess coordinates
  double psi_min = 1.0e300;
  double guess[3];

  // Calculate pixel coordinates
  for (int ipix = 0; ipix < npix; ipix++) {
    double vec[3];
    pix2vec_ring(NSIDE, ipix, vec);
    pixel_coords[ipix][0] = vec[0];
    pixel_coords[ipix][1] = vec[1];
    pixel_coords[ipix][2] = vec[2];
  }

  // Iterate over each pixel to find the minimum psi
  for (int ipix = 0; ipix < npix; ipix++) {
    double sum_arccos_squared = 0.0;

    // Calculate the sum of squared arccosines
    for (int pos_mv = 0; pos_mv < ell; pos_mv++) {
      double dot_product = (pixel_coords[ipix][0] * x[pos_mv]) +
                           (pixel_coords[ipix][1] * y[pos_mv]) +
                           (pixel_coords[ipix][2] * z[pos_mv]);

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

  // Convert the guess coordinates to spherical coordinates
  double theta_frechet[1], phi_frechet[1];
  vec2ang(guess, theta_frechet, phi_frechet);

  // Store the result in the output array
  s[0] = theta_frechet[0];
  s[1] = phi_frechet[0];

  free(pixel_coords);
}


double frechet_pol_min_single(unsigned dim, const double *x, double *grad,
                       void *data) {
  // Check for NULL pointers
  if (x == NULL || data == NULL) {
    return -1.0; // Return an error value
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
                    // Check for NULL pointers and valid n_vecs
    if (theta == NULL || phi == NULL || out_theta == NULL || out_phi == NULL) {
        return;
    }
    if (n_vecs <= 0) {
        fprintf(stderr, "fvs_single: n_vecs must be positive (got %d)\n", n_vecs);
        return;
    }

    /* Convert to Cartesian */
    double *x = malloc(sizeof(double) * n_vecs);
    double *y = malloc(sizeof(double) * n_vecs);
    double *z = malloc(sizeof(double) * n_vecs);
    double f;

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

    double s[2] = {0.0, 0.0};
    double lb[2] = {0, 0};
    double ub[2] = {M_PI, 2.0 * M_PI };
    double min_f = 1.0e300;

    nlopt_opt opt = nlopt_create(NLOPT_LN_NELDERMEAD, 2);

    nlopt_set_lower_bounds(opt, lb);
    nlopt_set_upper_bounds(opt, ub);
    nlopt_set_min_objective(opt, frechet_pol_min_single, &sd);
    nlopt_set_xtol_rel(opt, 1.0e-7);

    guess_single(n_vecs, x, y, z, s);

    if (nlopt_optimize(opt, s, &f) < 0) {
        fprintf(stderr, "fvs_single: nlopt failed\n");
    } else {
        if (f < min_f) {
            min_f = f;
            *out_theta = s[0];
            *out_phi   = s[1];
        }
    }
    if (*out_theta > (M_PI / 2)) {
        *out_theta = M_PI - *out_theta;
        *out_phi   = M_PI + *out_phi;

        if (*out_phi > (2.0 * M_PI)) {
            *out_phi -= (2.0 * M_PI);
        }
    }

    nlopt_destroy(opt);
    free(x); free(y); free(z);
}
