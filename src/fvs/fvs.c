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

typedef struct _FrechetData {
  int l;
  double *restrict x;
  double *restrict y;
  double *restrict z;
} FrechetData;

/*This function calculates the psi function of MVs for each ipix and returns the
 * initial guess for minimize */
void guess(int ell, const double *x, const double *y, const double *z,
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

      sum_arccos_squared +=
          (acos_val * acos_val) + ((M_PI - acos_val) * (M_PI - acos_val));
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

double frechet_pol_min(unsigned n, const double *x, double *grad,
                       void *my_func_data) {
  // Check for NULL pointers
  if (x == NULL || my_func_data == NULL) {
    return -1.0; // Return an error value
  }

  FrechetData *fd = (FrechetData *)my_func_data;
  const int twol = 2 * fd->l;
  double frechet_mean = 0.0;
  double x_v, y_v, z_v;

  // Extract theta and phi from input array
  const double theta = x[0];
  const double phi = x[1];

  // Calculate Cartesian coordinates
  x_v = sin(theta) * cos(phi);
  y_v = sin(theta) * sin(phi);
  z_v = cos(theta);

  // Calculate Frechet mean
  for (int i = 0; i < twol; i++) {
    const double ab = (x_v * fd->x[i]) + (y_v * fd->y[i]) + (z_v * fd->z[i]);
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

    for (int i = 0; i < twol; i++) {
      const double ab = (x_v * fd->x[i]) + (y_v * fd->y[i]) + (z_v * fd->z[i]);
      const double ab_clamped = fmax(-1.0, fmin(1.0, ab));
      const double acos_ab = acos(ab_clamped);
      const double sqrt_1mab = sqrt(1.0 - ab_clamped * ab_clamped);
+     if (fabs(ab_clamped) == 1.0) {
+       // Gradient is ill-defined here; handle carefully or skip
+       continue;
+     }
      grad[0] -=
          2.0 * acos_ab *
          (dxdtheta * fd->x[i] + dydtheta * fd->y[i] + dzdtheta * fd->z[i]) /
          sqrt_1mab;
      grad[1] -=
          2.0 * acos_ab * (dxdphi * fd->x[i] + dydphi * fd->y[i]) / sqrt_1mab;
    }
  }

  return frechet_mean;
}
void frechet_pol(int l, double *restrict theta, double *restrict phi,
                 double *frechet_vec_theta, double *frechet_vec_phi) {
  // Check for NULL pointers
  if (theta == NULL || phi == NULL || frechet_vec_theta == NULL ||
      frechet_vec_phi == NULL) {
    return;
  }

  const int twol = 2 * l;
  double x[twol], y[twol], z[twol];
  double x_ant[l], y_ant[l], z_ant[l];
  double f;
  int i;
  FrechetData fd = {l, x, y, z};

  *frechet_vec_theta = 0.0;
  *frechet_vec_phi = 0.0;

  double ant_theta[l], ant_phi[l];

  int muda = 0;
  for (i = 0; i < twol; i++) {
    if (theta[i] < M_PI / 2) {
      if (muda < l) {
        ant_theta[muda] = theta[i];
        ant_phi[muda] = phi[i];
        muda++;
      } else {
        printf("First Advisor - Antipodes Number not correspond\n");
        return;
      }
    }
  }

  if (muda != l) {
    printf("Second Advisor - Antipodes Number not correspond\n");
    return;
  }

  for (i = 0; i < l; i++) {
    x_ant[i] = sin(ant_theta[i]) * cos(ant_phi[i]);
    y_ant[i] = sin(ant_theta[i]) * sin(ant_phi[i]);
    z_ant[i] = cos(ant_theta[i]);
  }

  for (i = 0; i < twol; i++) {
    x[i] = sin(theta[i]) * cos(phi[i]);
    y[i] = sin(theta[i]) * sin(phi[i]);
    z[i] = cos(theta[i]);
  }

  {
    double lb[2] = {0, 0};
    double ub[2] = {M_PI, 2.0 * M_PI};
    double s[2] = {0.0, 0.0};
    double min_f = 1.0e300;

    nlopt_opt opt = nlopt_create(NLOPT_LN_NELDERMEAD, 2);

    nlopt_set_lower_bounds(opt, lb);
    nlopt_set_upper_bounds(opt, ub);
    nlopt_set_min_objective(opt, &frechet_pol_min, &fd);
    nlopt_set_xtol_rel(opt, 1.0e-7);

    guess(l, x_ant, y_ant, z_ant, s);

    if (nlopt_optimize(opt, s, &f) < 0) {
      printf("nlopt failed!\n");
    } else {
      if (f < min_f) {
        min_f = f;
        *frechet_vec_theta = s[0];
        *frechet_vec_phi = s[1];
      }
    }

    if (*frechet_vec_theta > (M_PI / 2)) {
      *frechet_vec_theta = M_PI - *frechet_vec_theta;
      *frechet_vec_phi = M_PI + *frechet_vec_phi;

      if (*frechet_vec_phi > (2 * M_PI)) {
        *frechet_vec_phi -= (2 * M_PI);
      }
    }

    nlopt_destroy(opt);
  }
}

void frechet_vec(double *coef_theta, double *coef_phi, double *frechet_theta,
                 double *frechet_phi, int LMAX) {
// Parallel loop for calculations
#pragma omp parallel for ordered
  for (int l = 2; l <= LMAX; l++) {

    double frechet_vec_theta, frechet_vec_phi;

    int lim_min = ((pow((l - 1), 2)) + (l - 1) - 2);
    int lim_max = ((pow(l, 2)) + l - 2);

    double theta_each_l[2 * l], phi_each_l[2 * l];
    int n = 0;

    for (int pos = lim_min; pos < lim_max; pos++) {
      theta_each_l[n] = coef_theta[pos];
      phi_each_l[n] = coef_phi[pos];

      n++;
    }

    // Calculate the Frechet mean
    frechet_pol(l, theta_each_l, phi_each_l, &frechet_vec_theta,
                &frechet_vec_phi);

    frechet_theta[l - 2] = frechet_vec_theta;
    frechet_phi[l - 2] = frechet_vec_phi;
  }
}
