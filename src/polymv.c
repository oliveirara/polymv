#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <gmp.h>
#include <mps/mps.h>
#include <nlopt.h>
#include <chealpix.h>
#include <omp.h>
#include <hdf5.h>

#define NSIDE 64
#define BUFFER_SIZE 8192
#define FILENAME_SIZE 400
#define ACOS(ab) acos (((ab) > 1.0) ? 1.0 : (((ab) < -1.0) ? -1.0 : (ab)))
#define USAGE_MESSAGE "Usage: %s <LMAX> <input_filename> [output_filename]\n"

typedef struct _FrechetData
{
  int l;
  double * restrict x;
  double * restrict y;
  double * restrict z;
} FrechetData;

void print_progress_bar(int current, int total) {
    int bar_width = 50;
    float progress = (float)current / total;
    int pos = bar_width * progress;

    printf("[");
    for (int i = 0; i < bar_width; ++i) {
        if (i < pos) printf("=");
        else if (i == pos) printf(">");
        else printf(" ");
    }
    printf("] %d%%\r", (int)(progress * 100));
    fflush(stdout);
}

/* Function to calculate the Polynomial Coefficients*/
void coefi_pol (int l, mpf_t al_real[], mpf_t al_imag[], mpf_t coef_real[], mpf_t coef_imag[], int LMAX)
{
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

    for (int m = 1; m <= l; m++)
    {
        index = ((m * ((2 * LMAX) + 1 - m)) / 2) + l;

        // Initialize coefficient arrays for current m
        mpf_inits(coef_real[l - m], coef_imag[l - m], coef_real[l + m], coef_imag[l + m], NULL);

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

/* Function to extract the roots */
void raizes_pol(int l, mpf_t coef_real[], mpf_t coef_imag[], double *raiz_real, double *raiz_imag) {
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
    mps_context_select_algorithm(s, MPS_OUTPUT_GOAL_APPROXIMATE); // Original was MPS_ALGORITHM_SECULAR_GA
    mps_context_set_output_prec (s, 53);

    // Set polynomial coefficients
    for (int i = 0; i < (2 * l) + 1; i++) {
        mps_monomial_poly_set_coefficient_q(s, p, i, rat_coef_real[i], rat_coef_imag[i]);
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

/* Function to find the coordinates of multipole vectors */
void coord_pol(int l, double raiz_real[], double raiz_imag[], double *theta, double *phi) {
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
        phi[i] = carg(z[i]);         // Phi domain is (0, 2*pi)
        theta[i] = 2 * atan(1 / R[i]); // Theta domain is (0, pi)
    }
}

/*This function calculates the psi function of MVs for each ipix and returns the initial guess for minimize */
void guess(int ell, const double *x, const double *y, const double *z, double *s) {
    // Check for NULL pointers
    if (x == NULL || y == NULL || z == NULL || s == NULL) {
        return;
    }

    // Define constants
    const int npix = (12 * NSIDE * NSIDE) / 2; // Total number of pixels divided by 2

    // Allocate arrays for pixel coordinates
    double pixel_coords[npix][3];

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
            double acos_val = acos(dot_product);
            sum_arccos_squared += (acos_val * acos_val) + ((M_PI - acos_val) * (M_PI - acos_val));
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
}

double frechet_pol_min(unsigned n, const double *x, double *grad, void *my_func_data) {
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
        const double acos_ab = acos(ab);
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
            const double acos_ab = acos(ab);
            const double sqrt_1mab = sqrt(1.0 - ab * ab);

            grad[0] -= 2.0 * acos_ab * (dxdtheta * fd->x[i] + dydtheta * fd->y[i] + dzdtheta * fd->z[i]) / sqrt_1mab;
            grad[1] -= 2.0 * acos_ab * (dxdphi * fd->x[i] + dydphi * fd->y[i]) / sqrt_1mab;
        }
    }

    return frechet_mean;
}

void frechet_pol2(int l, double *restrict theta, double *restrict phi, double *frechet_vec_theta, double *frechet_vec_phi, double *min_f) {
    // Check for NULL pointers
    if (theta == NULL || phi == NULL || frechet_vec_theta == NULL || frechet_vec_phi == NULL || min_f == NULL) {
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
                printf("Antipodes Number not correspond\n");
                return;
            }
        }
    }

    if (muda != l) {
        printf("Antipodes Number not correspond\n");
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

        nlopt_opt opt = nlopt_create(NLOPT_LN_NELDERMEAD, 2);

        nlopt_set_lower_bounds(opt, lb);
        nlopt_set_upper_bounds(opt, ub);
        nlopt_set_min_objective(opt, &frechet_pol_min, &fd);
        nlopt_set_xtol_rel(opt, 1.0e-7);

        guess(l, x_ant, y_ant, z_ant, s);

        if (nlopt_optimize(opt, s, &f) < 0) {
            printf("nlopt failed!\n");
        } else {
            *min_f = f;
            *frechet_vec_theta = s[0];
            *frechet_vec_phi = s[1];
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

void multipol_vec(mpf_t al_real[], mpf_t al_imag[], int LMAX, const char *output_filename) {
    // Define constants
    const int MVS_NUMERO = (LMAX * (LMAX + 1)) - 2;

    // Allocate memory for mvs_theta and mvs_phi
    double *mvs_theta = (double *)malloc(MVS_NUMERO * sizeof(double));
    double *mvs_phi = (double *)malloc(MVS_NUMERO * sizeof(double));
    if (mvs_theta == NULL || mvs_phi == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        free(mvs_theta);
        free(mvs_phi);
        return;
    }

    // Allocate memory for frechet_vec_theta, frechet_vec_phi, and min_f
    double *frechet_vec_theta = (double *)malloc((LMAX - 1) * sizeof(double));
    double *frechet_vec_phi = (double *)malloc((LMAX - 1) * sizeof(double));
    double *min_f = (double *)malloc((LMAX - 1) * sizeof(double));
    if (frechet_vec_theta == NULL || frechet_vec_phi == NULL || min_f == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        free(mvs_theta);
        free(mvs_phi);
        free(frechet_vec_theta);
        free(frechet_vec_phi);
        free(min_f);
        return;
    }

    // Parallel loop for calculations
    #pragma omp parallel for ordered
    for (int l = 2; l <= LMAX; l++) {
        mpf_t coef_real[(2 * l) + 1], coef_imag[(2 * l) + 1];
        double raiz_real[2 * l], raiz_imag[2 * l];
        double theta[2 * l], phi[2 * l];
        double frechet_theta, frechet_phi, min_f_val;

        // Calculate polynomial coefficients
        coefi_pol(l, al_real, al_imag, coef_real, coef_imag, LMAX);

        // Find polynomial roots
        raizes_pol(l, coef_real, coef_imag, raiz_real, raiz_imag);

        // Clear polynomial coefficients
        for (int i = 0; i < ((2 * l) + 1); i++) {
            mpf_clears(coef_real[i], coef_imag[i], NULL);
        }

        // Convert roots to polar coordinates
        coord_pol(l, raiz_real, raiz_imag, theta, phi);

        // Calculate Frechet mean
        frechet_pol2(l, theta, phi, &frechet_theta, &frechet_phi, &min_f_val);

        // Store theta and phi values
        int k = 0;
        for (int j = 0; j < 2 * l; j++) {
            k = (pow((l - 1), 2) + (l - 1) - 2) + j;
            mvs_theta[k] = theta[j];
            mvs_phi[k] = phi[j];
        }

        // Store Frechet mean values and min_f
        frechet_vec_theta[l - 2] = frechet_theta;
        frechet_vec_phi[l - 2] = frechet_phi;
        min_f[l - 2] = min_f_val;

        #pragma omp critical
        {
            print_progress_bar(l - 1, LMAX - 1);
        }
    }

    // Create HDF5 file
    hid_t file_id = H5Fcreate(output_filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file_id < 0) {
        fprintf(stderr, "Failed to create HDF5 file\n");
        free(mvs_theta);
        free(mvs_phi);
        free(frechet_vec_theta);
        free(frechet_vec_phi);
        free(min_f);
        return;
    }

    // Create datasets for mvs_theta and mvs_phi
    hsize_t dims[1] = {MVS_NUMERO};
    hid_t dataspace_id = H5Screate_simple(1, dims, NULL);
    hid_t dataset_id = H5Dcreate(file_id, "mvs_theta", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, mvs_theta);
    H5Dclose(dataset_id);

    dataset_id = H5Dcreate(file_id, "mvs_phi", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, mvs_phi);
    H5Dclose(dataset_id);

    // Create datasets for frechet_vec_theta, frechet_vec_phi, and min_f
    dims[0] = LMAX - 1;
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate(file_id, "frechet_vec_theta", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, frechet_vec_theta);
    H5Dclose(dataset_id);

    dataset_id = H5Dcreate(file_id, "frechet_vec_phi", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, frechet_vec_phi);
    H5Dclose(dataset_id);

    dataset_id = H5Dcreate(file_id, "min_f", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, min_f);
    H5Dclose(dataset_id);

    // Close HDF5 file
    H5Fclose(file_id);

    // Free memory
    free(mvs_theta);
    free(mvs_phi);
    free(frechet_vec_theta);
    free(frechet_vec_phi);
    free(min_f);

    printf("\n"); // Move to the next line after the progress bar is complete
}

int main(int argc, char *argv[]) {
    // Check for correct number of arguments
    if (argc < 3) {
        fprintf(stderr, USAGE_MESSAGE, argv[0]);
        return 1;
    }

    // Parse command-line arguments
    int LMAX = atoi(argv[1]);
    char *input_filename = argv[2];
    char *output_filename = (argc > 3) ? argv[3] : "results.h5";

    // Calculate the number of lines
    int ALMS_NUMERO_LINHAS = ((LMAX + 1) * (LMAX + 2)) / 2;

    // Allocate memory for real and imaginary parts
    mpf_t *al_real = (mpf_t *)malloc(ALMS_NUMERO_LINHAS * sizeof(mpf_t));
    mpf_t *al_imag = (mpf_t *)malloc(ALMS_NUMERO_LINHAS * sizeof(mpf_t));
    if (al_real == NULL || al_imag == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        free(al_real);
        free(al_imag);
        return 1;
    }

    // Open the file for reading
    FILE *file = fopen(input_filename, "r");
    if (!file) {
        perror("Error opening file");
        free(al_real);
        free(al_imag);
        return 1;
    }

    // Initialize and read values from the file
    for (int j = 0; j < ALMS_NUMERO_LINHAS; j++) {
        mpf_inits(al_real[j], al_imag[j], NULL);
        if (gmp_fscanf(file, "%Ff %Ff\n", &al_real[j], &al_imag[j]) != 2) {
            fprintf(stderr, "Error reading file\n");
            fclose(file);
            for (int k = 0; k <= j; k++) {
                mpf_clears(al_real[k], al_imag[k], NULL);
            }
            free(al_real);
            free(al_imag);
            return 1;
        }
    }

    // Close the file
    fclose(file);

    // Perform calculations
    multipol_vec(al_real, al_imag, LMAX, output_filename);

    // Clear and free memory
    for (int al_num = 0; al_num < ALMS_NUMERO_LINHAS; al_num++) {
        mpf_clears(al_real[al_num], al_imag[al_num], NULL);
    }
    free(al_real);
    free(al_imag);

    return 0;
}