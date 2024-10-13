#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <gmp.h>
#include <mps/mps.h>
#include <nlopt.h>
#include <chealpix.h>

#define LMAX 2000
#define NSIDE 64
#define ALMS_NUMERO_LINHAS ((LMAX+1) * (LMAX+2))/2
#define MVS_NUMERO (LMAX * (LMAX+1)) - 2 
#define ACOS(ab) acos (((ab) > 1.0) ? 1.0 : (((ab) < -1.0) ? -1.0 : (ab)))

typedef struct _FrechetData
{
  int l;
  double * restrict x;
  double * restrict y;
  double * restrict z;
} FrechetData;

/* Function to calculate the Polynomial Coefficients*/
void
coefi_pol (int l, mpf_t al_real[], mpf_t al_imag[], mpf_t coef_real[], mpf_t coef_imag[])
{
  mpz_t int_binom;
  mpf_t float_binom;
  mpf_t raiz_binom;
  mpf_t coef_negativ;
  int index; // This number can explode to big lists

  mpz_inits (int_binom, NULL);
  mpf_inits (float_binom, raiz_binom, coef_negativ, NULL);

  mpf_inits (coef_real[l], coef_imag[l], NULL);

  mpz_bin_uiui (int_binom, (2 * l), l); // Calculating the binomial term
  mpf_set_z (float_binom, int_binom);   // Convert integer to float for extract the roots
  mpf_sqrt (raiz_binom, float_binom);   // Extract the roots

  mpf_mul (coef_real[l], raiz_binom, al_real[l]); // Multipling the multipol moments for the binomial root
  mpf_mul (coef_imag[l], raiz_binom, al_imag[l]);

  for (int m = 1; m <= l; m++)
  {
    index = ((m * ((2 * LMAX) + 1 - m)) / 2) + l; // Calculating the position of alms
    mpf_inits (coef_real[l - m], coef_imag[l - m], coef_real[l + m], coef_imag[l + m], NULL);

    mpz_bin_uiui (int_binom, (2 * l), (m + l));   // Calculating the binomial term
    mpf_set_z (float_binom, int_binom);           // Convert integer to float for extract the roots
    mpf_sqrt (raiz_binom, float_binom);           // Extract the roots

    mpf_mul (coef_real[l - m], raiz_binom, al_real[index]);
    mpf_set_si (coef_negativ, pow (-1, m));
    mpf_mul (coef_real[l - m], coef_real[l - m], coef_negativ);

    mpf_mul (coef_imag[l - m], raiz_binom, al_imag[index]);
    mpf_set_si (coef_negativ, pow (-1, m) * (-1));
    mpf_mul (coef_imag[l - m], coef_imag[l - m], coef_negativ);

    mpf_mul (coef_real[l + m], raiz_binom, al_real[index]);
    mpf_mul (coef_imag[l + m], raiz_binom, al_imag[index]);
  }

  mpz_clears (int_binom, NULL);
  mpf_clears (float_binom, raiz_binom, NULL);
}

/* Function to extract the roots */
void
raizes_pol (int l, mpf_t coef_real[], mpf_t coef_imag[], double *raiz_real, double *raiz_imag)
{
  mpq_t rat_coef_real[(2 * l) + 1], rat_coef_imag[(2 * l) + 1];

  for (int i = 0; i < (2 * l) + 1; i++)
  { /* Convert float in rational */
    mpq_inits (rat_coef_real[i], rat_coef_imag[i], NULL);
    mpq_set_f (rat_coef_real[i], coef_real[i]);
    mpq_set_f (rat_coef_imag[i], coef_imag[i]);
  }

  mps_monomial_poly *p;
  mps_context *s;

  s = mps_context_new ();
  p = mps_monomial_poly_new (s, 2 * l);
  mps_context_select_algorithm (s, MPS_ALGORITHM_SECULAR_GA);

  for (int i = 0; i < ((2 * l) + 1); i++)
  {
    mps_monomial_poly_set_coefficient_q (s, p, i, rat_coef_real[i], rat_coef_imag[i]);
  }

  mps_context_set_input_poly (s, MPS_POLYNOMIAL (p));

  cplx_t *results = cplx_valloc (2 * l);

  mps_mpsolve (s);
  mps_context_get_roots_d (s, &results, NULL);

  for (int i = 0; i < (2 * l); i++)
  {
    raiz_real[i] = cplx_Re (results[i]);
    raiz_imag[i] = cplx_Im (results[i]);
  }

  for (int i = 0; i < ((2 * l) + 1); i++)
  {
    mpq_clears (rat_coef_real[i], rat_coef_imag[i], NULL);
  }

  mps_polynomial_free (s, MPS_POLYNOMIAL (p));
  mps_context_free (s);
  cplx_vfree (results);
}

/* Function to find the coordinates of multipole vectors */
void
coord_pol (int l, double raiz_real[], double raiz_imag[], double *theta, double *phi)
{

  double R[2 * l];
  double complex z[2 * l];

  for (int i = 0; i < (2 * l); i++)
  {
    z[i]     = raiz_real[i] + (raiz_imag[i] * I);
    R[i]     = cabs (z[i]);
    phi[i]   = carg (z[i]);         // Phi domain is (0, 2*pi)
    theta[i] = 2 * atan (1 / R[i]); // Theta domain is (0, pi)
  }
}

/*This function calculates the psi function of MVs for each ipix and returns the initial guess for minimize */
void 
guess(int ell, const double *x, const double *y, const double *z, double *s)
{
    int npix = (12 * NSIDE * NSIDE) / 2; // Total number of pixels divide per 2

    double pixel_coords[npix][3];        // ipix that we work
    double vec[3];

    double psi_min = 1.0e300;
    double chute_x, chute_y, chute_z;
    double psi[npix];
    double guess[3];

    for (int ipix = 0; ipix < npix; ipix++) {
        double vec[3];
        pix2vec_ring(NSIDE, ipix, vec);
        pixel_coords[ipix][0] = vec[0];
        pixel_coords[ipix][1] = vec[1];
        pixel_coords[ipix][2] = vec[2];
    }

    for(int ipix = 0; ipix < npix; ipix ++){
        double sum_arccos_squared = 0.0;
        double dot_product = 0.0;

        for(int pos_mv = 0; pos_mv < (ell); pos_mv ++){
            dot_product =  (pixel_coords[ipix][0] * x[pos_mv]) + (pixel_coords[ipix][1] * y[pos_mv]) + (pixel_coords[ipix][2] * z[pos_mv]);
            sum_arccos_squared += (ACOS(dot_product) * ACOS(dot_product)) + ((M_PI - ACOS(dot_product)) * (M_PI - ACOS(dot_product)));
        }

        psi[ipix] = sum_arccos_squared;

        if(sum_arccos_squared < psi_min){
            psi_min = sum_arccos_squared;

            guess[0] = pixel_coords[ipix][0];
            guess[1] = pixel_coords[ipix][1];
            guess[2] = pixel_coords[ipix][2];
        }
    }

    double theta_frechet[1], phi_frechet[1];
    vec2ang(guess, theta_frechet, phi_frechet);

    s[0] = theta_frechet[0];
    s[1] = phi_frechet[0];
}

double 
frechet_pol_min (unsigned n, const double *x, double *grad, void *my_func_data)
{
  FrechetData *fd = (FrechetData *) my_func_data;
  const int twol = 2 * fd->l;
  double fechet_mean = 0.0;  
  double x_v, y_v, z_v;
  int i;

  {
    const double theta = x[0];
    const double phi   = x[1];

    x_v = sin (theta) * cos (phi);
    y_v = sin (theta) * sin (phi);
    z_v = cos (theta);

    for (i = 0; i < twol; i++)
    {
      const double ab      = (x_v * fd->x[i]) + (y_v * fd->y[i]) + (z_v * fd->z[i]);
      const double acos_ab = ACOS (ab);
      fechet_mean += acos_ab * acos_ab;
    }
    
    if (grad != NULL)
    {
      const double dxdtheta = cos (theta) * cos (phi);
      const double dydtheta = cos (theta) * sin (phi);
      const double dzdtheta = - sin (theta);
    
      const double dxdphi = - sin (theta) * sin (phi);
      const double dydphi = + sin (theta) * cos (phi);
    
      grad[0] = 0.0;
      grad[1] = 0.0;
    
      for (i = 0; i < twol; i++)
      {
        const double ab      = (x_v * fd->x[i]) + (y_v * fd->y[i]) + (z_v * fd->z[i]);
        const double acos_ab = ACOS (ab);
        const double sqrt_1mab = sqrt (1.0 - ab * ab);

        grad[0] -= 2.0 * acos_ab * (dxdtheta * fd->x[i] + dydtheta * fd->y[i] + dzdtheta * fd->z[i]) / sqrt_1mab;
        grad[1] -= 2.0 * acos_ab * (dxdphi   * fd->x[i] + dydphi   * fd->y[i]                      ) / sqrt_1mab;
      }
    }
  }

  return fechet_mean;
}

void
frechet_pol2 (int l, double * restrict theta, double * restrict phi, double *frechet_vec_theta, double *frechet_vec_phi)
{
  const int twol = 2 * l;
  double x[2 * l], y[2 * l], z[2 * l];
  double x_ant[2 * l], y_ant[2 * l], z_ant[2 * l];
  double f;
  int i;
  FrechetData fd = {l, x, y, z};

  *frechet_vec_theta = 0.0;
  *frechet_vec_phi = 0.0;

  double ant_theta[l], ant_phi[l];

  int muda = 0;
  for (i = 0; i < twol; i++)
  {
    if(theta[i] < M_PI/2)
    {
      ant_theta[muda] = theta[i];
      ant_phi[muda] = phi[i];
      muda++;
    }
  }

  if(muda != l){
    printf("Antipodes Number not correspond\n");
  }

  for (i = 0; i < l; i++)
  {
    x_ant[i] = sin (ant_theta[i]) * cos (ant_phi[i]);
    y_ant[i] = sin (ant_theta[i]) * sin (ant_phi[i]);
    z_ant[i] = cos (ant_theta[i]);
  }

  for (i = 0; i < twol; i++)
  {
    x[i] = sin (theta[i]) * cos (phi[i]);
    y[i] = sin (theta[i]) * sin (phi[i]);
    z[i] = cos (theta[i]);
  }

  {
    double lb[2] = {0,  0};
    double ub[2] = {M_PI,  + 2.0 * M_PI};
    double s[2]  = {0.0, 0.0};
    double min_f = 1.0e300;

    nlopt_opt opt;
    
    opt = nlopt_create (NLOPT_LN_NELDERMEAD, 2);

    nlopt_set_lower_bounds (opt, lb);
    nlopt_set_upper_bounds (opt, ub);
    nlopt_set_min_objective(opt, &frechet_pol_min, &fd);
    nlopt_set_xtol_rel(opt, 1.0e-7);

    guess(l, x_ant, y_ant, z_ant, s);

    if (nlopt_optimize(opt, s, &f) < 0) 
    {
      printf("nlopt failed!\n");
    }
    else 
    {
      if (f < min_f)
      {
        min_f = f;
        *frechet_vec_theta = s[0];
        *frechet_vec_phi   = s[1];
      }
    }

    if (*frechet_vec_theta > (M_PI / 2))
    {
    *frechet_vec_theta = M_PI - *frechet_vec_theta;
    *frechet_vec_phi   = M_PI + *frechet_vec_phi;

    if (*frechet_vec_phi > (2 * M_PI))
    {
        *frechet_vec_phi = *frechet_vec_phi - (2 * M_PI);
    }
    }

    nlopt_destroy (opt);
  }
}

void
multipol_vec (int i, mpf_t al_real[], mpf_t al_imag[])
{
  {
    static  double mvs_theta[MVS_NUMERO];
    static double mvs_phi[MVS_NUMERO];

    char filename[400];
    sprintf(filename,"mvs_fullsky_smica_noise_%d_.dat", i);

    FILE *FVs_theta_phi = fopen ("FVs_theta_phi.dat", "a");

//    char buff0[8192];
//    char buff1[8192];
    
//    memset (buff0, '\0', sizeof (buff0));
//    memset (buff1, '\0', sizeof (buff1));

//    setvbuf (MVs_theta_phi, buff0, _IOFBF, 8192);
//    setvbuf (FVs_theta_phi, buff2, _IOFBF, 8192);

    for (int l = 2; l <= LMAX; l++)
    {
      mpf_t coef_real[(2 * l) + 1], coef_imag[(2 * l) + 1];
      double raiz_real[2 * l], raiz_imag[2 * l];
      double theta[2 * l], phi[2 * l];
      coefi_pol (l, al_real, al_imag, coef_real, coef_imag);
      raizes_pol (l, coef_real, coef_imag, raiz_real, raiz_imag);

      for (int i = 0; i < ((2 * l) + 1); i++)
      {
        mpf_clears (coef_real[i], coef_imag[i], NULL);
      } 

      coord_pol (l, raiz_real, raiz_imag, theta, phi);
      frechet_pol (l, theta, phi,  &frechet_vec_theta, &frechet_vec_phi);

      int k = 0;
      for(int j = 0; j < 2 * l; j++)
      {
        k = (pow((l-1),2) + (l-1) - 2) + j;
        mvs_theta[k] = theta[j];
        mvs_phi[k] = phi[j];
      }

      printf ("MC %i -- Calculado ell %i\r\r\r\r\r\r\r\r\r\r\r\r\r", i, l);
      fprintf (FVs_theta_phi, "%f %f\n", frechet_vec_theta, frechet_vec_phi);
    }

    FILE *MVs_theta_phi = fopen (filename, "a"); 
    for (int j = 0; j < MVS_NUMERO; j++)
    {
      fprintf (MVs_theta_phi, "%.15lf %.15lf\n", mvs_theta[j], mvs_phi[j]);

    }
    fclose (MVs_theta_phi);
    fclose (FVs_theta_phi);
  }
}

int
main (int argc, char *argv[ ])
{
  int mc = atoi(argv[1]);

  char filename[300];
  
  static mpf_t al_real[ALMS_NUMERO_LINHAS];
  static mpf_t al_imag[ALMS_NUMERO_LINHAS];

  sprintf(filename,"./alm_fullsky_smica_noise_%d_.dat", mc);

  FILE *file = fopen (filename, "r");
  for (int j = 0; j < ALMS_NUMERO_LINHAS; j++)
  {
    mpf_inits (al_real[j], al_imag[j], NULL);
    gmp_fscanf (file, "%Ff %Ff\n", &al_real[j], &al_imag[j]);
  }
  
  fclose (file);
  printf ("Importado MC %i\n", mc);

  multipol_vec (mc, al_real, al_imag);

  for (int al_num = 0; al_num < ALMS_NUMERO_LINHAS; al_num++)
  {
    mpf_clears (al_real[al_num], al_imag[al_num], NULL);
  }

  printf ("Calculado MC %i\n", mc);
}
