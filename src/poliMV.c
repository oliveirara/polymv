#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>
#include <mps/mps.h>
#include <nlopt.h>
#include <complex.h>
#include <time.h>
#define o 1127251
#define als_numero 1 /* Numero de als gerados */
#define lmax 1500
#define NPIX 196608 /* NSIDE = 128 */

/* Função para calcular os coeficientes */
void
coefi_pol (int l, mpf_t al_real[], mpf_t al_imag[], mpf_t coef_real[], mpf_t coef_imag[])
{
  mpz_t int_binom;
  mpf_t float_binom;
  mpf_t raiz_binom;
  mpf_t coef_negativ;
  int index; /* levar em conta que esse valor do index pode explodir para listas grandes */

  mpz_inits (int_binom, NULL);
  mpf_inits (float_binom, raiz_binom, coef_negativ, NULL);

  mpf_inits (coef_real[l], coef_imag[l], NULL);

  mpz_bin_uiui (int_binom, (2 * l), l); /* Calculando Binomio */
  mpf_set_z (float_binom, int_binom);   /* Convertendo inteiro em float para tirar a raiz */
  mpf_sqrt (raiz_binom, float_binom);   /* Tirando a raiz */

  mpf_mul (coef_real[l], raiz_binom, al_real[l]); /* Multiplicando o al pela raiz do binom */
  mpf_mul (coef_imag[l], raiz_binom, al_imag[l]);


  for (int m = 1; m <= l; m++)
  {
    index = ((m * ((2 * lmax) + 1 - m)) / 2) + l; /*calculando a posição dos alms */

    mpf_inits (coef_real[l - m], coef_imag[l - m], coef_real[l + m], coef_imag[l + m], NULL);

    mpz_bin_uiui (int_binom, (2 * l), (m + l)); /* Calculando o binômio */
    mpf_set_z (float_binom, int_binom);         /* Convertendo inteiro em float para tirar a raiz */
    mpf_sqrt (raiz_binom, float_binom);         /* Tirando a raiz */

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

/* Função para extrair as raízes */
void
raizes_pol (int l, mpf_t coef_real[], mpf_t coef_imag[], double *raiz_real, double *raiz_imag)
{
  mpq_t rat_coef_real[(2 * l) + 1], rat_coef_imag[(2 * l) + 1];

  for (int i = 0; i < (2 * l) + 1; i++)
  { /*convertendo float em rational (pois até então o algoritmo só resolve com racional) */
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

/* Função para encontrar as coordenadas dos vetores multipolares */
void
coord_pol (int l, double raiz_real[], double raiz_imag[], double *theta, double *phi)
{
/*    double theta[2*l], phi[2*l], R[2*l]; */
  double R[2 * l];
  double complex z[2 * l];

  for (int i = 0; i < (2 * l); i++)
  {
    z[i]     = raiz_real[i] + (raiz_imag[i] * I);
    R[i]     = cabs (z[i]);
    phi[i]   = carg (z[i]); /* Theta domain is (0, pi)        Phi domain is (0, 2*pi)  */
    theta[i] = 2 * atan (1 / R[i]);
  }
}

void
eta_phi_pol (int l, double theta[], double phi[], double *eta, double *varphi)
{
  for (int i = 0; i < (2 * l); i++)
  {
    eta[i]    = (1 - cos (theta[i]));
    varphi[i] = phi[i] / (2 * M_PI);
  }
}

typedef struct _FrechetData
{
  int l;
  double * restrict x;
  double * restrict y;
  double * restrict z;
} FrechetData;

#define ACOS(ab) acos (((ab) > 1.0) ? 1.0 : (((ab) < -1.0) ? -1.0 : (ab)))

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
    //printf (".", fechet_mean);
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

/*Função para encontrar as coordenadas dos vetores de Fréchet */
void
frechet_pol (int l, double * restrict theta, double * restrict phi, double * restrict theta_esfera, double * restrict phi_esfera, 
             double * restrict x_esfera, double * restrict y_esfera, double * restrict z_esfera,
             double *frechet_vec_theta, double *frechet_vec_phi)
{
  FILE *psi = fopen ("psis(l=1500).dat", "a");
//  double frechet_vec_theta;
//  double frechet_vec_phi;
  const int twol = 2 * l;
  double x[2 * l], y[2 * l], z[2 * l];
  double f, fechet_mean;
  

  for (int i = 0; i < twol; i++)
  {
    x[i] = sin (theta[i]) * cos (phi[i]);
    y[i] = sin (theta[i]) * sin (phi[i]);
    z[i] = cos (theta[i]);
  }

  f = 0.0;
  for (int j = 0; j < twol; j++)
  {
    const double ab      = (x[j] * x_esfera[0]) + (y[j] * y_esfera[0]) + (z[j] * z_esfera[0]);
    const double acos_ab = ACOS (ab);
    f += acos_ab * acos_ab;
  }

  fprintf (psi, "%f %f %f\n", f, theta_esfera[0], phi_esfera[0]);
  *frechet_vec_theta = theta_esfera[0];
  *frechet_vec_phi   = phi_esfera[0];

  for (int i = 1; i < NPIX; i++)
  {
    fechet_mean = 0.0;

    for (int j = 0; j < twol; j++)
    {
      const double ab      = (x[j] * x_esfera[i]) + (y[j] * y_esfera[i]) + (z[j] * z_esfera[i]);
      const double acos_ab = ACOS (ab);
      fechet_mean += acos_ab * acos_ab;
      
    }

    fprintf (psi, "%f %f %f\n", fechet_mean, theta_esfera[i], phi_esfera[i]);

    if (f > fechet_mean)
    {
      f                  = fechet_mean;
      *frechet_vec_theta = theta_esfera[i];
      *frechet_vec_phi   = phi_esfera[i];
    }
  }

  if (*frechet_vec_theta > (M_PI / 2))
  {
    *frechet_vec_theta = M_PI - *frechet_vec_theta;
    *frechet_vec_phi   = M_PI + *frechet_vec_phi;

    if (*frechet_vec_phi > (2 * M_PI))
      *frechet_vec_phi = *frechet_vec_phi - (2 * M_PI);
  }

  fclose (psi);

//  printf ("M1 % 22.15g % 22.15g % 22.15g\n", *frechet_vec_theta, *frechet_vec_phi, f);


//    *frechet_vec_eta = 1-cos(frechet_vec_theta); 
//    *frechet_vec_varphi = frechet_vec_phi/(2*M_PI); 
}


void
frechet_pol2 (int l, double * restrict theta, double * restrict phi, double *frechet_vec_theta, double *frechet_vec_phi)
{
//  double frechet_vec_theta;
//  double frechet_vec_phi;
  const int twol = 2 * l;
  double x[2 * l], y[2 * l], z[2 * l];
  double f;
  int i;
  FrechetData fd = {l, x, y, z};
  
  *frechet_vec_theta = 0.0;
  *frechet_vec_phi = 0.0;

  for (i = 0; i < twol; i++)
  {
    x[i] = sin (theta[i]) * cos (phi[i]);
    y[i] = sin (theta[i]) * sin (phi[i]);
    z[i] = cos (theta[i]);
  }

  {
//    double lb[2] = {-1.0 * M_PI,  -2.0 * M_PI};
//    double ub[2] = {+1.0 * M_PI,  +2.0 * M_PI};
    double lb[2] = {0,  0};
    double ub[2] = {M_PI,  +2.0 * M_PI};
    double s[2]  = {0.0, 0.0};
    double min_f = 1.0e300;
    int ntheta = 8;
    int nphi   = 8;

    nlopt_opt opt;
    
    //opt = nlopt_create (NLOPT_LD_SLSQP, 2);
    opt = nlopt_create (NLOPT_LN_NELDERMEAD, 2);

    nlopt_set_lower_bounds (opt, lb);
    nlopt_set_upper_bounds (opt, ub);
    nlopt_set_min_objective(opt, &frechet_pol_min, &fd);
    nlopt_set_xtol_rel(opt, 1.0e-10);
    
    for (i = 0; i <= ntheta; i++)
    {
      int j;
      s[0] = 0.0;
      
      for (j = 0; j <= nphi; j++)
      {
        s[0] = 1.0 * M_PI / ((ntheta) + 1.0) * (i + 1.0);
        s[1] = 2.0 * M_PI / ((nphi)   + 1.0) * (j + 1.0);

        if (nlopt_optimize(opt, s, &f) < 0) 
        {
          printf("nlopt failed!\n");
        }
        else 
        {
//          printf("M2 % 22.15g % 22.15g % 22.15g\n", s[0], s[1], f);
          if (f < min_f)
          {
            min_f = f;

              *frechet_vec_theta = s[0];
              *frechet_vec_phi   = s[1];
          }
        }
      }
    }

    if (*frechet_vec_theta > (M_PI / 2))
  {
    *frechet_vec_theta = M_PI - *frechet_vec_theta;
    *frechet_vec_phi   = M_PI + *frechet_vec_phi;

    if (*frechet_vec_phi > (2 * M_PI))
      *frechet_vec_phi = *frechet_vec_phi - (2 * M_PI);
  }

//   printf ("M2 % 22.15g % 22.15g % 22.15g\n", *frechet_vec_theta, *frechet_vec_phi, min_f);
    nlopt_destroy (opt);

//    *frechet_vec_eta = 1-cos(frechet_vec_theta); 
//    *frechet_vec_varphi = frechet_vec_phi/(2*M_PI);
  } 
}


void
multipol_vec (mpf_t al_real[], mpf_t al_imag[])
{
/*
  double theta_esfera[NPIX], phi_esfera[NPIX], x_esfera[NPIX], y_esfera[NPIX], z_esfera[NPIX];

  {
    FILE *theta_phi = fopen ("ThetaPhi(NSiDE128).dat", "r");

    for (int i = 0; i < NPIX; i++)
    {
      fscanf (theta_phi, "%lf %lf\n", &theta_esfera[i], &phi_esfera[i]);
    }

    fclose (theta_phi);
  }

  {
    FILE *x_y_z = fopen ("x_y_z(NSIDE128).dat", "r");

    for (int i = 0; i < NPIX; i++)
    {
      fscanf (x_y_z, "%lf %lf %lf\n", &x_esfera[i], &y_esfera[i], &z_esfera[i]);
      //printf ("TESTE: (%f %f %f) %f\n", x_esfera[i], y_esfera[i], z_esfera[i], sqrt (z_esfera[i]*z_esfera[i] + x_esfera[i]*x_esfera[i] + y_esfera[i]*y_esfera[i]));      
    }
    fclose (x_y_z);
  }
*/
  {
      

    FILE *MVs = fopen ("MVs(Smica2015).dat", "a");
    FILE *FVs = fopen ("FVs(Smica2015).dat", "a");
    char buff0[8192];
    char buff1[8192];
    
    memset (buff0, '\0', sizeof (buff0));
    memset (buff1, '\0', sizeof (buff1));

    setvbuf (MVs, buff0, _IOFBF, 8192);
    setvbuf (FVs, buff1, _IOFBF, 8192);


    for (int l = 2; l <= lmax; l++)
    {
//      clock_t t;

//      t = clock ();

      mpf_t coef_real[(2 * l) + 1], coef_imag[(2 * l) + 1];
      double raiz_real[2 * l], raiz_imag[2 * l];
      double theta[2 * l], phi[2 * l];
      double frechet_vec_theta, frechet_vec_phi;

      coefi_pol (l, al_real, al_imag, coef_real, coef_imag);
      raizes_pol (l, coef_real, coef_imag, raiz_real, raiz_imag);

      for (int i = 0; i < ((2 * l) + 1); i++)
      {
        mpf_clears (coef_real[i], coef_imag[i], NULL);
      }

      coord_pol (l, raiz_real, raiz_imag, theta, phi);
/*      eta_phi_pol(l, theta, phi, eta, varphi); */
//      frechet_pol (l, theta, phi, theta_esfera, phi_esfera, x_esfera, y_esfera, z_esfera, &frechet_vec_theta, &frechet_vec_phi);
//     printf ("M1: %f %f\n", frechet_vec_theta, frechet_vec_phi);
//      frechet_pol2 (l, theta, phi,  &frechet_vec_theta, &frechet_vec_phi);
//     printf ("M2: %f %f\n", frechet_vec_theta, frechet_vec_phi);
      //printf ("% 22.15g % 22.15g\n", frechet_vec_theta, frechet_vec_phi);

      for (int j = 0; j < 2 * l; j++)
      {
        fprintf (MVs, "%f %f\n", theta[j], phi[j]);
      }

      fprintf (FVs, "%f %f\n", frechet_vec_theta, frechet_vec_phi);
      printf ("Calculado l = %i\n", l);

//      t = clock () - t;
//      printf ("Tempo para cada l: %lf\n", ((double) t) / ((CLOCKS_PER_SEC)));
    }

    

    fclose (MVs);
    fclose (FVs);
  }
}

int
main ()
{
  clock_t t;

  t = clock ();

  FILE *file;

  static mpf_t al_real[o];
  static mpf_t al_imag[o];


  char filename[400];

  for(int i = 0; i < als_numero; i++)
  {
//    sprintf(filename,"/home/ricardo/Mestrado/Programas/Testes/TestesPoli2/histogramas/ALMs/AL%d.dat", i);
		file = fopen("Alms(Smica2015).dat", "r");

    for (int j = 0; j < o; j++)
    {
      mpf_inits (al_real[j], al_imag[j], NULL);
      gmp_fscanf (file, "%Ff %Ff\n", &al_real[j], &al_imag[j]);

    }
  
    fclose (file);
  
    multipol_vec (al_real, al_imag);

    for (int i = 0; i < o; i++)
    {
      mpf_clears (al_real[i], al_imag[i], NULL);
    }
  }
  t = clock () - t;
  printf ("Tempo de execucao: %lf\n", ((double) t) / ((CLOCKS_PER_SEC)));
}

