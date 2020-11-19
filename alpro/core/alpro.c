/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief ALPRO: Axion Like PROpagation
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "alpro.h"
#include "prototypes.h"
int main(int argc, char const *argv[]);

int main(int argc, char const *argv[])
{
  double mass, energy, M, B, omega_pl;
  double alpha, distance, phi, logenergy, g_a, L, ne;
  double *A_out, *A_init, *phis;
  double *A_temp;
  int i;
  gsl_vector_complex *A_new;

  mass = 1.0;
  energy = 1.0;
  M = 1.0;
  B = 1.0;
  omega_pl = 1.0;
  distance = 1.0;
  A_temp = calloc (sizeof (double), 6);
  A_new = gsl_vector_complex_alloc (3);

  PropagateOne (A_new, mass, energy, M, B, omega_pl, phi, distance, A_temp);

  return (1.0);
}


void
PropagateOne (gsl_vector_complex * A_new, double mass, double energy,
	      double M, double B, double omega_pl, double phi,
	      double distance, double *A_data)
{
  gsl_vector_complex_view A;
  gsl_matrix_complex_view U0_matrix, T1, T2, T3;
  gsl_matrix_complex *U0;
  gsl_complex exp_term;
  double EVarray[3];
  double Deltas[NDELTAS], alpha;


  A = gsl_vector_complex_view_array (A_data, 3);
  U0 = gsl_matrix_complex_alloc (3, 3);

  // A_new = gsl_vector_complex_alloc (3);
  /* get deltas and eigenvalues */
  get_deltas (Deltas, mass, energy, M, B, omega_pl);
  get_eigenvalues (EVarray, Deltas);
  /* calculate T matrices from mixing angle, alpha (eq 3 */
  // alpha = 0.5 * atan (2.0 * Deltas[AG] / (Deltas[PL] - Deltas[AA]));
  // printf("ALPHA is %8.4e Deltas %8.4e %8.4e %8.4e\n", alpha, Deltas[AG], Deltas[PL], Deltas[AA]);
  alpha = 0.5 * atan2 (2.0 * Deltas[AG], (Deltas[PL] - Deltas[AA]));
  // printf("ALPHA is %8.4e Deltas %8.4e %8.4e %8.4e\n", alpha, Deltas[AG], Deltas[PL], Deltas[AA]);
  //alpha = 0.5 * atan(B / M * 2.0 * energy / ((mass * mass) - (omega_pl * omega_pl))); 
  //printf("ALPHA is %8.4e Deltas %8.4e %8.4e %8.4e\n", alpha, Deltas[AG], Deltas[PL], Deltas[AA])
  get_T_matrices (alpha, EVarray, Deltas, &T1, &T2, &T3);
  /* alternative approach commented out */
  //get_T_matrices2 (EVarray, Deltas, &T1, &T2, &T3);
  distance = distance * UNIT_DELTA;
  /* construct the transfer matrix for the idealised parallel situation */
  get_U0 (U0, EVarray, &T1.matrix, &T2.matrix, &T3.matrix, distance);

  /* apply the rotation matrix */
  apply_rotation_matrix (phi, U0);

  /* multiply the state vector (A) by the transfer matrix (U0) 
     result is stored in A_new */
  exp_term = gsl_complex_polar (1.0, energy * distance);
  gsl_matrix_complex_scale (U0, exp_term);

  gsl_blas_zgemv (CblasNoTrans, GSL_COMPLEX_ONE, U0, &A.vector,
		  GSL_COMPLEX_ZERO, A_new);

  gsl_matrix_complex_free (U0);
}










