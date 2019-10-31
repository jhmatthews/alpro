/*gsl matrix solver header files */
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_complex_math.h>
/* alpro.c */
void get_U0(gsl_matrix_complex *U0_new, double EVarray[3], gsl_matrix_complex *T1, gsl_matrix_complex *T2, gsl_matrix_complex *T3, double distance);
void get_deltas(double Deltas[], double mass, double energy, double M, double B, double omega_pl);
int get_eigenvalues(double EVarray[3], double Deltas[4]);
void get_T_matrices(double alpha, gsl_matrix_complex_view *T1, gsl_matrix_complex_view *T2, gsl_matrix_complex_view *T3);
void get_T_matrices2(double EVarray[3], double Deltas[4], gsl_matrix_complex_view *T1, gsl_matrix_complex_view *T2, gsl_matrix_complex_view *T3);
int apply_rotation_matrix(double phi, gsl_matrix_complex *U0);
