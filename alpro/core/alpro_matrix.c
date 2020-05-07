/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief ALPRO: Axion Like PROpagation
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "alpro.h"
#include "prototypes.h"

void
get_U0 (gsl_matrix_complex * U0_new, double EVarray[3],
	gsl_matrix_complex * T1, gsl_matrix_complex * T2,
	gsl_matrix_complex * T3, double distance)
{
  gsl_complex EVexp1, EVexp2, EVexp3, xx;
  int i, j;

  /* get complex coefficients (e^i lambda_j d) */
  EVexp1 = gsl_complex_polar (1.0, EVarray[0] * distance);
  EVexp2 = gsl_complex_polar (1.0, EVarray[1] * distance);
  EVexp3 = gsl_complex_polar (1.0, EVarray[2] * distance);

  /* multiply complex coefficients by the matrices */
  gsl_matrix_complex_scale (T1, EVexp1);
  gsl_matrix_complex_scale (T2, EVexp2);
  gsl_matrix_complex_scale (T3, EVexp3);

  /* add the terms together and copy into the new matrix */
  /* note that after this process T1 T2 and T3 have been messed with! */
  gsl_matrix_complex_add (T1, T2);
  gsl_matrix_complex_add (T1, T3);
  gsl_matrix_complex_memcpy (U0_new, T1);
}


void
get_deltas (double Deltas[], double mass, double energy, double M, double B,
	    double omega_pl)
{
  double x;
  Deltas[PL] = -(omega_pl * omega_pl) / 2.0 / energy * UNIT_DELTA;
  Deltas[AA] = -(mass * mass) / 2.0 / energy * UNIT_DELTA;
  Deltas[AG] = B / 2.0 / M * UNIT_DELTA;

  /* now get Delta_osc */
  x = (Deltas[PL] - Deltas[AA]) * (Deltas[PL] - Deltas[AA]);
  Deltas[OSC] = sqrt (x + (4.0 * Deltas[AG] * Deltas[AG])) * UNIT_DELTA;

  // int i;
  // for (i == 0; i < NDELTAS; i++)
  // {
  //   Deltas[i] = Deltas[i]*1e30;
  // }

}

int
get_eigenvalues (double EVarray[3], double Deltas[NDELTAS])
{
  /* get eigenvalues of mixing matrix */
  EVarray[0] = Deltas[PL] / UNIT_DELTA;
  EVarray[1] = 0.5 * (Deltas[PL] + Deltas[AA] - Deltas[OSC]) / UNIT_DELTA;
  EVarray[2] = 0.5 * (Deltas[PL] + Deltas[AA] + Deltas[OSC]) / UNIT_DELTA;
  return 0;
}

void
get_T_matrices (double alpha, double EVarray[3], double Deltas[4], gsl_matrix_complex_view * T1,
		gsl_matrix_complex_view * T2, gsl_matrix_complex_view * T3)
{
  int i, j;
  int nrows = 3;
  double T1_temp[nrows][nrows], T2_temp[nrows][nrows], T3_temp[nrows][nrows];
  double T1_data[2 * nrows * nrows], T2_data[2 * nrows * nrows],
    T3_data[2 * nrows * nrows];
  double t, u, v, s;
  t = Deltas[PL];
  v = Deltas[AG];
  u = Deltas[AA];
  /* zero the matrixes */
  for (i = 0; i < nrows; i++)
    {
      for (j = 0; j < nrows; j++)
	{
	  T1_temp[i][j] = T2_temp[i][j] = T3_temp[i][j] = 0.0;
	}
    }

  /* first matrix, eq 35 */
  T1_temp[0][0] = 1.0;

  /* second and third matrixes, eqs 36 and 37 */
  T2_temp[1][1] = T3_temp[2][2] = sin (alpha) * sin (alpha);
  T2_temp[2][1] = T2_temp[1][2] = -sin (alpha) * cos (alpha);
  T3_temp[2][1] = T3_temp[1][2] = sin (alpha) * cos (alpha);
  T2_temp[2][2] = T3_temp[1][1] = cos (alpha) * cos (alpha);


  // printf("T2 Matrix a:\n");
  // print_matrix(T2_temp, 3);

  //   /* zero the matrixes */
  // for (i = 0; i < nrows; i++)
  //   {
  //     for (j = 0; j < nrows; j++)
  // {
  //   T1_temp[i][j] = T2_temp[i][j] = T3_temp[i][j] = 0.0;
  // }
  //   }


  //   /* second and third matrixes, eqs 36 and 37 */
  // T2_temp[1][1] = T3_temp[2][2] =
  //   (EVarray[2] - t) / (EVarray[2] - EVarray[1]);
  // T2_temp[2][1] =
  //   (EVarray[2] - t) * (EVarray[1] - t) / v / (EVarray[2] - EVarray[1]);
  // T2_temp[1][2] = -v / (EVarray[2] - EVarray[1]);
  // T2_temp[2][2] = T3_temp[1][1] =
  //   -(EVarray[1] - t) / (EVarray[2] - EVarray[1]);
  // T3_temp[2][1] = -T2_temp[2][1];
  // T3_temp[1][2] = -T2_temp[1][2];

  // printf("T2 Matrix b:\n");
  // print_matrix(T2_temp, 3);

  // exit(0);

  /* populate the matrixes */
  for (i = 0; i < nrows; i++)
    {
      for (j = 0; j < nrows; j++)
	{
	  T1_data[(2 * i * nrows) + (2 * j)] = T1_temp[i][j];
	  T2_data[(2 * i * nrows) + (2 * j)] = T2_temp[i][j];
	  T3_data[(2 * i * nrows) + (2 * j)] = T3_temp[i][j];
	  T1_data[1 + (2 * i * nrows) + (2 * j)] = 0.0;
	  T2_data[1 + (2 * i * nrows) + (2 * j)] = 0.0;
	  T3_data[1 + (2 * i * nrows) + (2 * j)] = 0.0;
	}
    }
  *T1 = gsl_matrix_complex_view_array (T1_data, nrows, nrows);
  *T2 = gsl_matrix_complex_view_array (T2_data, nrows, nrows);
  *T3 = gsl_matrix_complex_view_array (T3_data, nrows, nrows);
}

void
get_T_matrices2 (double EVarray[3], double Deltas[4],
		 gsl_matrix_complex_view * T1, gsl_matrix_complex_view * T2,
		 gsl_matrix_complex_view * T3)
{
  int i, j;
  int nrows = 3;
  double T1_temp[nrows][nrows], T2_temp[nrows][nrows], T3_temp[nrows][nrows];
  double T1_data[2 * nrows * nrows], T2_data[2 * nrows * nrows],
    T3_data[2 * nrows * nrows];

  double t, u, v, s;
  t = Deltas[PL];
  v = Deltas[AG];
  u = Deltas[AA];

  /* zero the matrixes */
  for (i = 0; i < nrows; i++)
    {
      for (j = 0; j < nrows; j++)
	{
	  T1_temp[i][j] = T2_temp[i][j] = T3_temp[i][j] = 0.0;
	}
    }

  /* first matrix, eq 35 */
  T1_temp[0][0] = 1.0;

  /* second and third matrixes, eqs 36 and 37 */
  T2_temp[1][1] = T3_temp[2][2] =
    (EVarray[2] - t) / (EVarray[2] - EVarray[1]);
  T2_temp[2][1] =
    (EVarray[2] - t) * (EVarray[1] - t) / v / (EVarray[2] - EVarray[1]);
  T2_temp[1][2] = -v / (EVarray[2] - EVarray[1]);
  T2_temp[2][2] = T3_temp[1][1] =
    -(EVarray[1] - t) / (EVarray[2] - EVarray[1]);
  T3_temp[2][1] = -T2_temp[2][1];
  T3_temp[1][2] = -T2_temp[1][2];


  /* zero the matrixes */
  for (i = 0; i < nrows; i++)
    {
      for (j = 0; j < nrows; j++)
	{
	  T1_data[(2 * i * nrows) + (2 * j)] = T1_temp[i][j];
	  T2_data[(2 * i * nrows) + (2 * j)] = T2_temp[i][j];
	  T3_data[(2 * i * nrows) + (2 * j)] = T3_temp[i][j];
	  T1_data[1 + (2 * i * nrows) + (2 * j)] = 0.0;
	  T2_data[1 + (2 * i * nrows) + (2 * j)] = 0.0;
	  T3_data[1 + (2 * i * nrows) + (2 * j)] = 0.0;
	}
    }
  *T1 = gsl_matrix_complex_view_array (T1_data, nrows, nrows);
  *T2 = gsl_matrix_complex_view_array (T2_data, nrows, nrows);
  *T3 = gsl_matrix_complex_view_array (T3_data, nrows, nrows);
}



int
apply_rotation_matrix (double phi, gsl_matrix_complex * U0)
{

  int i, j;
  gsl_matrix_complex *C, *v_transpose;	/* coefficient matrix A'A  */
  gsl_matrix_complex_view v_matrix;
  int nrows = 3;
  double v_data[nrows * nrows * 2], v[nrows][nrows];

  /* allocate matrices */
  C = gsl_matrix_complex_alloc (nrows, nrows);	/* Data matrix */
  v_transpose = gsl_matrix_complex_alloc (nrows, nrows);

  /* first create rotation matrix. Zero elements */
  for (i = 0; i < nrows; i++)
    {
      for (j = 0; j < nrows; j++)
	{
	  v[i][j] = 0.0;
	}
    }

  /* populate elements of rotation matrix */
  v[0][0] = cos (phi);
  v[0][1] = -sin (phi);
  v[1][0] = sin (phi);
  v[1][1] = cos (phi);
  v[2][2] = 1.0;

  /* flatten the 2D array into a 1D array */
  for (i = 0; i < nrows; i++)
    {
      for (j = 0; j < nrows; j++)
	{
	  v_data[(2 * i * nrows) + (2 * j)] = v[i][j];

	  /* set the complex part to zero */
	  v_data[1 + (2 * i * nrows) + (2 * j)] = 0.0;
	}
    }

  /* create matrix objects */
  /* need to fix this so the array of doubles is right!! */
  v_matrix = gsl_matrix_complex_view_array (v_data, nrows, nrows);

  /* calculate the transpose of the rotation matrix */
  gsl_matrix_complex_transpose_memcpy (v_transpose, &v_matrix.matrix);

  /* do the actual matrix multiplication */
  /* first multiply U0 by V and store in C */
  gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE,
		  &v_matrix.matrix, U0, GSL_COMPLEX_ZERO, C);

  /* now multiply V transpose by C, and put back in U0 */
  gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, C, v_transpose,
		  GSL_COMPLEX_ZERO, U0);

  /* now we have the required matrix stored in U0_matrix and we can return */
  gsl_matrix_complex_free (C);
  gsl_matrix_complex_free (v_transpose);

  return 0;
}

/* Create 1D Carray from PyArray
   Assumes PyArray is contiguous in memory. Credit Lou Pecora     */
double *
pyvector_to_Carrayptrs (PyArrayObject * arrayin)
{
  int i, n;

  n = arrayin->dimensions[0];
  return (double *) arrayin->data;	/* pointer to arrayin data as double */
}


void print_matrix(double matrix[3][3], int nrows)
{
  int i, j;
  for (i = 0; i < nrows; i++)
    {
      for (j = 0; j < nrows; j++)
  {
    printf("%8.4e ", matrix[i][j]);
  }
    printf("\n");
    }
}

