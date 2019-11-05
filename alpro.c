/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief ALPRO: Axion Like PROpagation
*/
/* ///////////////////////////////////////////////////////////////////// */
#include <math.h>
#include "prototypes.h"
#include "const.h"
#include <Python.h>
#include <numpy/arrayobject.h>


/*gsl matrix solver header files */
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_complex_math.h>

/* some define statements for Deltas array 
   -- act like subscript indices */
#define PL 0
#define AA 1
#define AG 2
#define OSC 3
#define NDELTAS 4

static PyObject* get_P(PyObject* self, PyObject* args)
{

  gsl_matrix_complex_view U0_matrix, T1, T2, T3;
  gsl_vector_complex_view A;
  gsl_vector_complex *A_new;
  gsl_matrix_complex *U0;
  gsl_complex exp_term;
  double mass, energy, M, B, omega_pl;
  double EVarray[3];
  double Deltas[NDELTAS];
  double alpha, distance, phi, logenergy, g_a, L;
  int ienergy, i;

  /* Variables for calling from Python */
  PyArrayObject *in_array;
  PyObject      *out_array;
  NpyIter *in_iter;
  NpyIter *out_iter;
  NpyIter_IterNextFunc *in_iternext;
  NpyIter_IterNextFunc *out_iternext;

  /*  parse single numpy array argument */
  if (!PyArg_ParseTuple(args, "O!ddddd", &PyArray_Type, &in_array, &phi, &B, &L, &g_a, &mass))
      return NULL;

  B *= UNIT_GAUSS;

  /*  construct the output array, like the input array */
  out_array = PyArray_NewLikeArray(in_array, NPY_ANYORDER, NULL, 0);
  if (out_array == NULL)
      return NULL;

  /*  create the iterators */
  in_iter = NpyIter_New(in_array, NPY_ITER_READONLY, NPY_KEEPORDER,
                           NPY_NO_CASTING, NULL);
  if (in_iter == NULL)
      goto fail;

  out_iter = NpyIter_New((PyArrayObject *)out_array, NPY_ITER_READWRITE,
                        NPY_KEEPORDER, NPY_NO_CASTING, NULL);
  if (out_iter == NULL) {
      NpyIter_Deallocate(in_iter);
      goto fail;
  }

  in_iternext = NpyIter_GetIterNext(in_iter, NULL);
  out_iternext = NpyIter_GetIterNext(out_iter, NULL);
  if (in_iternext == NULL || out_iternext == NULL) {
      NpyIter_Deallocate(in_iter);
      NpyIter_Deallocate(out_iter);
      goto fail;
  }
  double ** in_dataptr = (double **) NpyIter_GetDataPtrArray(in_iter);
  double ** out_dataptr = (double **) NpyIter_GetDataPtrArray(out_iter);

  /*  iterate over the arrays */
  do {
    double A_data[] = { 0.0, 1.0, 0.0, 0.0, 0.0, 0.0 };

    /* mass and coupling constant */
    M = 1.0 / g_a;

    /* set plasma frequency for n = 1e-2 */
    omega_pl = sqrt (4.0 * const_PI * const_ECHARGE * const_ECHARGE * 1e-2 / const_MELEC) * const_HEV;
    // logenergy = (3 + ((double) ienergy * dlogenergy));
    energy = **in_dataptr;

    /* 1.0 kpc */
    distance = L * 1000.0 * const_PC * UNIT_LENGTH;

    A = gsl_vector_complex_view_array (A_data, 3);
    U0 = gsl_matrix_complex_alloc (3, 3);
    A_new = gsl_vector_complex_alloc (3);


    /* get deltas and eigenvalues */
    get_deltas (Deltas, mass, energy, M, B, omega_pl);
    get_eigenvalues (EVarray, Deltas);

    /* calculate T matrices from mixing angle, alpha (eq 34) */
    alpha = 0.5 * atan (2.0 * Deltas[AG] / (Deltas[PL] - Deltas[AA]));
    get_T_matrices (alpha, &T1, &T2, &T3);
    // get_T_matrices2 (EVarray, Deltas, &T1, &T2, &T3);


    /* construct the transfer matrix for the idealised parallel situation */
    get_U0 (U0, EVarray, &T1.matrix, &T2.matrix, &T3.matrix, distance);

    /* apply the rotation matrix */
    apply_rotation_matrix (phi, U0);

    /* multiply the state vector (A) by the transfer matrix (U0) 
       results is stored in A_new */
    exp_term = gsl_complex_polar (1.0, energy * distance);
    gsl_matrix_complex_scale (U0, exp_term);
    gsl_blas_zgemv (CblasNoTrans, GSL_COMPLEX_ONE, U0, &A.vector, GSL_COMPLEX_ZERO, A_new);


    // printf ("Vectors: %8.4e ", energy);
    for (i = 0; i < 3; i++)
    {
      gsl_complex x = gsl_vector_complex_get (A_new, i);
      double abs2 = gsl_complex_abs2 (x);
      // printf ("%8.4e ", abs2);

      if (i == 2) **out_dataptr = abs2;
      // printf()
    }
    // printf ("\n");

    gsl_vector_complex_free (A_new);
    gsl_matrix_complex_free (U0);
  }
  while (in_iternext (in_iter) && out_iternext (out_iter));

  /*  clean up and return the result */
  NpyIter_Deallocate(in_iter);
  NpyIter_Deallocate(out_iter);
  Py_INCREF(out_array);
  return out_array;

  /*  in case bad things happen */
  fail:
    Py_XDECREF(out_array);
    return NULL;
}

/*  define functions in module */
static PyMethodDef AlproMethods[] =
{
     {"get_P", get_P, METH_VARARGS,
         "evaluate the photon to axion probability"},
     {NULL, NULL, 0, NULL}
};


#if PY_MAJOR_VERSION >= 3
/* module initialization */
/* Python version 3*/
static struct PyModuleDef cModPyDem = {
    PyModuleDef_HEAD_INIT,
    "alpro_module", "Some documentation",
    -1,
    AlproMethods
};
PyMODINIT_FUNC PyInit_alpro(void) {
    PyObject *module;
    module = PyModule_Create(&cModPyDem);
    if(module==NULL) return NULL;
    /* IMPORTANT: this must be called */
    import_array();
    if (PyErr_Occurred()) return NULL;
    return module;
}

#else
/* module initialization */
/* Python version 2 */
PyMODINIT_FUNC initalpro(void) {
    PyObject *module;
    module = Py_InitModule("alpro", AlproMethods);
    if(module==NULL) return;
    /* IMPORTANT: this must be called */
    import_array();
    return;
}

#endif

void
get_U0 (gsl_matrix_complex * U0_new, double EVarray[3],
        gsl_matrix_complex * T1, gsl_matrix_complex * T2, gsl_matrix_complex * T3, double distance)
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
get_deltas (double Deltas[], double mass, double energy, double M, double B, double omega_pl)
{
  double x;
  Deltas[PL] = -(omega_pl * omega_pl) / 2.0 / energy;
  Deltas[AA] = -(mass * mass) / 2.0 / energy;
  Deltas[AG] = B / 2.0 / M;

  /* now get Delta_osc */
  x = (Deltas[PL] - Deltas[AA]) * (Deltas[PL] - Deltas[AA]);
  Deltas[OSC] = sqrt (x + (4.0 * Deltas[AG] * Deltas[AG]));
}

int
get_eigenvalues (double EVarray[3], double Deltas[NDELTAS])
{
  /* get eigenvalues of mixing matrix */
  EVarray[0] = Deltas[PL];
  EVarray[1] = 0.5 * (Deltas[PL] + Deltas[AA] - Deltas[OSC]);
  EVarray[2] = 0.5 * (Deltas[PL] + Deltas[AA] + Deltas[OSC]);
  return 0;
}

void
get_T_matrices (double alpha, gsl_matrix_complex_view * T1, gsl_matrix_complex_view * T2, gsl_matrix_complex_view * T3)
{
  int i, j;
  int nrows = 3;
  double T1_temp[nrows][nrows], T2_temp[nrows][nrows], T3_temp[nrows][nrows];
  double T1_data[2 * nrows * nrows], T2_data[2 * nrows * nrows], T3_data[2 * nrows * nrows];

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

  /* zero the matrixes */
  for (i = 0; i < nrows; i++)
  {
    for (j = 0; j < nrows; j++)
    {
      T1_data[2 * i * nrows + j] = T1_temp[i][j];
      T2_data[2 * i * nrows + j] = T2_temp[i][j];
      T3_data[2 * i * nrows + j] = T3_temp[i][j];
      T1_data[1 + (i * nrows) + j] = 0.0;
      T2_data[1 + (i * nrows) + j] = 0.0;
      T3_data[1 + (i * nrows) + j] = 0.0;
    }
  }
  *T1 = gsl_matrix_complex_view_array (T1_data, nrows, nrows);
  *T2 = gsl_matrix_complex_view_array (T2_data, nrows, nrows);
  *T3 = gsl_matrix_complex_view_array (T3_data, nrows, nrows);
}

void
get_T_matrices2 (double EVarray[3], double Deltas[4], gsl_matrix_complex_view * T1, gsl_matrix_complex_view * T2,
                 gsl_matrix_complex_view * T3)
{
  int i, j;
  int nrows = 3;
  double T1_temp[nrows][nrows], T2_temp[nrows][nrows], T3_temp[nrows][nrows];
  double T1_data[2 * nrows * nrows], T2_data[2 * nrows * nrows], T3_data[2 * nrows * nrows];

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
  T2_temp[1][1] = T3_temp[2][2] = (EVarray[2] - t) / (EVarray[2] - EVarray[1]);
  T2_temp[2][1] = (EVarray[2] - t) * (EVarray[1] - t) / v / (EVarray[2] - EVarray[1]);
  T2_temp[1][2] = -v / (EVarray[2] - EVarray[1]);
  T2_temp[2][2] = T3_temp[1][1] = -(EVarray[1] - t) / (EVarray[2] - EVarray[1]);
  T3_temp[2][1] = -T2_temp[2][1];
  T3_temp[1][2] = -T2_temp[1][2];


  /* zero the matrixes */
  for (i = 0; i < nrows; i++)
  {
    for (j = 0; j < nrows; j++)
    {
      T1_data[2 * i * nrows + j] = T1_temp[i][j];
      T2_data[2 * i * nrows + j] = T2_temp[i][j];
      T3_data[2 * i * nrows + j] = T3_temp[i][j];
      T1_data[1 + (i * nrows) + j] = 0.0;
      T2_data[1 + (i * nrows) + j] = 0.0;
      T3_data[1 + (i * nrows) + j] = 0.0;
    }
  }
  *T1 = gsl_matrix_complex_view_array (T1_data, nrows, nrows);
  *T2 = gsl_matrix_complex_view_array (T2_data, nrows, nrows);
  *T3 = gsl_matrix_complex_view_array (T3_data, nrows, nrows);
}

// void print_complex_matrix(gsl_matrix_complex matrix)
// {
//   for (i = 0; i < nrows; i++)
//     {
//       print ()
//     }
// }


int
apply_rotation_matrix (double phi, gsl_matrix_complex * U0)
{

  int i, j;
  gsl_matrix_complex *C, *v_transpose;  /* coefficient matrix A'A  */
  gsl_matrix_complex_view v_matrix;
  int nrows = 3;
  double v_data[nrows * nrows * 2], v[nrows][nrows];

  /* allocate matrices */
  C = gsl_matrix_complex_alloc (nrows, nrows);  /* Data matrix */
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
      v_data[2 * (i * nrows + j)] = v[i][j];

      /* set the complex part to zero */
      v_data[1 + (i * nrows + j)] = 0.0;
    }
  }

  /* create matrix objects */
  /* need to fix this so the array of doubles is right!! */
  v_matrix = gsl_matrix_complex_view_array (v_data, nrows, nrows);

  /* calculate the transpose of the rotation matrix */
  gsl_matrix_complex_transpose_memcpy (v_transpose, &v_matrix.matrix);

  /* do the actual matrix multiplication */
  /* first multiply U0 by V and store in C */
  gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, &v_matrix.matrix, U0, GSL_COMPLEX_ZERO, C);

  /* now multiply V transpose by C, and put back in U0 */
  gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, C, v_transpose, GSL_COMPLEX_ZERO, U0);

  /* now we have the required matrix stored in U0_matrix and we can return */

  /* free memory? should I make this more efficient? */
  gsl_matrix_complex_free (C);
  gsl_matrix_complex_free (v_transpose);

  return 0;
}
