/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief ALPRO: Axion Like PROpagation
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "alpro.h"
#include "prototypes.h"

/* improve this to get single and multiple probabilities */
static PyObject *
get_P (PyObject * self, PyObject * args)
{
  double mass, energy, M, B, omega_pl;
  double alpha, distance, phi, logenergy, g_a, L, ne;
  double *A_out, *A_init;
  double *A_temp;
  int i;
  gsl_vector_complex *A_new;

  /* Variables for calling from Python */
  PyArrayObject *in_array, *A_array;
  PyObject *out_array, *Aout_array;
  NpyIter *in_iter;
  NpyIter *out_iter;
  NpyIter_IterNextFunc *in_iternext;
  NpyIter_IterNextFunc *out_iternext;
  A_temp = calloc (sizeof (double), 6);

  /*  parse single numpy array argument with a series of double variables after it */
  if (!PyArg_ParseTuple
      (args, "O!O!dddddd", &PyArray_Type, &in_array, &PyArray_Type, &A_array,
       &phi, &B, &L, &g_a, &mass, &ne))
    return NULL;

  /* B should be given in Gauss, so convert to natural units */
  B *= UNIT_GAUSS;

  /* L kpc, converted to natural units */
  distance = L * 1000.0 * const_PC * UNIT_LENGTH;

  /* coupling constant in units used by De Angelis et al. */
  M = 1.0 / g_a;

  /* set plasma frequency for n = 1e-2 */
  omega_pl =
    sqrt (4.0 * const_PI * const_ECHARGE * const_ECHARGE * ne /
	  const_MELEC) * const_HBAR_EV;

  /*  output arrays to match the input arrays */
  out_array = PyArray_NewLikeArray (in_array, NPY_ANYORDER, NULL, 0);
  if (out_array == NULL)
    return NULL;
  Aout_array = PyArray_NewLikeArray (A_array, NPY_ANYORDER, NULL, 0);
  if (Aout_array == NULL)
    return NULL;

  /* make the iterators */
  in_iter = NpyIter_New (in_array, NPY_ITER_READONLY, NPY_KEEPORDER,
			 NPY_NO_CASTING, NULL);
  if (in_iter == NULL)
    goto fail;

  // A_iter = NpyIter_New (A_array, NPY_ITER_READONLY, NPY_KEEPORDER,
  //      NPY_NO_CASTING, NULL);

  out_iter = NpyIter_New ((PyArrayObject *) out_array, NPY_ITER_READWRITE,
			  NPY_KEEPORDER, NPY_NO_CASTING, NULL);
  if (out_iter == NULL)
    {
      NpyIter_Deallocate (in_iter);
      goto fail;
    }

  in_iternext = NpyIter_GetIterNext (in_iter, NULL);
  out_iternext = NpyIter_GetIterNext (out_iter, NULL);
  if (in_iternext == NULL || out_iternext == NULL)
    {
      NpyIter_Deallocate (in_iter);
      NpyIter_Deallocate (out_iter);
      goto fail;
    }
  double **in_dataptr = (double **) NpyIter_GetDataPtrArray (in_iter);
  double **out_dataptr = (double **) NpyIter_GetDataPtrArray (out_iter);
  /* done making iterators for energy */

  /* initial and output A_arrays */
  A_new = gsl_vector_complex_alloc (3);
  A_init = pyvector_to_Carrayptrs (A_array);
  A_out = pyvector_to_Carrayptrs (Aout_array);
  int ienergy = 0;

  // for (i = 0; i < 6; i++)
  //   A_out[i] = A_init[i];

  /*  iterate over the energy arrays */
  do
    {
      for (i = 0; i < 6; i++)
      {
	      A_temp[i] = A_init[(ienergy * 6) + i];
        printf ("%8.4e %d\n", A_temp[i], (ienergy * 6) + i);
      }
      printf("\n");

      /* get energy from array */
      energy = **in_dataptr;

      /* do the propagation through one domain */
      PropagateOne (A_new, mass, energy, M, B, omega_pl, phi, distance,
		    A_temp);

      /* copy probabilities to output array */
      /* need to modify to output array as well */
      for (i = 0; i < 3; i++)
	{
	  gsl_complex x = gsl_vector_complex_get (A_new, i);
	  A_out[(ienergy * 6) + i] = x.dat[0];
	  A_out[(ienergy * 6) + 3 + i] = x.dat[1];

	  double abs2 = gsl_complex_abs2 (x);

	  if (i == 2)
	    **out_dataptr = abs2;
	}
      for (i = 0; i < 6; i++)
      {
        printf ("JM: %8.4e %8.4e\n", A_out[(ienergy * 6) + i], **out_dataptr);
      }
      ienergy++;
    }
  while (in_iternext (in_iter) && out_iternext (out_iter));

  /*  clean up and return the result */
  gsl_vector_complex_free (A_new);
  NpyIter_Deallocate (in_iter);
  NpyIter_Deallocate (out_iter);
  Py_INCREF (out_array);
  return Py_BuildValue("OO", out_array, Aout_array);

  /*  in case bad things happen */
fail:
  Py_XDECREF (out_array);
  return NULL;
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
  alpha = 0.5 * atan (2.0 * Deltas[AG] / (Deltas[PL] - Deltas[AA]));
  get_T_matrices (alpha, &T1, &T2, &T3);
  /* alternative approach commented out */
  // get_T_matrices2 (EVarray, Deltas, &T1, &T2, &T3);


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















/* Python interfacing stuff */

/*  define functions in module */
static PyMethodDef AlproMethods[] = {
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

PyMODINIT_FUNC
PyInit_alpro (void)
{
  PyObject *module;
  module = PyModule_Create (&cModPyDem);
  if (module == NULL)
    return NULL;
  /* IMPORTANT: this must be called */
  import_array ();
  if (PyErr_Occurred ())
    return NULL;
  return module;
}

#else
/* module initialization */
/* Python version 2 */
PyMODINIT_FUNC
initalpro (void)
{
  PyObject *module;
  module = Py_InitModule ("alpro", AlproMethods);
  if (module == NULL)
    return;
  /* IMPORTANT: this must be called */
  import_array ();
  return;
}

#endif
