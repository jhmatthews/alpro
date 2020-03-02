/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief ALPRO: Axion Like PROpagation header file
*/
/* ///////////////////////////////////////////////////////////////////// */
#include <math.h>
// #include "prototypes.h"
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

// typedef struct domain 
// {
// 	double B;
// 	double ne, omega_pl;
// 	double L;
// 	double M, g_a;
// 	double distance;
// 	double mass;
// 	double phi;
// } domain_dummy, *DomainPtr;

// typedef struct test 
// {
// 	double A;
// } test_dummy, *TestPtr;
