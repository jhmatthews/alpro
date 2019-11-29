from distutils.core import setup, Extension
import numpy

libraries = ['m', 'gsl', 'gslcblas']

# this should be a symbolic link 
library_dirs = ['gsl_lib']
include_dir = 'gsl_include'

# define the extension module
alpro_module = Extension('alpro',
                    include_dirs=[numpy.get_include(), "/Users/matthews/winds/python/include"],
                    libraries = libraries,
                    library_dirs = library_dirs,
                    sources = ['alpro.c', 'alpro_matrix.c'])

# run the setup
setup(ext_modules=[alpro_module])  