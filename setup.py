from distutils.core import setup, Extension
#from Cython.Build import cythonize
import numpy

libraries = ['m', 'gsl', 'gslcblas']
library_dirs = ['/Users/matthews/winds/python/lib']

# define the extension module
alpro_module = Extension('alpro',
                    include_dirs=[numpy.get_include(), "/Users/matthews/winds/python/include"],
                    libraries = libraries,
                    library_dirs = library_dirs,
                    sources = ['alpro.c'])

# run the setup
setup(ext_modules=[alpro_module])  