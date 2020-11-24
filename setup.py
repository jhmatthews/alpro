from distutils.core import setup, Extension
import numpy
import os

#os.environ["CC"] = "gcc-10"
libraries = ['m', 'gsl', 'gslcblas']

# IMPROVE this should be a symbolic link / or an env variable
library_dirs = ['gsl_lib']
include_dir = 'gsl_include'
os.environ["CFLAGS"] = "-O0"

# define the extension module
alpro_module = Extension('alprocore',
                    include_dirs=[numpy.get_include(), include_dir],
                    libraries = libraries,
                    library_dirs = library_dirs,
                    sources = ['alpro/core/alpro.c', 'alpro/core/alpro_matrix.c'])

# data files 
data_files = ["data/*.dat"]

# run the setup
setup(name = 'alpro',
	  version = '1.0',
	  packages = ["alpro"],
	  description = 'Propagation of axion-like particles through magnetic fields',
	  author_email = 'matthews@ast.cam.ac.uk',
	  author = 'James Matthews',
	  ext_modules=[alpro_module], 
	  py_modules=["alpro"],
	  package_data = {'alpro': data_files},
	  include_package_data=True)  
