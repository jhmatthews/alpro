from distutils.core import setup, Extension
import numpy

libraries = ['m', 'gsl', 'gslcblas']

# IMPROVE this should be a symbolic link / or an env variable
library_dirs = ['../gsl_lib']
include_dir = '../gsl_include'

# define the extension module
# define the extension module
alpro_module = Extension('alprocore',
                    include_dirs=[numpy.get_include(), include_dir],
                    libraries = libraries,
                    library_dirs = library_dirs,
                    sources = ['test_gsl/core/alpro.c', 'test_gsl/core/alpro_matrix.c'])

# data files 
data_files = ["data/*.dat"]

# run the setup
setup(name = 'test_gsl',
	  version = '1.0',
	  packages = ["test_gsl"],
	  description = 'Testing GSL',
	  author_email = 'matthews@ast.cam.ac.uk',
	  author = 'James Matthews',
	  ext_modules=[alpro_module], 
	  py_modules=["test_gsl"],
	  package_data = {'test_gsl': data_files},
	  include_package_data=True)  
