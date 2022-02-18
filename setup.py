from distutils.core import setup, Extension
import numpy
import os

# data files 
data_files = ["data/*.dat"]

# run the setup
setup(name = 'alpro',
	  version = '1.0',
	  packages = ["alpro"],
	  description = 'Propagation of axion-like particles through magnetic fields',
	  author_email = 'matthews@ast.cam.ac.uk',
	  author = 'James Matthews',
	  py_modules=["alpro"],
	  package_data = {'alpro': data_files},
	  include_package_data=True)  
