from .alpro_numba import *
from .main import *
import alpro.util as util
import alpro.test as test
import alpro.pure as pure
import alpro.fourier as fourier
import alpro.fourier_core as fourier_core
import importlib

# test for mpi4py before including the parallel section of alpro
found_mpi = True 
try:
  import mpi4py
except ModuleNotFoundError:
  found_mpi = False 

if found_mpi:
    import alpro.parallel as parallel

__author__ = "James Matthews"
__email__ = "matthews@ast.cam.ac.uk"
__description__ = "Propagation of axion-like particles through magnetic fields"
__version__ = "1.1"

# Paper citation
__citation__ = __bibtex__ = """@ARTICLE{Matthews2022,
       author = {{Matthews}, James H. and {Reynolds}, Christopher S. and {Marsh}, M.~C. David and {Sisk-Reyn{\'e}s}, J{\'u}lia and {Rodman}, Payton E.},
        title = "{How do Magnetic Field Models Affect Astrophysical Limits on Light Axion-like Particles? An X-ray Case Study with NGC 1275}",
      journal = {arXiv e-prints},
     keywords = {Astrophysics - High Energy Astrophysical Phenomena, Astrophysics - Cosmology and Nongalactic Astrophysics, High Energy Physics - Phenomenology},
         year = 2022,
        month = feb,
          eid = {arXiv:2202.08875},
        pages = {arXiv:2202.08875},
archivePrefix = {arXiv},
       eprint = {2202.08875},
 primaryClass = {astro-ph.HE},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2022arXiv220208875M},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}"""

# Zenodo DOI citation - note that this cites "all versions" 
__zenodo__ = """@software{{alpro,
author       = {{James Matthews}},
title        = {{{{alpro: Axion-Like PROpagation}}}},
month        = feb,
year         = 2022,
doi          = {{10.5281/zenodo.6137185}},
version      = {{ {} }},
publisher    = {{Zenodo}},
url          = {{https://doi.org/10.5281/zenodo.6137185}}
}}""".format(__version__)

