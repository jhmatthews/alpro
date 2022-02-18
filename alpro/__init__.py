from .alpro_numba import *
from .main import *
import alpro.util as util
import alpro.test as test
import alpro.pure as pure
import alpro.fourier as fourier
import alpro.fourier_core as fourier_core
import importlib
mpi_loader = importlib.find_loader('mpi4py')
found_mpi = mpi_loader is not None
if found_mpi:
    import alpro.parallel as parallel

__author__ = "James Matthews"
__email__ = "matthews@ast.cam.ac.uk"
__description__ = "Propagation of axion-like particles through magnetic fields"
__version__ = "1.1"
__citation__ = __bibtex__ = __zenodo__ = """@misc{{alpro,
author       = {{James Matthews}},
title        = {{{{alpro: Axion-Like PROpagation}}}},
month        = feb,
year         = 2022,
doi          = {{10.5281/zenodo.6079444}},
version      = {{ {} }},
publisher    = {{Zenodo}},
url          = {{https://doi.org/10.5281/zenodo.6137185}}
}}""".format(__version__)
