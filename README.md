## ALPro

[![Documentation Status](https://readthedocs.org/projects/alpro/badge/?version=latest)](https://alpro.readthedocs.io/en/latest/?badge=latest)
[![ALPro tests](https://github.com/jhmatthews/alpro/actions/workflows/alpro-tests.yml/badge.svg)](https://github.com/jhmatthews/alpro/actions/workflows/alpro-tests.yml)
[![DOI](https://zenodo.org/badge/218793537.svg)](https://zenodo.org/badge/latestdoi/218793537)
[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![PEP8](https://img.shields.io/badge/code%20style-pep8-orange.svg)](https://www.python.org/dev/peps/pep-0008/)

Alpro (Axion-Like PROpagation) is a code for solving the equation of motion for a photon-ALP beam with mixing between the photon and ALP states. 

### Installation

alpro can be installed by cloning the repository and running

```
python setup.py install
```

in the directory. You may need a sudo, or a --user flag, depending on your system privileges. 

### Prerequisites

I have tested with Python 3.7 and later. Module requirements are relatively light, essentially ```matplotlib,numba,numpy,scipy```. Full details can be found in the requirements.txt file. A minimal running environment can be loaded by starting a clean virtual environment with, e.g., ```python3 -m venv /path/to/virtual/env/``` and then running ```pip install -r requirements.txt``` before installing alpro. 

### Documentation

Documentation is hosted on [ReadTheDocs](https://alpro.readthedocs.io/en/latest/). 

### Release Paper and Contact

The code will be described in a paper, currently in preparation. If you are interested in using the code before the paper comes out, please contact me on matthews [at] ast.cam.ac.uk. 
