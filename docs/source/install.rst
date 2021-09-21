Installation
-----------------------------

To install alpro, run

.. code:: bash

    python setup.py install

Depending on your system, you may need to run ``sudo python setup.py install`` or ``python setup.py install --user``. You will need a number of standard python modules and a working installation of numba.

I have tested with Python 3.7 and later. Module requirements are relatively light, essentially ```matplotlib,numba,numpy,scipy```. Full details can be found in the requirements.txt file that ships with the repository. A minimal running environment can be loaded by starting a clean virtual environment with, e.g., ```python3 -m venv /path/to/virtual/env/``` and then running ```pip install -r requirements.txt``` before installing alpro. 

Running ALPro
====================================

To check ALPro has installed, try importing it as a module 

.. code:: bash

    python -c "import alpro"


