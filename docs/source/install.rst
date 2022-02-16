Installation
-----------------------------

To install alpro, run

.. code:: bash

    python setup.py install

Depending on your system, you may need to run ``sudo python setup.py install`` or ``python setup.py install --user``. You will need a number of standard python modules and a working installation of numba.

I have tested with Python 3.7 and later. Module requirements are relatively light, essentially ``matplotlib,numba,numpy,scipy``. Full details can be found in the requirements.txt file that ships with the repository. A minimal running environment can be loaded by starting a clean virtual environment with, e.g., 

.. code:: bash

    python -m venv /path/to/virtual/env/ 

and then running 

.. code:: bash

    pip install -r requirements.txt

before installing alpro. 

Running ALPro
====================================

To check ALPro has installed, try importing it as a module 

.. code:: bash

    python -c "import alpro"

A basic test of the code can be run using 

.. code:: bash

    python -m unittest alpro.test

if all goes well you should see output like 

	Ran 4 tests in 4.447s

	OK

Parallelisation 
====================================

alpro includes a parallelisation routine. This is only imported in the ``__init__.py`` file if the `mpi4py module <https://mpi4py.readthedocs.io/en/stable/>`_ is found, otherwise it is assumed the user does not wish to use this routine. If ``mpi4py`` was present at the time of installation, the parallel routine can be tested with 

.. code:: bash

    mpirun -n 4 python -m unittest alpro.parallel

and the user can check the speedup by changing the integer after ``-n``. If this test produces an ``ImportError``, try installing ``mpi4py`` and then reinstalling alpro from scratch. Instructions on using the parallel routines are given in :ref:`Examples <examples>`.

