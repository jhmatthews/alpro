Tests
--------------

Alpro has been well-tested in a range of regimes against a few different analytical and numerical implementations. A basic test of alpro can be run using 

.. code:: bash

    python -m unittest alpro.test

if all goes well you should see output like 

	Ran 7 tests in 8.385s

	OK

This test will also create a figure with a filename like `Test_16-02-2022.png`, which corresponds to the date of the test. 

The notebooks below provide further details of some these tests with descriptions and instructions on how to run them individually. 

.. toctree::
   :glob:

   tests/*
