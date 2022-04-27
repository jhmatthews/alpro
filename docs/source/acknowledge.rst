Acknowledging ALPro
-------------------------------

If you use ALPro in your work, please acknowledge it and cite both the `Zenodo DOI <https://doi.org/10.5281/zenodo.6079445>`_ and the following paper:

	`J H Matthews et al 2022 </#>`_
	“How do magnetic field models affect astrophysical limits on light axion-like particles? An X-ray case study with NGC 1275”, Accepted to ApJ. 

Bibtex for the paper can be accessed in python via 

.. code:: python

     import alpro; print (alpro.__citation__)

which prints out::

	@ARTICLE{Matthews2022,
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
	}



Bibtex for the Zenodo DOI can be accessed in python via 

.. code:: python

     print (alpro.__zenodo__)

which prints out::

	@software{alpro,
		author       = {James Matthews},
		title        = {{{{alpro: Axion-Like PROpagation}}}},
		month        = feb,
		year         = 2022,
		publisher    = {Zenodo},
		version      = {v1.1},
		doi          = {10.5281/zenodo.6137185},
		url          = {https://doi.org/10.5281/zenodo.6137185}
	}


If you use the `pruning` scheme (an adaptive treatment for conservatively dealing with ALP resonances), please cite the paper where it is described:

	`J Sisk Reynes et al 2021 <https://ui.adsabs.harvard.edu/abs/2022MNRAS.510.1264S/abstract>`_ 
	"New constraints on light axion-like particles using Chandra transmission grating spectroscopy of the powerful cluster-hosted quasar H1821+643", MNRAS, 510, 64

