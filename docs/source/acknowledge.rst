Acknowledging ALPro
-------------------------------

If you use ALPro in your work, please acknowledge it and cite both the `Zenodo DOI <https://doi.org/10.5281/zenodo.6079445>`_ and the following paper:

	`J H Matthews et al 2022 </#>`_
	“How do magnetic field models affect astrophysical limits on light axion-like particles? An X-ray case study with NGC 1275”, Accepted to ApJ. 

Bibtex for the Zenodo DOI can be accessed in python via 

.. code:: python

     import alpro; print (alpro.__citation__)

which prints out::

    @misc{alpro,
	author       = {James Matthews},
	title        = {{alpro: Axion-Like PROpagation}},
	month        = feb,
	year         = 2022,
	doi          = {110.5281/zenodo.6079445},
	version      = { 1.0 },
	publisher    = {Zenodo},
	url          = {https://doi.org/10.5281/zenodo.6079445}
	}


If you use the `pruning` scheme (an adaptive treatment for conservatively dealing with ALP resonances), please cite the paper where it is described:

	`J Sisk Reynes et al 2021 <https://ui.adsabs.harvard.edu/abs/2022MNRAS.510.1264S/abstract>`_ 
	"New constraints on light axion-like particles using Chandra transmission grating spectroscopy of the powerful cluster-hosted quasar H1821+643", MNRAS, 510, 64

