The Fourier Formalism
--------------------------------------

The Fourier formalism is described by `Marsh et al 2022 <https://ui.adsabs.harvard.edu/abs/2021arXiv210708040M/abstract>`_. 
Marsh et al. showed that, to leading order in :math:`g_{a\gamma}`, the conversion probability :math:`P_{\gamma_i\to a}` (with :math:`i\in [x,y]`) can be related to :math:`\Delta_i(z)`, or equivalently, the magnetic field profile along the line of sight $:math:`B_i(z)`, using Fourier-like transforms. Although this treatment breaks down when the conversion probabilities exceed 5-10 per cent, it is an extremely useful framework when considering how different magnetic field treatments affect the conversion probability. 

In the massive ALP (:math:`m_a \gg \omega_{\rm pl}`) case and focusing on the $x$-component only, the conversion probability can be written as 

.. math::

	P_{\gamma_x \to a}(\eta) = {\cal F}_s( \Delta_{x})^2 + {\cal F}_c( \Delta_{x})^2 

where 

.. math::

	{\cal F}_c = \int^\infty_0 f(z) \cos(\eta z) dz

and 

.. math::

	{\cal F}_s = \int^\infty_0 f(z) \sin(\eta z) dz

denote the cosine and sine transforms using a conjugate variable :math:`\eta=m_a^2/2E`. By applying the `Wiener Khinchin theorem <https://en.wikipedia.org/wiki/Wiener%E2%80%93Khinchin_theorem>`_, the conversion probability can also be expressed in terms of a cosine transform of the autocorrelation function of the line-of-sight magnetic field, :math:`c_{B_x}`, as 

.. math::

	P_{\gamma_x \to a}(\eta) = \frac{g_{a\gamma}^2}{2} {\cal F}_c \Big( c_{B_x}(L)\Big).

In the massless case (:math:`m_a \ll \omega_{\rm pl}`), the same formalism applies if we transform to new variables.  Specifically, :math:`\Delta_x` is replaced by the function :math:`G=2 \Delta_x/\omega_{\rm pl}^2`, the line of sight distance coordinate is replaced by a phase factor proportional to the electron column density

.. math::

    \varphi=\frac{1}{2}\int_0^z{\rm d}z^\prime\omega_{\rm pl}^2(z^\prime),

and the conjugate Fourier variable is :math:`\lambda=1/E`. With these transformations,
:math:`\Delta_x` is replaced by the function :math:`G=2 \Delta_x/\omega_{\rm pl}` with a coordinate change from :math:`z` to a phase :math:`\varphi`, and a conjugate variable  :math:`\lambda=1/E`. Nevertheless, the basic principle is similar to the massive ALP case, as the conversion probability can still be expressed as a simple transform of a function of the line of sight perpendicular magnetic field. Although the forms above are for a polarized beam, the unpolarised survival probability can always be obtained directly from 

.. math::

	P_{\gamma\gamma} = (1 - P_{\gamma a}) = 1 - \frac{1}{2} \left( P_{\gamma_{x}\rightarrow a} + P_{\gamma_{y}\rightarrow a} \right),

ALPro comes with routines for applying the Fourier formalism in the massless and massive regimes. More information can be found under :doc:`Examples & Tutorials <../examples>`.  