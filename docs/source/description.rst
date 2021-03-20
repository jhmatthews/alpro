Code Description and Method
--------------------------------------

The discussion here follows roughly the discussions given by 
`de Angelis et al 2011 <https://ui.adsabs.harvard.edu/abs/2011PhRvD..84j5030D/abstract>`_
and `Marsh et al 2017 <https://ui.adsabs.harvard.edu/abs/2017JCAP...12..036M/abstract>`_. The method
in the code most closely follows the former, but notation often matches the latter.

Background theory and equations
====================================
Alpro deals with relativistic ALPs, such that the relevant propagation equation (in the z direction)
is that of a reduced Schrodinger-like equation, which is first-order and given by

.. math::

	\left( i \frac{d}{d z} + \omega + {\cal M}_0 \right) \, \psi (z) = 0.

where :math:`\omega` is the beam energy, :math:`\psi (z)` is the beam photon-ALP state vector and :math:`{\cal M}_0` is the mixing matrix, given by

.. math::

	{\cal M}_0 = 
	\begin{bmatrix}
	\Delta_{xx} & \Delta_{xy} & \Delta^{x}_{ a \gamma} \\
	\Delta_{yx} & \Delta_{yy} & \Delta^{y}_{a \gamma} \\
	\Delta^{x}_{a \gamma} & \Delta^{y}_{\rm a \gamma} & \Delta_{a a} \\
	\end{bmatrix}

where we neglect the Faraday rotation terms :math:`\Delta_{yx}` and :math:`\Delta_{xy}`, and the other terms are given by (for the :math:`x` component)

.. math::

	\Delta^{x}_{ a \gamma} = \frac{g_a B_x}{2} \\
    \Delta_{a a} = -\frac{m_a^2}{2 \omega} \\
    \Delta_{xx} = -\frac{\omega_p^2}{2\omega} \\

where :math:`m_a` and :math:`g_a` are the ALP mass and coupling constant, and are the fundamental ALP parameters in the problem. :math:`\omega_p` is the normal plasma frequency. The beam equation above resembles a Schrodinger equation with Hamiltonian :math:`-(\omega + {\cal M}_0)`. 

Code operation
=====================
The code follows more or less exactly the above formalism. 
