Background theory and equations
--------------------------------------

The discussion here follows roughly the discussions given by 
`de Angelis et al 2011 <https://ui.adsabs.harvard.edu/abs/2011PhRvD..84j5030D/abstract>`_
and `Marsh et al 2017 <https://ui.adsabs.harvard.edu/abs/2017JCAP...12..036M/abstract>`_. The method
in the code most closely follows the former, but notation often matches the latter.

Alpro deals with relativistic ALPs, such that the relevant propagation equation (in the z direction)
is that of a reduced Schrodinger-like equation, which is first-order and given by

.. math::

	\left( i \frac{d}{d z} + \omega + {\cal M}_0 \right) \, \psi (z) = 0.

where :math:`\omega` is the beam energy, :math:`\psi (z) = (A_x, A_y, a)` is the beam photon-ALP state vector and :math:`{\cal M}_0` is the mixing matrix, given by

.. math::

	{\cal M}_0 = 
	\begin{bmatrix}
	\Delta_{xx} & \Delta_{xy} & \Delta^{x}_{ a \gamma} \\
	\Delta_{yx} & \Delta_{yy} & \Delta^{y}_{a \gamma} \\
	\Delta^{x}_{a \gamma} & \Delta^{y}_{\rm a \gamma} & \Delta_{a a} \\
	\end{bmatrix}

where we neglect the Faraday rotation terms :math:`\Delta_{yx}` and :math:`\Delta_{xy}`, and the other terms are given by (here quoted for the :math:`x` component)

.. math::

	\Delta^{x}_{ a \gamma} = \frac{g_a B_x}{2} \\
    \Delta_{a a} = -\frac{m_a^2}{2 \omega} \\
    \Delta_{xx} = \Delta_{p} = -\frac{\omega_p^2}{2\omega} \\

where :math:`m_a` and :math:`g_a` are the ALP mass and coupling constant, and are the fundamental ALP parameters in the problem. :math:`\omega_p` is the normal plasma frequency. The beam equation above resembles a Schrodinger equation with Hamiltonian :math:`-(\omega + {\cal M}_0)`. We now 
define the transfer matrix :math:`U_0(z, z_0)` which is the solution of the beam equation with initial condition :math:`U_0=1`. Then, the state vector at a given distance :math:`z` is given by

.. math::

	\psi (z) = U_0 (z, z_0) \, \psi (z_0) 

Idealised Case: B Aligned with y axis
===============================================
To obtain :math:`U_0`, we first consider the idealised case in which we choose the magnetic field to be aligned with the (we will generalise this result later). In this case the off-diagonal :math:`\Delta^{x}_{ a \gamma}` terms disappear and the mixing matrix becomes 

.. math::

	{\cal M}_0 = 
	\begin{bmatrix}
	\Delta_{p} & 0 & 0 \\
	0 & \Delta_{p} & \Delta^{y}_{a \gamma} \\
	0 & \Delta^{y}_{\rm a \gamma} & \Delta_{a a} \\
	\end{bmatrix}

This matrix can be diagonalised and it's Eigenvalues and Eigenvectors found (see e.g. Appendix A of `de Angelis et al 2011 <https://ui.adsabs.harvard.edu/abs/2011PhRvD..84j5030D/abstract>`_). The Eigenvalues are 

.. math::

	\lambda_1 &=& \Delta_{p} \\
    \lambda_2 &=& \frac{1}{2} (\Delta_{p} + \Delta_{a a} - \Delta_{\rm osc}) \\
    \lambda_3 &=& \frac{1}{2} (\Delta_{p} + \Delta_{a a} + \Delta_{\rm osc}) \\ 

where 

.. math::

    \Delta_{\rm osc} = \left[ (\Delta_{p} - \Delta_{a a})^2 + 4 (\Delta^{y}_{\rm a \gamma})^2 \right]^{1/2}

the Eigenvectors are most conveniently expressed in terms of an ALP mixing angle, :math:`\alpha`, defined as 

.. math::

    \alpha = \frac{1}{2} \rm{arctan} \left( \frac{2 \Delta^{y}_{\rm a \gamma}} {\Delta_p - \Delta_{a a}} \right)

which is a function of :math:`\omega, B, g, m_a` and :math:`\omega_{pl}`. In this case the Eigenvectors :math:`T_0, T_1, T_2` are written purely in terms of :math:`\sin \alpha` and  :math:`\cos \alpha` and resemble rotation matrices (they are given in the `Eigenvectors`_ subsection at the bottom of this document).

.. have a pure polarisation state initial condition :math:`\psi (z) = (1, 0, 0)` and we choose the magnetic field to be aligned with the (we will generalise this result later). 

General Case
===============================================
In the general case, the magnetic field is not aligned with either axis. If we define the angle :math:`\phi` as the angle the magentic field :math:`\vec{B}` makes with the y axis, then the general mixing matrix is a rotated version of the idealised mixing matrix, i.e. 

.. math::

	{\cal M}_0 = {\cal R}^T (\phi) {\cal M}^0_0 {\cal R} (\phi)

where 

.. math::

	{\cal R} (\phi) &=
	\begin{bmatrix}
	\cos \phi & \sin \phi & 0 \\
	\sin \phi & \cos \phi & 0 \\
	0 & 0 & 0
	\end{bmatrix}\\

Similarly, the generalised transfer matrix $U_0$ is given by 

.. math::

	{\cal U}_0 = {\cal R}^T (\phi) {\cal U}^0_0 {\cal R} (\phi) 



In this way, the propagation problem can be described completely by just the two angles :math:`\alpha` and :math:`\phi` and solved using a series of matrix operations. 


Code operation
=====================
The code follows more or less exactly the above formalism. 




Eigenvectors
=====================

.. \begin{equation}
.. \label{mravvq2abcappZX1}
.. {\cal U}_0 (y, y_0; 0) = e^{i {\lambda}_1 (y - y_0)} \, T_{0,1} (0) + e^{i {\lambda}_2 (y - y_0)} \, T_{0,2} (0) + e^{i {\lambda}_3 (y - y_0)} \, T_{0,3} (0)~, 
.. \end{equation}
.. where the matrices $T_{0,1} (0)$, $T_{0,2} (0) $ and $T_{0,3} (0)$ are just those defined by Eqs. (\ref{mravvq1app}), (\ref{mravvq2app}) and (\ref{mravvq3app}) as specialized to the present situation. Actually, a simplification is brought about by introducing the photon-ALP mixing angle
.. \begin{equation}
.. \label{a15}
.. \alpha = \frac{1}{2} \, {\rm arctg} \left(\frac{2 \, \Delta_{a \gamma} }{\Delta_{\rm pl} - \Delta_{a a}} \right) = 
.. \frac{1}{2} \, {\rm arctg} \left[\left( \frac{B}{M} \right) \left(\frac{2 E}{m^2 - {\omega}_{\rm pl}^2} \right) \right]~,
.. \end{equation}
.. since then simple trigonometric manipulations allow us to express the above matrices in the simpler form


.. math::

	T_{0} &=
	\begin{bmatrix}
	1 & 0 & 0 \\
	0 & 0 & 0 \\
	0 & 0 & 0
	\end{bmatrix}\\
	T_{1} &=
	\begin{bmatrix}
	0 & 0 & 0 \\ 
	0 & \sin^2 \alpha & - \sin \alpha \cos \alpha \\
	0 & - \sin \alpha \cos \alpha & \cos^2 \alpha 
	\end{bmatrix}\\
	T_{2} &=
	\begin{bmatrix}
	0 & 0 & 0 \\
	0 & \cos^2 \alpha  & \sin \alpha \cos \alpha \\
	0 & \sin \alpha \cos \alpha  & \sin^2 \alpha 
	\end{bmatrix}

