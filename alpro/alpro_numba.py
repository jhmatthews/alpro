import numpy as np
import scipy.linalg as linalg
from numba import jit
UNIT_LENGTH = 50676.79373667135
UNIT_GAUSS = 0.01953548032
HBAR_EV = 6.582119569e-16
MELEC = 9.10938356e-28
PARSEC = 3.0857e18
E = 4.8032045057134676e-10


@jit(nopython=True)
def get_P_matrix(energies, initial_state, phi, B, L, g_a, mass, ne, cmb_qed_term=True):

    # convert things to natural units
    M = 1.0 / g_a
    omega_pl = np.sqrt(4.0 * np.pi * E * E * ne / MELEC) * HBAR_EV
    distance = L * 1000.0 * PARSEC * UNIT_LENGTH
    B *= UNIT_GAUSS

    # array for storing probabilities
    Pnew = np.zeros_like(energies)

    Deltas = get_deltas(mass, energies, M, B, omega_pl, cmb_qed_term=cmb_qed_term)
    EVarray = get_eigenvalues(Deltas)
    alpha = 0.5 * np.arctan2(2.0 * Deltas[2], (Deltas[0] - Deltas[1]))

    # set up the initial state
    # get pure polarization vectors
    polx, poly, pola = get_pure_pol_vectors()
    pol = initial_state

    for i in range(len(energies)):
        energy = energies[i]
        one_pol = PropagateOnePolMatrix(
            alpha[i], EVarray[i], phi[i], pol[i, :, :], distance, energy)
        pol[i] = one_pol
        Pnew[i] = np.real(np.sum(np.diag(np.dot(pola, one_pol))))

    # what we return depends on the mode and is consistent with the
    # supplied initial_state

    new_state = pol

    return (Pnew, new_state)


@jit(nopython=True)
def get_P(energies, initial_state, phi, B, L, g_a, mass, ne, cmb_qed_term=True):

    # convert things to natural units
    M = 1.0 / g_a
    omega_pl = np.sqrt(4.0 * np.pi * E * E * ne / MELEC) * HBAR_EV
    distance = L * 1000.0 * PARSEC * UNIT_LENGTH
    B *= UNIT_GAUSS

    # array for storing probabilities
    Pnew = np.zeros_like(energies)

    Deltas = get_deltas(mass, energies, M, B, omega_pl, cmb_qed_term=cmb_qed_term)
    EVarray = get_eigenvalues(Deltas)
    alpha = 0.5 * np.arctan2(2.0 * Deltas[2], (Deltas[0] - Deltas[1]))

    # set up the initial state
    # convert to complex form if necessary.
    # make sure we return the same format
    # if initial_state.dtype == np.float:
    #     return_float = True
    #     initial_state = initial_state + 0.0j
    # else:
    #     return_float = False

    Ainit = initial_state
    A_new = np.zeros_like(Ainit)

    for i in range(len(energies)):
        energy = energies[i]

        # this will be a 3-vector in complex form
        A = Ainit[i]
        A2 = PropagateOne(alpha[i], EVarray[i], phi[i], A, distance, energy)

        A_new[i, :] = A2
        axion_amp = A2[2]
        Pnew[i] = np.abs(axion_amp) ** 2

    # what we return depends on the mode and is consistent with the
    # supplied initial_state
    # if return_float:
    #     new_state = A_new.real
    # else:
    new_state = A_new

    return (Pnew, new_state)


@jit(nopython=True)
def get_pure_pol_vectors():
    x = np.zeros((3, 3), dtype=np.complex128)
    y = np.zeros((3, 3), dtype=np.complex128)
    a = np.zeros((3, 3), dtype=np.complex128)
    x[0, 0] = 1.0 + 0.0j
    y[1, 1] = 1.0 + 0.0j
    a[2, 2] = 1.0 + 0.0j
    return (x, y, a)


@jit(nopython=True)
def PropagateOnePolMatrix(alpha, EVarray, phi, pol, distance, energy):
    # get deltas and eigenvalues
    #Deltas = get_deltas (mass, energy, M, B, omega_pl)
    #EVarray = get_eigenvalues (Deltas)

    # calculate T matrices from mixing angle, alpha (eq 3)
    #alpha = 0.5 * np.arctan2 (2.0 * Deltas["AG"], (Deltas["PL"] - Deltas["AA"]))
    T1, T2, T3 = get_T_matrices(alpha)

    #distance = distance
    # construct the transfer matrix for the idealised parallel situation */
    U0 = get_U0(EVarray, T1, T2, T3, distance)

    # # apply the rotation matrix
    U0 = apply_rotation_matrix(phi, U0)

    # # multiply the state vector (A) by the transfer matrix (U0)
    # # result is stored in A_new
    exp_term = 1.0 * np.exp(energy * 1j * distance)
    U0 = exp_term * U0
    #A_new = T1

    pol_new = np.dot(U0, np.dot(pol, U0.transpose().conjugate()))
    return (pol_new)


@jit(nopython=True)
def PropagateOne(alpha, EVarray, phi, A, distance, energy):

    # get deltas and eigenvalues
    #Deltas = get_deltas (mass, energy, M, B, omega_pl)
    #EVarray = get_eigenvalues (Deltas)

    # calculate T matrices from mixing angle, alpha (eq 3)
    #alpha = 0.5 * np.arctan2 (2.0 * Deltas["AG"], (Deltas["PL"] - Deltas["AA"]))
    T1, T2, T3 = get_T_matrices(alpha)

    #distance = distance
    # construct the transfer matrix for the idealised parallel situation */
    U0 = get_U0(EVarray, T1, T2, T3, distance)

    # # apply the rotation matrix
    U0 = apply_rotation_matrix(phi, U0)

    # # multiply the state vector (A) by the transfer matrix (U0)
    # # result is stored in A_new
    exp_term = 1.0 * np.exp(energy * 1j * distance)
    U0 = exp_term * U0
    #A_new = T1

    A_new = np.dot(A, U0)
    return (A_new)


@jit(nopython=True)
def apply_rotation_matrix(phi, U0):
    """
    Applies a 3x3 rotation matrix to the input matrix `U0`.

    Parameters:
    -----------
    phi : float
        Rotation angle in radians.
    U0 : np.ndarray
        A 3x3 complex matrix to which the rotation will be applied.

    Returns:
    --------
    np.ndarray
        A 3x3 complex matrix resulting from the rotation of `U0`.
    """
    nrows = 3
    v = np.zeros((nrows, nrows), dtype=np.complex128)

    # populate elements of rotation matrix
    v[0, 0] = np.cos(phi) + 0j
    v[0, 1] = -np.sin(phi) + 0j
    v[1, 0] = np.sin(phi) + 0j
    v[1, 1] = np.cos(phi) + 0j
    v[2, 2] = 1.0 + 0j
    vt = np.transpose(v)

    answer1 = np.dot(U0, v)
    answer = np.dot(vt, answer1)

    return (answer)


@jit(nopython=True)
def get_T_matrices(alpha):
    """
    Generates three specific 3x3 transformation matrices based on an input angle `alpha`. 
    These are the T matrices in de Angelis+ 2011. 

    Parameters:
    -----------
    alpha : float
        Angle in radians used to compute the transformation matrices.

    Returns:
    --------
    tuple of np.ndarray
        Three 3x3 complex transformation matrices (T1, T2, T3). 
    """
    nrows = 3
    T1 = np.zeros((nrows, nrows), dtype=np.complex128)
    T2 = np.zeros((nrows, nrows), dtype=np.complex128)
    T3 = np.zeros((nrows, nrows), dtype=np.complex128)

    # first matrix, eq 35
    T1[0, 0] = 1.0 + 0j

    # second and third matrixes, eqs 36 and 37
    T2[1, 1] = T3[2, 2] = np.sin(alpha) * np.sin(alpha) + 0j
    T2[2, 1] = T2[1, 2] = -np.sin(alpha) * np.cos(alpha) + 0j
    T3[2, 1] = T3[1, 2] = np.sin(alpha) * np.cos(alpha) + 0j
    T2[2, 2] = T3[1, 1] = np.cos(alpha) * np.cos(alpha) + 0j

    return (T1, T2, T3)


@jit(nopython=True)
def get_deltas(mass, energy, M, B, omega_pl, cmb_qed_term=True):
    """
    Computes the Delta terms for photon-ALP oscillation in the presence of a magnetic field and plasma.

    Parameters:
    -----------
    mass : float
        Particle mass.
    energy : np.ndarray
        Array of particle energy values.
    M : float
        Coupling constant.
    B : float
        Magnetic field strength.
    omega_pl : float
        Plasma frequency.
    cmb_qed_term : bool 
        whether to include the CMB and QED term or not 

    Returns:
    --------
    np.ndarray
        A 4xN array of Delta terms for each energy value, where N is the length of `energy`.
    """

    Deltas = np.zeros((4, len(energy)))

    # account for QED vacuum polarization and CMB photon dispersion
    Deltas[0] = (-(omega_pl * omega_pl) / 2.0 / energy) 

    if cmb_qed_term:
        Bcrit = 4.414e13 * UNIT_GAUSS
        w2_term = ((1.42e-4 * (B/Bcrit)**2) + (0.522e-42)) * energy 
        Deltas[0] += w2_term

    Deltas[1] = -(mass * mass) / 2.0 / energy
    Deltas[2] = B / 2.0 / M

    # now get Delta_osc
    x = (Deltas[0] - Deltas[1]) * (Deltas[0] - Deltas[1])
    Deltas[3] = np.sqrt(x + (4.0 * Deltas[2] * Deltas[2]))
    #print ("DELTAS", Deltas)
    return (Deltas)


@jit(nopython=True)
def get_eigenvalues(Deltas):
    """
    Calculates eigenvalues of the mixing matrix using the provided Delta terms.

    Parameters:
    -----------
    Deltas : np.ndarray
        A 4xN array of Delta terms where N is the number of energy values.

    Returns:
    --------
    np.ndarray
        An Nx3 array of eigenvalues for each set of Delta terms, where N is the length of Deltas[0].
    """
    EVarray = np.zeros((3, len(Deltas[0])))
    EVarray[0, :] = Deltas[0]
    EVarray[1, :] = 0.5 * (Deltas[0] + Deltas[1] - Deltas[3])
    EVarray[2, :] = 0.5 * (Deltas[0] + Deltas[1] + Deltas[3])
    return (EVarray.T)


@jit(nopython=True)
def get_U0(EVarray, T1, T2, T3, distance):
    """
    Computes the mixing matrix U0 based on eigenvalues and transformation matrices.

    Parameters:
    -----------
    EVarray : np.ndarray
        Array of eigenvalues with shape (N, 3), where N is the number of energy values.
    T1, T2, T2 : np.ndarray
        The 3x3 complex transformation matrices.
    distance : float
        Distance over which oscillation effects are calculated.

    Returns:
    --------
    np.ndarray
        A 3x3 complex matrix `U0` resulting from the eigenvalue and transformation matrix combination.
    """
    # get complex coefficients (e^i lambda_j d) */
    EVexp1 = 1.0 * np.exp(EVarray[0] * 1j * distance)
    EVexp2 = 1.0 * np.exp(EVarray[1] * 1j * distance)
    EVexp3 = 1.0 * np.exp(EVarray[2] * 1j * distance)

    # multiply complex coefficients by the matrices and add together */
    U0 = (EVexp1 * T1) + (EVexp2 * T2) + (EVexp3 * T3)
    # /* add the terms together and copy into the new matrix */

    return U0
