import numpy as np 
import scipy.linalg as linalg
from numba import jit
UNIT_LENGTH= 50676.79373667135
UNIT_GAUSS = 0.01953548032
HBAR_EV = 6.582119569e-16
MELEC = 9.10938356e-28
PARSEC = 3.0857e18
E = 4.8032045057134676e-10

@jit(nopython=True)      
def get_P2(energies, Ainit, phi, B, L, g_a, mass, ne):
    M = 1.0 / g_a 
    omega_pl = np.sqrt(4.0 * np.pi * E * E * ne / MELEC) * HBAR_EV
    distance = L * 1000.0 * PARSEC * UNIT_LENGTH
    B *= UNIT_GAUSS
    A_new = np.zeros_like(Ainit)
    Pnew = np.zeros_like(energies)

    Deltas = get_deltas (mass, energies, M, B, omega_pl)
    EVarray = get_eigenvalues (Deltas)
    alpha = 0.5 * np.arctan2 (2.0 * Deltas[2], (Deltas[0] - Deltas[1]))

    # U0_all = []

    for i in range(len(energies)):
        Atemp = Ainit[i].reshape((3,2))

        # write in complex form 
        energy = energies[i]
        A = Atemp[:,0] + 1j * Atemp[:,1]
        A2 = PropagateOne(alpha[i], EVarray[i], phi[i], A, distance, energy)

        A_new[i,::2] = A2.real
        A_new[i,1::2] = A2.imag
        #print (A2.shape)
        mod = A2[2]
        Pnew[i] = np.abs(mod) **2

    #Pnew = 1
    return (Pnew, A_new)

@jit(nopython=True) 
def PropagateOne (alpha, EVarray, phi, A, distance, energy):
    # get deltas and eigenvalues 
    #Deltas = get_deltas (mass, energy, M, B, omega_pl)
    #EVarray = get_eigenvalues (Deltas)

    # calculate T matrices from mixing angle, alpha (eq 3)
    #alpha = 0.5 * np.arctan2 (2.0 * Deltas["AG"], (Deltas["PL"] - Deltas["AA"]))
    T1, T2, T3 = get_T_matrices (alpha)

    #distance = distance 
    # construct the transfer matrix for the idealised parallel situation */
    U0 = get_U0 (EVarray, T1, T2, T3, distance)

    # # apply the rotation matrix 
    U0 = apply_rotation_matrix (phi, U0)

    # # multiply the state vector (A) by the transfer matrix (U0) 
    # # result is stored in A_new 
    exp_term = 1.0 * np.exp(energy * 1j * distance)
    U0 = exp_term * U0
    #A_new = T1

    A_new = np.dot(A,U0)
    return (A_new)

@jit(nopython=True) 
def apply_rotation_matrix (phi, U0):
    nrows = 3
    v = np.zeros( (nrows,nrows), dtype=np.complex128)

    # populate elements of rotation matrix 
    v[0,0] = np.cos (phi) + 0j
    v[0,1] = -np.sin (phi) + 0j
    v[1,0] = np.sin (phi) + 0j
    v[1,1] = np.cos (phi) + 0j
    v[2,2] = 1.0+0j
    vt = v.T

    answer = np.dot(U0,v)
    answer = np.dot(vt,U0)

    return (answer)

@jit(nopython=True) 
def get_T_matrices (alpha):
    nrows = 3
    T1 = np.zeros((nrows,nrows), dtype=np.complex128)
    T2 = np.zeros((nrows,nrows), dtype=np.complex128)
    T3 = np.zeros((nrows,nrows), dtype=np.complex128)

    # first matrix, eq 35 
    T1[0,0] = 1.0 + 0j

    # second and third matrixes, eqs 36 and 37 
    T2[1,1] = T3[2,2] = np.sin (alpha) * np.sin (alpha) + 0j
    T2[2,1] = T2[1,2] = -np.sin (alpha) * np.cos (alpha) + 0j
    T3[2,1] = T3[1,2] = np.sin (alpha) * np.cos (alpha) + 0j
    T2[2,2] = T3[1,1] = np.cos (alpha) * np.cos (alpha) + 0j

    return (T1, T2, T3)

@jit(nopython=True) 
def get_deltas (mass, energy, M, B, omega_pl):
    Deltas = np.zeros((4,len(energy)))
    Deltas[0] = -(omega_pl * omega_pl) / 2.0 / energy
    Deltas[1] = -(mass * mass) / 2.0 / energy
    Deltas[2] = B / 2.0 / M 

    # now get Delta_osc 
    x = (Deltas[0] - Deltas[1]) * (Deltas[0] - Deltas[1])
    Deltas[3] = np.sqrt (x + (4.0 * Deltas[2] * Deltas[2]))
    return (Deltas)

@jit(nopython=True) 
def get_eigenvalues (Deltas):
    # get eigenvalues of mixing matrix
    EVarray = np.zeros( (3,len(Deltas[0])))
    EVarray[0,:] = Deltas[0] 
    EVarray[1,:] = 0.5 * (Deltas[0] + Deltas[1] - Deltas[3])
    EVarray[2,:] = 0.5 * (Deltas[0] + Deltas[1] + Deltas[3])
    return (EVarray.T)


@jit(nopython=True) 
def get_U0 (EVarray, T1, T2, T3, distance):
    # get complex coefficients (e^i lambda_j d) */
    EVexp1 = 1.0 * np.exp(EVarray[0] * 1j * distance)
    EVexp2 = 1.0 * np.exp(EVarray[1] * 1j * distance)
    EVexp3 = 1.0 * np.exp(EVarray[2] * 1j * distance)

    # multiply complex coefficients by the matrices and add together */
    U0 = (EVexp1 * T1) + (EVexp2 * T2) + (EVexp3 * T3)
    #/* add the terms together and copy into the new matrix */

    return U0



