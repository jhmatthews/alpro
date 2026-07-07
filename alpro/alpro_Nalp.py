import numpy as np
from numba import jit
import alpro.models as models

# constants reused from the original module
UNIT_LENGTH = 50676.79373667135
UNIT_GAUSS = 0.01953548032
HBAR_EV = 6.582119569e-16
MELEC = 9.10938356e-28
PARSEC = 3.0857e18
E = 4.8032045057134676e-10

# ----------------------------------------------------------------------------------------------------------------------------- #
# Extending ALP propagation to N ALPs
# ----------------------------------------------------------------------------------------------------------------------------- #

@jit(nopython=True)
def apply_rotation_matrix_Nalp(phi, U0, N):
    """ Apply (N+2)x(N+2) rotation to input U0 """

    # construct empty (N+2)x(N+2) matrix
    R = np.zeros((N+2, N+2), dtype=np.complex128)
    c, s = np.cos(phi), np.sin(phi)

    # photon block
    R[0, 0] = c
    R[0, 1] = -s
    R[1, 0] = s
    R[1, 1] = c

    # ALP block
    for i in range(N):
        R[i+2, i+2] = 1.0

    # apply the rotation matrix to the input matrix U0
    U_rot = R.T @ U0 @ R

    return U_rot

@jit(nopython=True)
def construct_mixing_matrix_Nalp(masses, couplings, energy, B, omega_pl, cmb_qed_term=True):
    """ constructs the idealised (N+2)x(N+2) mixing matrix given input masses, couplings etc """

    N = len(masses)
    dim = N + 2

    # construct empty (N+2)x(N+2) matrix
    M0 = np.zeros((dim, dim), dtype=np.complex128)

    # photon-photon block
    delta_p = -omega_pl**2 / (2.0 * energy)
    M0[0, 0] = delta_p
    M0[1, 1] = delta_p

    # QED vacuum polarisation correction term
    if cmb_qed_term:
        B_crit = 4.414e13 * UNIT_GAUSS
        qed_term = ((1.42e-4 * (B/B_crit)**2) + 0.522e-42) * energy
        M0[0, 0] += qed_term
        M0[1, 1] += qed_term

    # ALP-ALP block
    for i in range(N):
        M0[i+2, i+2] = -masses[i]**2 / (2.0 * energy)

    # photon-ALP mixing block
    for i in range(N):
        delta_agi_y = couplings[i] * B / 2.0
        M0[1, i+2] = delta_agi_y
        M0[i+2, 1] = delta_agi_y

    return M0

@jit(nopython=True)
def get_U0_Nalp(masses, couplings, energy, B, omega_pl, distance, cmb_qed_term=True):
    """ calculates the idealised transfer matrix """

    # retrieve the idealised mixing matrix
    M0_ideal = construct_mixing_matrix_Nalp(masses, couplings, energy, B, omega_pl, cmb_qed_term=cmb_qed_term) 

    # find eigenvalues and eigenvectors for the idealised mixing matrix
    eigenvalues, eigenvectors = np.linalg.eigh(M0_ideal)

    # contruct the transfer matrix in the diagonal basis
    diag_U0 = np.diag(np.exp(1j * eigenvalues * distance))

    # transform back to the original basis
    U0 = eigenvectors @ diag_U0 @ eigenvectors.T

    return U0

@jit(nopython=True)
def PropagateOne_Nalp(masses, couplings, energy, B, omega_pl, phi, A, distance, cmb_qed_term=True):
    """ propagates state vector through single domain cell """

    N = len(masses)

    # retrieve the transfer matrix
    U0 = get_U0_Nalp(masses, couplings, energy, B, omega_pl, distance, cmb_qed_term=cmb_qed_term)

    # apply the rotation matrix for the general case
    U0 = apply_rotation_matrix_Nalp(phi, U0, N)

    # apply back the global phase factor
    phase = 1.0 * np.exp(energy * 1j * distance)
    U0 = phase * U0

    # propagate the state vector A using the transfer matrix U0
    A_new = np.dot(A, U0)

    return (A_new)

@jit(nopython=True)
def get_P_Nalp(energies, initial_state, phi, B, L, couplings, masses, ne, cmb_qed_term=True):

    # convert to natural units
    omega_pl = np.sqrt(4.0 * np.pi * E**2 * ne / MELEC) * HBAR_EV
    distance = L * 1000.0 * PARSEC * UNIT_LENGTH
    B *= UNIT_GAUSS

    # arrays for storing the results
    Pnew = np.zeros_like(energies, dtype=np.float64)
    Ainit = initial_state
    A_new = np.zeros_like(Ainit, dtype=np.complex128)

    for i in range(len(energies)):
        energy = energies[i]

        A = Ainit[i, :]
        A2 = PropagateOne_Nalp(masses, couplings, energy, B, omega_pl, phi[i], A, distance, cmb_qed_term=cmb_qed_term)

        A_new[i, :] = A2
        # now sum over all N ALP states to get conversion probability (i.e. a_1^2 + ... + a_n^2)
        Pnew[i] = np.sum(np.abs(A2[2:])**2)

    new_state = A_new
    return (Pnew, new_state)

def set_params_Nalp(self, couplings, masses):

    if len(couplings) != len(masses):
        raise ValueError("couplings and masses must have the same length")
    
    self.couplings = np.array(couplings)
    self.masses = np.array(masses)
    self.N_alps = len(masses)

def propagate_Nalp(self, domain=None, energies=None, pol="both"):

    if domain is None:
        domain = self.domain
    if energies is None:
        energies = self.energies

    dim = self.N_alps + 2

    # decide which polarization states to compute
    xpol, ypol = False, False
    if pol == "both":
        xpol, ypol = True, True
    elif pol == "x":
        xpol = True
    elif pol == "y":
        ypol = True
    else:
        raise ValueError("pol keyword must be 'x', 'y' or 'both'")
    
    # initialize the state vectors
    if ypol:
        init_y = np.zeros((len(energies), dim), dtype=np.complex128)
        init_y[:, 1] = 1.0
    if xpol:
        init_x = np.zeros((len(energies), dim), dtype=np.complex128)
        init_x[:, 0] = 1.0

    self.P_radial = np.zeros((len(domain.deltaL), len(energies)))

    # loop over the domain
    for i in range(len(domain.deltaL)):
        L = domain.deltaL[i]
        B = domain.B[i]
        phi = domain.phi[i] * np.ones_like(energies)
        ne = domain.ne[i]

        if ypol:
            P_y, new_y = get_P_Nalp(energies,init_y,phi,B,L,self.couplings,self.masses,ne)
            init_y = new_y
        if xpol:
            P_x, new_x = get_P_Nalp(energies,init_x,phi,B,L,self.couplings,self.masses,ne)
            init_x = new_x

        if xpol and ypol:
            Ptot = 0.5 * (P_y + P_x)
        elif xpol:
            Ptot = P_x
        elif ypol:
            Ptot = P_y

        self.P_radial[i, :] = Ptot

    self.P = Ptot
    self.energies = energies

    return (self.P, self.P_radial)

# ----------------------------------------------------------------------------------------------------------------------------- #
# Dealing with the issue of resonances in N ALP scenarios (turbulent domains)
# ----------------------------------------------------------------------------------------------------------------------------- #


def get_resonant_conditions(domain, energies, couplings, scale_factor):
    """ calculates conditions for a cell being resonant for each field parameter, energy, and ALP parameter.
    scale factor defines limit, i.e. 0.01 -> the 1% level """

    # get necessary quantities and convert to natural units
    domain_sizes = domain.deltaL
    L_natural_units = domain_sizes * models.unit.kpc * models.unit.unit_length_natural # eV^-1
    B_perp = np.sqrt(domain.Bx**2 + domain.By**2) # in units of Gauss
    B_perp_natural = B_perp * models.unit.unit_gauss_natural # Gauss to eV^2

    # 1. ensure Delta is small 
    energy_thresholds = scale_factor * energies[:, np.newaxis] / L_natural_units[np.newaxis, :]

    # 2. ensure Theta is large 
    magnetic_thresholds_per_alp = scale_factor * np.abs(couplings)[:, np.newaxis, np.newaxis] * energies[np.newaxis, :, np.newaxis] * B_perp_natural[np.newaxis, np.newaxis, :]
    
    # 3. ensure m_eff is small due to a real resonance and not simply due to small m_a and ω_p 
    plasma_thresholds = scale_factor * domain.omega_p**2

    # 4. ensure Theta * Delta < 1
    angle_condition_per_alp = 0.5 * np.abs(couplings[:, np.newaxis]) * (B_perp_natural * L_natural_units)[np.newaxis, :]

    return energy_thresholds, magnetic_thresholds_per_alp, plasma_thresholds, angle_condition_per_alp

def find_resonances(domain,energies,masses,couplings,scale_factor):
    """ finds which cells are resonant in the domain.
    A boolean is used to track resonances. True = Resonant, False = not Resonant.
    Cell must satisfy all conditions to be considered resonant. 
    Depends on the ALP parameters so each ALP has their own profile."""

    # get energy-dependent resonance conditions
    energy_thresholds, magnetic_thresholds_per_alp, plasma_thresholds, angle_condition_per_alp = get_resonant_conditions(domain,energies,couplings,scale_factor)
    per_alp_diagnostics = {}

    # for each ALP find the resonances
    for m in range(len(masses)):
        mass_name = f"ALP_{m+1}"
        m_eff_sq = np.abs(masses[m]**2 - domain.omega_p**2)
        angle_condition = angle_condition_per_alp[m, :] < 1.0
        plasma_condition = m_eff_sq < plasma_thresholds

        # find cells satisfying all conditions
        alp_resonant_mask = np.zeros(len(domain.deltaL),dtype=bool)
        for e in range(len(energies)):
            energy_condition = m_eff_sq < energy_thresholds[e, :]
            magnetic_condition = m_eff_sq < magnetic_thresholds_per_alp[m, e, :]
            combined_condition = energy_condition & magnetic_condition & plasma_condition & angle_condition
            alp_resonant_mask |= combined_condition
        alp_resonant_indices = np.where(alp_resonant_mask)[0]

        per_alp_diagnostics[mass_name] = {'mass':masses[m].item(),'coupling':couplings[m].item(),'n_resonant_cells':len(alp_resonant_indices),'resonant_indices': alp_resonant_indices.tolist()}
        
    return per_alp_diagnostics

def classify_resonances(domain,energies,masses,couplings,scale_factor,verbose=False):
    """ classifies resonances as 'isolated' or 'clustered'. 
    'isolated' denotes resonant cells with neighbouring cells that are NOT resonant.
    'clustered' denotes resonant cells with neibouring cells that ARE also resonant."""
     
    # fetch resonances from find_resonances
    per_alp_diagnostics = find_resonances(domain,energies,masses,couplings,scale_factor)
    per_alp_clustering = {}
    n_cells = len(domain.deltaL)

    for m in range(len(masses)):

        alp_resonant_indices = np.array(per_alp_diagnostics[f"ALP_{m+1}"]['resonant_indices'],dtype=int)
        alp_resonant_indices_bool = np.zeros(n_cells, dtype=bool)
        alp_resonant_indices_bool[alp_resonant_indices] = True
        alp_n_resonances = len(alp_resonant_indices)

        # classify resonant cells for this ALP and store indices 
        # done by checking values in adjacent cells of the resonance tracking boolean array.
        has_adjacent = np.zeros(alp_n_resonances, dtype=bool)
        for i, idx in enumerate(alp_resonant_indices):
            if idx > 0 and alp_resonant_indices_bool[idx - 1]:
                has_adjacent[i] = True
            if idx < n_cells - 1 and alp_resonant_indices_bool[idx + 1]:
                has_adjacent[i] = True

        clustered_cell_indices = alp_resonant_indices[has_adjacent]
        isolated_cell_indices = alp_resonant_indices[~has_adjacent]
        mass_name = f"ALP_{m+1}"
        per_alp_clustering[mass_name] = {'n_clustered':len(clustered_cell_indices),'n_isolated':len(isolated_cell_indices),'clustered_cell_indices':clustered_cell_indices.tolist(),'isolated_cell_indices':isolated_cell_indices.tolist()}
        
    # print summary diagnostics
    if verbose:
        for m in range(len(masses)):
            mass_name = f"ALP_{m+1}"
            mass_value = per_alp_diagnostics[mass_name]['mass']
            coupling = per_alp_diagnostics[mass_name]['coupling']
            n_resonances = per_alp_diagnostics[mass_name]['n_resonant_cells']
            n_clustered = per_alp_clustering[mass_name]['n_clustered']
            n_isolated = per_alp_clustering[mass_name]['n_isolated']
            resonant_indices = per_alp_diagnostics[mass_name]['resonant_indices']
            print(f"{mass_name}, mass: {mass_value}, coupling: {coupling}, n_resonances: {n_resonances}, n_isolated: {n_isolated}, n_clustered: {n_clustered}, Resonant Indices: {np.array(sorted(resonant_indices))}")
        
    return per_alp_diagnostics, per_alp_clustering

def refine_domain(domain,cluster,energies,masses,couplings,scale_factor,res_class):
    """ subdivides resonant cells in the domain. Gives new cells appropriate parametes.
    User can choose to subdivide all cells, only 'isolated' or only 'clustered' cells."""

    if res_class not in ["isolated", "clustered", "all"]:
        raise ValueError("res_class must be 'isolated', 'clustered', or 'all'")

    # 1. FINDS RESONANT CELLS IN THE DOMAIN

    # store initial cell count
    initial_cell_count = len(domain.deltaL)

    # find cells to refine
    _, per_alp_clustering = classify_resonances(domain,energies,masses,couplings,scale_factor)

    resonances_to_refine = set()
    for m in range(len(masses)):

        # Add isolated cells 
        if res_class in ["isolated", "all"]:
            resonances_to_refine.update(per_alp_clustering[f"ALP_{m+1}"]['isolated_cell_indices'])
            
        # Add clustered cells
        if res_class in ["clustered", "all"]:
            resonances_to_refine.update(per_alp_clustering[f"ALP_{m+1}"]['clustered_cell_indices'])
            
    cells_to_refine = sorted(list(resonances_to_refine))
    
    # if no cells to refine, return early
    if len(cells_to_refine) == 0:
        return domain
    
    # 2. SUBDIVIDES THOSE RESONANT CELLS

    # begin the cell subdivision process (split resonant cells into two subcells)
    cell_mask = np.zeros(len(domain.deltaL), dtype=bool)
    cell_mask[cells_to_refine] = True

    # create new arrays for the refined domain
    new_cell_count = initial_cell_count + len(cells_to_refine) 
    new_r_edges = np.zeros(new_cell_count + 1)
    parent_cell_map = np.zeros(new_cell_count, dtype=int)
    edge_idx, cell_idx = 0, 0

    for i in range(len(domain.deltaL)):
        
        new_r_edges[edge_idx] = domain.r[i]
        edge_idx += 1

        # check if this cell is resonant
        if cell_mask[i]:
            
            midpoint = domain.r[i] + 0.5 * domain.deltaL[i]
            new_r_edges[edge_idx] = midpoint
            edge_idx += 1

            parent_cell_map[cell_idx] = i
            parent_cell_map[cell_idx + 1] = i
            cell_idx += 2
        else:
            # keep non-resonant cell as is
            parent_cell_map[cell_idx] = i
            cell_idx += 1

    # add the end edge of the last cell
    new_r_edges[edge_idx] = domain.r[-1] + domain.deltaL[-1]
    new_r_edges = new_r_edges[:edge_idx + 1]    

    # finalise the new domain arrays
    new_r = new_r_edges[:-1]
    new_deltaL = np.diff(new_r_edges)
    new_rcen = new_r + 0.5 * new_deltaL
    new_cell_count = len(new_r)

    # 3. APPROPRIATELY ASSIGN NEW CELLS THEIR PARAMETERS

    # assign the new subcells their physical quantities
    new_ne = np.zeros(new_cell_count)
    new_Bx = np.zeros(new_cell_count)
    new_By = np.zeros(new_cell_count)
    new_Bz = np.zeros(new_cell_count)

    # calculate B-field unit vectors for parent cells
    parent_Btot = np.sqrt(domain.Bx[parent_cell_map]**2 + domain.By[parent_cell_map]**2 + domain.Bz[parent_cell_map]**2)
    Bx_hat = domain.Bx[parent_cell_map] / parent_Btot
    By_hat = domain.By[parent_cell_map] / parent_Btot
    Bz_hat = domain.Bz[parent_cell_map] / parent_Btot

    for i in range(new_cell_count):
    
        # re-evaluate quantities from cluster profile at new radius
        new_ne[i], B_tot_new = cluster.profile(new_rcen[i])
    
        # construct new B-field
        new_Bx[i] = B_tot_new * Bx_hat[i]
        new_By[i] = B_tot_new * By_hat[i]
        new_Bz[i] = B_tot_new * Bz_hat[i]

    # calculate derived quantities
    new_B = np.sqrt(new_Bx**2 + new_By**2)  
    new_phi = np.arctan2(new_Bx, new_By)   
    new_omega_p = models.omega_p(new_ne)    

    # update the domain with refined arrays
    domain.r = new_r
    domain.rcen = new_rcen
    domain.deltaL = new_deltaL
    domain.ne = new_ne
    domain.Bx = new_Bx
    domain.By = new_By
    domain.Bz = new_Bz
    domain.B = new_B
    domain.phi = new_phi
    domain.omega_p = new_omega_p

    return domain

def iterate_refine_domain(domain,cluster,energies,masses,couplings,scale_factor,n_iterations,res_class):
    """ iterates domain refinement procedure. User can choose number of iterations. """

    for _ in range(n_iterations):

        initial_cell_count = len(domain.deltaL)
        
        domain = refine_domain(domain,cluster,energies,masses,couplings,scale_factor,res_class)

        final_cell_count = len(domain.deltaL)

        if initial_cell_count == final_cell_count:
            break

    return domain

def propagate_with_refined_domain_Nalp(self,domain,cluster,energies,masses,couplings,scale_factor,pol="both",n_iterations=10,res_class='all'):
    """ propagation function which implements the domain refinement algorithm prior to calculation."""

    refined_domain = iterate_refine_domain(domain,cluster,energies,masses,couplings,scale_factor,n_iterations,res_class)

    P, P_radial = propagate_Nalp(self,refined_domain,energies,pol=pol)

    return P, P_radial