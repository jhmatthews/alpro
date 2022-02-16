from .fourier_core import *
import alpro
from scipy.integrate import simps
from scipy.interpolate import interp1d 

def resample_fourier(new_redge):
        '''
        Resample a box array on to a new 1D grid using 1d interpolation.
        Must be called after the Bx, By, r arrays are already populated.
        '''
        interp1d_kwargs={"kind": "nearest", "fill_value": 0}

        interp_array_r = np.concatenate( (self.r[0:1], self.rcen, self.r[-1:] + self.deltaL[-1:]))
        interp_Bx = np.concatenate( (self.Bx[0:1], self.Bx, self.Bx[-1:]))
        interp_By = np.concatenate( (self.By[0:1], self.By, self.By[-1:]))
        interp_x = interp1d(interp_array_r, interp_Bx, **interp1d_kwargs)
        interp_y = interp1d(interp_array_r, interp_By, **interp1d_kwargs)

        interp_Bz = np.concatenate( (self.Bz[0:1], self.Bz, self.Bz[-1:]))
        interp_z = interp1d(interp_array_r, interp_Bz, **interp1d_kwargs)

        # populate new values 
        self.rcen = 0.5 * (new_redge[1:] + new_redge[:-1])
        self.Bx = interp_x(self.rcen)
        self.By = interp_y(self.rcen)
        self.Bz = interp_z(self.rcen)
        self.r = new_redge[:-1]
        self.deltaL = new_redge[1:] - new_redge[:-1]
        self.B = np.sqrt(self.Bx**2  + self.By **2)
        self.phi = np.arctan2(self.Bx,self.By) 
        if profile:
            self.ne, _ = self.profile(self.rcen)
        else:
            interp_ne = np.concatenate( (self.ne[0:1], self.ne, self.ne[-1:]))
            interp_n = interp1d(interp_array_r, interp_ne, **interp1d_kwargs)
            self.ne = interp_n(self.rcen)


def get_phi_numerical(s):
    '''
    Calculate phi numerically from an alpro.Survival class instance 
    using Simpsons integration. Used in the massless ALP Fourier formalism

    Parameters
    -------------
    s       alpro.Survival 
            survival class instance to use for calculation 

    Returns
    -------------
    phi     array-like
            phi coordinate
    '''
    omegap_sq = s.domain.omega_p ** 2
    r_natural = s.domain.rcen * unit.unit_length_natural * unit.kpc
    phi = 0.5 * np.array([simps(omegap_sq[:i+1], x = r_natural[:i+1]) for i in range(0,len(r_natural)-1)])
    return phi


def pga_massive(s, Nsamples=200000, f_res = 4, pol = "both", return_autocorr=False):

    zmax = np.max(s.domain.r) + s.domain.deltaL[-1]
    pad_factor =  (Nsamples / f_res)

    new_z = np.linspace(0, pad_factor * zmax, Nsamples)
    new_z = np.linspace(0, (Nsamples / f_res) * zmax, Nsamples)
    zcen = 0.5 * (new_z[1:] + new_z[:-1])
    new_z_use = new_z[new_z<=zmax]

    domain_old = s.domain.resample_box(new_z_use, interp1d_kwargs={"kind": "nearest", "fill_value": 0, "bounds_error": False}, profile=False)

    Bxnew = np.concatenate( (s.domain.Bx, np.zeros(len(zcen)-len(s.domain.Bx))) )
    Bynew = np.concatenate( (s.domain.By, np.zeros(len(zcen)-len(s.domain.By))) )
    omegap = alpro.models.omega_p(s.domain.ne)
    omegap = np.concatenate( (s.domain.omega_p, np.zeros(len(zcen)-len(s.domain.omega_p))) )

    if pol == "x" or pol == "both":
        Ex, Px, kfreq, corr = get_pga_autocorr(zcen, s.g_a, Bxnew, mass=s.mass, mode="massive", return_autocorr=True, transform_type=1, coord_units="kpc")

    if pol == "y" or pol == "both":
        Ey, Py, kfreq, corr = get_pga_autocorr(zcen, s.g_a, Bynew, mass=s.mass, mode="massive", return_autocorr=True, transform_type=1, coord_units="kpc")

    if pol == "x":
        E = Ex 
        P = Px
    elif pol == "y":
        E = Ey 
        P = Py 
    elif pol == "both":
        E = Ex
        P = 0.5 * (Px + Py)

    #r_to_return = s.domain.rcen

    #s.domain = domain_old

    if return_autocorr:
        return (E, P, r_to_return, corr)
    else:
        return (E, P)




def physical_from_numerics(s, Nsamples=100000, pad_factor =10.0, pol="both"):
    
    phi = get_phi_numerical(s)
    start = 0
    end = 0

    PHIMAX = phi[-1]
    new_phi = np.linspace(0, pad_factor * PHIMAX, Nsamples)
    phicen = 0.5 * (new_phi[1:] + new_phi[:-1])
    Bxnew = Bynew = None

    N = len(phi) + end
    if pol == "x" or pol == "both":
        interp_func = interp1d(phi, s.domain.Bx[start:N], kind="nearest", fill_value=0, bounds_error=False)
        Bxnew = interp_func(phicen)
    if pol == "y" or pol == "both":
        interp_func = interp1d(phi, s.domain.By[start:N], kind="nearest", fill_value=0, bounds_error=False)
        Bynew = interp_func(phicen)

    interp_omega = interp1d(phi, s.domain.omega_p[start:N], kind="nearest", fill_value=1e-20, bounds_error=False)
    omegap = interp_omega (phicen)
    length_calibration = unit.kpc * unit.unit_length_natural

    return (phicen / length_calibration, omegap, Bxnew, Bynew)


def pga_massless(s, Nsamples=200000, f_res = 4, pol="both"):

    pad_factor =  (Nsamples / f_res)

    phicen, omegap, Bx, By = physical_from_numerics(s, Nsamples = Nsamples, pad_factor = pad_factor, pol = pol)

    if pol == "x" or pol == "both":
        E, Px, k = get_pga_autocorr(phicen, s.g_a, Bx, omegap=omegap, mode="massless")
    if pol == "y" or pol == "both":
        E, Py, k = get_pga_autocorr(phicen, s.g_a, By, omegap=omegap, mode="massless")

    if pol == "both":
        P = 0.5 * (Px + Py)
    elif pol == "x":
        P = Px
    elif pol == "y":
        P = Py

    return (E, P)