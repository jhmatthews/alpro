from alpro.models import unit
import scipy.fftpack as fftpack
import numpy as np

# wrappers to the scipy DCT and DST to make sure we always divide by 2
# this factor 2 is in the definition for the DCT/DST in scipy docs
# but conventions vary and we don't use it in our formalism


def dct(*args, **kwargs):
    return (fftpack.dct(*args, **kwargs) / 2.0)


def dst(*args, **kwargs):
    return (fftpack.dst(*args, **kwargs) / 2.0)


def autocorr(x):
    '''
    get the (unitless) autocorrelation function of the array "x"
    '''
    result = np.correlate(x, x, mode='full')
    return result[result.size // 2:]


def get_autocorr(dx, y):
    '''
    get the autocorrelation function in the right units for a spacing dx
    '''
    c_f = dx * autocorr(y)
    return (c_f)


def get_transforms(x, g, B, omegap=None, mass=None,
                   mode="massive", transform_type=1, coord_units="kpc"):
    '''
    calculate the cosine and sine transforms of B

    Parameters:
            x 				array-like
                                            the coordinate for the calculation - either physical space
                                            or phi

            g 				float
                                            coupling constant in eV^-1

            B 				array-like
                                            magnetic field in Gauss

            mass 			float or Nonetype
                                            ALP mass in eV

            mode 			str
                                            massive or massless

    Returns:
            kfreq 			array-like
                                            Fourier samples of eta or lambda
                                            (conjugate variable depends on regime)
            dsx				array-like
                                            the discrete sine transform
            dcx				array-like
                                            the discrete cosine transform

    '''
    N = len(x) + 1

    if coord_units == "kpc":
        xnat = x * unit.kpc * unit.unit_length_natural
    elif coord_units == "natural":
        xnat = x
    else:
        raise ValueError(
            "Did not understand unit {}. Must be 'kpc' or 'natural'".format(coord_units))

    if mode == "massive":
        if mass is None:
            raise ValueError("mass kwarg must be set in massive case!")

        delta_x_array = g * B * unit.unit_gauss_natural / 2.0

    elif mode == "massless":
        if omegap is None:
            raise ValueError(
                "get_transforms: omegap kwarg must be set in massless case!")

        delta_x_array = g * B * unit.unit_gauss_natural / omegap / omegap

    # phi or r must be uniformly spaced
    dx_all = np.diff(xnat)
    dx = dx_all[0]
    dcx = dx * dct(delta_x_array, type=transform_type)
    dsx = dx * dst(delta_x_array, type=transform_type)
    kfreq = np.pi * np.arange(0, N - 1, 1) / (dx * N)

    return (kfreq, dsx, dcx)


def get_kfreq(dx, N):
    kfreq = np.pi * np.arange(0, N - 1, 1) / (dx * N)
    return (kfreq)


def get_energy_samples(x, mass=None, omegap=None,
                       mode="massive", coord_units="kpc"):

    if coord_units == "kpc":
        xnat = x * unit.kpc * unit.unit_length_natural
    elif coord_units == "natural":
        xnat = x
    else:
        raise ValueError(
            "Did not understand unit {}. Must be 'kpc' or 'natural'".format(coord_units))

    # get number of points + 1
    N = len(x) + 1

    # phi or r must be uniformly spaced
    dx_all = np.diff(xnat)
    dx = dx_all[0]

    # check evenly spaced array within tolerance
    assert (np.all(np.fabs(dx_all - dx) / dx < 1e-10))

    kfreq = get_kfreq(dx, N)

    if mode == "massive":
        if mass is None:
            raise ValueError(
                "get_pgg_transforms: mass kwarg must be set in massive case!")
        with np.errstate(divide='ignore'):
            E = mass * mass / 2.0 / kfreq

    elif mode == "massless":
        if omegap is None:
            raise ValueError(
                "get_pgg_transforms: omegap kwarg must be set in massless case!")
        with np.errstate(divide='ignore'):
            E = 1.0 / kfreq

    else:
        raise ValueError(
            "Did not understand mode {}. Must be 'massive' or 'massless'".format(mode))

    return (E)


def get_pga_transforms(x, g, B, omegap=None, mass=None, mode="massive",
                       transform_type=1, coord_units="kpc", return_transforms=True):
    '''
    Calculate a conversion probability using the autocorrelation function

    Parameters:
            x 				array-like
                                            the coordinate for the calculation - either physical space
                                            or phi

            g 				float
                                            coupling constant in eV^-1

            B 				array-like
                                            magnetic field in Gauss

            omegap			array-like
                                            plasma frequency in eV

            mode 			string
                                            whether to use massive ALPs or massless ALPs formalism

            return_autocorr	bool
                                            whether to return the autocorrelation function

    Returns:
            E 				array-like
                                            energies in eV

            P 				array-like
                                            conversion probability at these energies

            kfreq			array-like
                                            conjugate variable - either lambda or eta depending on mode

            c_f 			array-like
                                            autocorrelation function
    '''

    # get number of points + 1
    N = len(x) + 1

    if coord_units == "kpc":
        xnat = x * unit.kpc * unit.unit_length_natural
    elif coord_units == "natural":
        xnat = x
    else:
        raise ValueError(
            "Did not understand unit {}. Must be 'kpc' or 'natural'".format(coord_units))

    # phi or r must be uniformly spaced
    dx_all = np.diff(xnat)
    dx = dx_all[0]

    # check evenly spaced array within tolerance
    assert (np.all(np.fabs(dx_all - dx) / dx < 1e-10)), "Array of Fourier samples is not evenly spaced!"

    # call function that calculates sine and cosine transforms
    # exact calculation depends on the mode
    kfreq, dsx, dcx = get_transforms(xnat, g, B, omegap=omegap, mass=mass,
                                     mode=mode, transform_type=transform_type, coord_units="natural")
    P = (dsx**2) + (dcx**2)

    if mode == "massive":
        if mass is None:
            raise ValueError(
                "get_pgg_transforms: mass kwarg must be set in massive case!")

        with np.errstate(divide='ignore'):
            E = mass * mass / 2.0 / kfreq

    elif mode == "massless":
        if omegap is None:
            raise ValueError(
                "get_pgg_transforms: omegap kwarg must be set in massless case!")
        with np.errstate(divide='ignore'):
            E = 1.0 / kfreq

    else:
        raise ValueError(
            "Did not understand mode {}. Must be 'massive' or 'massless'".format(mode))

    if return_transforms:
        return (E, P, kfreq, dcx, dsx)
    else:
        return (E, P, kfreq)


def get_pga_autocorr(x, g, B, omegap=None, mass=None, mode="massive",
                     return_autocorr=False, transform_type=1, coord_units="kpc"):
    '''
    Calculate a conversion probability using the cosine and sine transforms

    Parameters:
            x 				array-like
                                            the coordinate for the calculation - either physical space
                                            or phi

            g 				float
                                            coupling constant in eV^-1

            B 				array-like
                                            magnetic field in Gauss

            omegap			array-like
                                            plasma frequency in eV

            mode 			string
                                            whether to use massive ALPs or massless ALPs formalism

            return_autocorr	bool
                                            whether to return the autocorrelation function

    Returns:
            E 				array-like
                                            energies in eV

            P 				array-like
                                            conversion probability at these energies

            kfreq			array-like
                                            conjugate variable - either lambda or eta depending on mode

            c_f 			array-like
                                            autocorrelation function
    '''

    # get number of points + 1
    N = len(x) + 1

    if coord_units == "kpc":
        xnat = x * unit.kpc * unit.unit_length_natural
    elif coord_units == "natural":
        xnat = x
    else:
        raise ValueError(
            "Did not understand unit {}. Must be 'kpc' or 'natural'".format(coord_units))

    # phi or r must be uniformly spaced
    dx_all = np.diff(xnat)
    dx = dx_all[0]

    # get the conjugate variable points
    kfreq = np.pi * np.arange(0, N - 1, 1) / (dx * N)

    # check evenly spaced array within tolerance
    assert (np.all(np.fabs(dx_all - dx) / dx < 1e-10)), "Array of Fourier samples is not evenly spaced!"

    if mode == "massive":
        if mass is None:
            raise ValueError("mass kwarg must be set in massive case!")

        c_f = get_autocorr(dx, B * unit.unit_gauss_natural)

        # get the transform of the autocorrelation function
        # we divide by two here because the definition used by scipy
        # has an additional factor in that we don't want.
        cos_transform = dx * dct(c_f, type=transform_type)

        with np.errstate(divide='ignore'):
            E = mass * mass / 2.0 / kfreq
        P = g * g / 2.0 * cos_transform
        corr = c_f

    elif mode == "massless":
        if omegap is None:
            raise ValueError("omegap kwarg must be set in massless case!")

        # G is already in natural units.
        # need to multiply by dphi to account for fact that this
        # is an integral over phi
        Gphi = g * B * unit.unit_gauss_natural / omegap / omegap
        c_G = get_autocorr(dx, Gphi)

        # get the transform of the autocorrelation function
        # we divide by two here because the definition used by scipy
        # has an additional factor in that we don't want.
        cos_transform = dx * dct(c_G, type=transform_type)
        probability = 2.0 * cos_transform

        with np.errstate(divide='ignore'):
            E = 1.0 / kfreq
        P = probability
        corr = c_G

    else:
        raise ValueError(
            "Did not understand mode {}. Must be 'massive' or 'massless'".format(mode))

    if return_autocorr:
        return (E, P, kfreq, corr)
    else:
        return (E, P, kfreq)
