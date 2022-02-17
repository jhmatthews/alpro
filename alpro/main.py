import alpro.models as models
import numpy as np
import alpro
import alpro.util as util
import matplotlib.pyplot as plt
import alpro.fourier as fourier


class Survival:
    '''
    High-level class which interfaces with actual ALP
    calculation as well as cluster models to compute
    survival probability curves. This is the starting
    point for any calculation.
    '''

    def __init__(self, ModelType=None, implementation="numba",
                 pol_matrix=False, xmin=3.5, xmax=10.0):
        self.model = ModelType
        self.coherence_func = None

        self.available_models = [
            "1821",
            "1275a",
            "1275b",
            "custom",
            "file",
            "uniform",
            "libanov",
            None]

        if self.model == "1821":
            self.cluster = models.ClusterProfile(model="russell")
            pl = util.my_powerlaw(n=1.2, xmin=xmin, xmax=xmax)
            self.coherence_func = pl.rvs
            self.coherence_r0 = None

        elif self.model == "1275a":
            self.cluster = models.ClusterProfile(model="a")
            pl = util.my_powerlaw(n=1.2, xmin=xmin, xmax=xmax)
            self.coherence_func = pl.rvs
            self.coherence_r0 = None

        elif self.model == "1275b":
            self.cluster = models.ClusterProfile(model="b")
            pl = util.my_powerlaw(n=1.2, xmin=xmin, xmax=xmax)
            self.coherence_func = pl.rvs
            self.set_coherence_r0(50.0)

        elif self.model == "custom":
            self.cluster = models.ClusterProfile(model="custom")
            self.coherence_r0 = None

        elif self.model not in self.available_models:
            raise ValueError(
                "Model type not recognised, must be one of:",
                self.available_models)

        # allow implementation to be pure python
        if implementation == "python":
            self.get_P = alpro.pure.get_P
        elif implementation == "numba":
            self.get_P = alpro.get_P

        self.pol_matrix_bool = pol_matrix

    def set_coherence_pl(self, n=1.2, xmin=3.5, xmax=10.0, r0=None):
        '''
        setup a powerlaw distribution of coherence lengths / box sizes
        '''
        pl = util.my_powerlaw(n=n, xmin=xmin, xmax=xmax)
        self.coherence_func = pl.rvs
        self.set_coherence_r0(r0)

    def show_available_models(self):
        print(self.available_models)

    def set_coherence_r0(self, r0):
        '''
        Set the scale length over which the coherence length varies with radius.
        '''
        # this allows for variation of coherence length with radius in model B
        self.coherence_r0 = r0
        self.init_model()  # must re-initialise model after r0 is set, see issue #4

    def set_churazov_density(self):
        '''
        Manually override a density or B field function to the churazov profile
        '''
        self.cluster.density = self.cluster.churazov_density

    def init_model(self, lcorr=None):
        '''
        Initialise the ALP model. Sets up an appropriate domain.
        '''
        if self.coherence_func is None:
            self.coherence_func = lcorr

        if self.model == "libanov" or self.model == "uniform":
            self.domain = models.FieldModel(profile=None)
        else:
            self.domain = models.FieldModel(
                profile=self.cluster.profile,
                coherence_r0=self.coherence_r0)

    def setup_regular_model(self, B=None, L=None, ne=None, N=1, phi=0.0):
        '''
        set up a regular field model
        '''
        self.init_model()

        if self.model == "libanov":
            self.domain.create_libanov_field()
        elif self.model == "uniform":
            self.domain.single_box(phi, B, L, ne, N=1)

    def initialise_domain(self, profile=None):
        self.domain = models.FieldModel(profile=profile)

    # def setup_cell_model(self):

    def get_curve(self, energies, random_seed, L, r0=10.0,
                  radial_profile=False, rm_reject=np.inf, cell_centered=True):

        if self.model == "libanov":
            self.domain.create_libanov_field()
        else:
            self.domain.create_box_array(
                L,
                random_seed,
                self.coherence_func,
                r0=r0,
                cell_centered=cell_centered)

        propagation_func = self.propagate

        # only compute curve if RM Acceptable
        # if self.domain.rm < rm_reject:
        P, P_radial = propagation_func(self.domain, energies)
        # else:
        #	P, P_radial = None, None

        if radial_profile:
            return (P, P_radial)
        else:
            return (P)

    def set_params(self, g=1e-13 * 1e-9, mass=1e-13, invgev=False):
        '''
        Set the ALP coupling constant and mass

        Parameters:
                g 			float
                                        ALP coupling constant in units set by invgev boolean
                mass 		float
                                        ALP mass in eV
                invgev 		bool
                                        False = g is in eV^-1, True = g is in GeV^-1
        '''
        if invgev:
            g = g * 1e-9
        self.g_a = g
        self.mass = mass

    def propagate_with_pruning(self, domain, energies, pol="both", threshold=0.1,
                               refine=10, required_res=3, return_pruned_domain=False):
        #resonance_prune(self, mass, threshold = 0.1, refine = 50, required_res = 3)
        domain_pruned = self.domain.resonance_prune(self.mass,
                                                    threshold=threshold, refine=refine, required_res=required_res)
        P, P_radial = self.propagate(domain_pruned, energies, pol=pol)
        if return_pruned_domain:
            return (P, P_radial, domain_pruned)
        else:
            return (P, P_radial)

    def propagate(self, domain=None, energies=None,
                  pol="both", pol_matrix=None):
        '''
        Propagate an unpolarised beam through a domain and
        calculate conversion into axion-like particles

        Parameters:
                domain 		object
                                        an alpro.models.FieldModel instance containing
                                        magnetic field components, individual cell sizes
                                        and densities. default is
                                        the domain stored in the class instance.

                energies	array-like
                                        array of photon energies in electron volts. default is
                                        the array stored in the class instance.

                pol 		string (optional)
                                        which polarization state to consider.
                                        must be 'x', 'y' or 'both'

                pol_matrix 	bool or Nonetype
                                        whether to use the polarization matrix pr
                                        just average the two polarization states.
                                        if ==None, then use the value stored in self.pol_matrix_bool
                                        Both produce the same results but using the polarization
                                        matrix could be used to future to predict polarization
                                        signatures.

        Returns:
                P 			array-like
                                        conversion probability as a function of energy

                Pradial 	array-like (2D)
                                        conversion probability as a function of energy and distance
                                        (evaluated at cell centers)
        '''
        if domain is None:
            domain = self.domain
        if energies is None:
            energies = self.energies

        # decide which polarization states to compute based on pol kwarg
        if pol == "both":
            ypol = True
            xpol = True
        elif pol == "x":
            ypol = False
            xpol = True
        elif pol == "y":
            ypol = True
            xpol = False
        else:
            raise ValueError("pol keyword must be 'x', 'y' or 'both'")

        if pol_matrix is None:
            pol_matrix = self.pol_matrix_bool

        # set up the initial state vectors depending on the mode
        if pol_matrix:
            init_x, init_y, _ = pure_pol_vectors_like(energies)
            calculate_P = alpro.get_P_matrix
        else:
            if ypol:
                init_y = np.zeros((len(energies), 3), dtype=np.complex128)
                init_y[:, 1] = 1.0

            if xpol:
                init_x = np.zeros((len(energies), 3), dtype=np.complex128)
                init_x[:, 0] = 1.0

            calculate_P = self.get_P

        self.P_radial = np.zeros((len(domain.deltaL), len(energies)))

        # loop over the domain
        for i in range(len(domain.deltaL)):
            L = domain.deltaL[i]
            B = domain.B[i]
            phi = domain.phi[i] * np.ones_like(energies)
            ne = domain.ne[i]

            if ypol:
                P_y, new_y = calculate_P(
                    energies, init_y, phi, B, L, self.g_a, self.mass, ne)
                init_y = new_y

            if xpol:
                P_x, new_x = calculate_P(
                    energies, init_x, phi, B, L, self.g_a, self.mass, ne)
                init_x = new_x

            if xpol and ypol:
                Ptot = 0.5 * (P_y + P_x)
            elif xpol:
                Ptot = P_x
            elif ypol:
                Ptot = P_y

            # also store the radial profile
            self.P_radial[i, :] = Ptot

        self.P = Ptot
        self.energies = energies

        return (self.P, self.P_radial)

    def default_plot(self, plot_kwargs={}, mode="conversion",
                     theory=None, theory_kwargs={}):
        '''
        Make a basic plot of the survival probability against energy
        using most recent calculation stored in Survival class.

        Parameters:
                plot_kwargs 	dict
                                                dictionary of kwargs to pass to plt.plot

                mode 			str
                                                where to plot "conversion" or "survival" probability

                theory 			function
                                                function to overlay on plot, e.g. s.analytic_pga

                theory_kwargs 	dict
                                                dictionary of kwargs to pass to plt.plot for theory overlay

        Returns:
                fig 			matplotlib.pyplot.figure
                                                figure instance
        '''
        fig = plt.figure()
        if mode == "conversion":
            plt.ylabel(r"$P_{\gamma a}$", fontsize=18)
            P = self.P
        elif mode == "survival":
            plt.ylabel(r"$P_{\gamma\gamma}$", fontsize=18)
            P = 1.0 - self.P

        plt.plot(self.energies, P, **plot_kwargs)
        if theory is not None:
            plt.plot(
                self.energies,
                theory(
                    self.energies,
                    **theory_kwargs),
                **plot_kwargs,
                ls="--")
        plt.semilogx()
        plt.xlabel("$E$ (eV)", fontsize=18)
        return (fig)

    def analytic_pga(self, energies):
        '''
        wrapper to util.P_gamma_analytic that passes the arguments stored
        in the Survival class
        '''
        P = util.P_gamma_analytic(
            energies,
            self.domain.B,
            self.domain.ne,
            self.mass,
            self.g_a,
            self.domain.deltaL)
        return P

    def propagate_fourier(self, domain=None, pol="both",
                          mode="auto", N=100000, f_res=2000, pad_factor=40.0):
        '''
        Propagate an unpolarised beam through a domain and
        calculate conversion into axion-like particles using the fourier formalism

        Parameters:
                domain 		object
                                        an alpro.models.FieldModel instance containing
                                        magnetic field components, individual cell sizes
                                        and densities. default is
                                        the domain stored in the class instance.

                pol 		string (optional)
                                        which polarization state to consider.
                                        must be 'x', 'y' or 'both'

                mode 		string (optional)
                                        How to decide which Fourier formalism to use. The default, "auto" is recommended,
                                        in which case the decision is made by comparing the plasma frequency and the mass.
                                        However, you can force a calculation in the massive or massless regime.

                N 			int
                                        How to decide which Fourier formalism to use. The default, "auto" is recommended,
                                        in which case the decision is made by comparing the plasma frequency and the mass.
                                        However, you can force a calculation in the massive or massless regime.
                                        Recommended value around 1e4 to 1e5.

                f_res 		int
                                        factor by which to resample the resolution of the grid compared to the original resolution.
                                        The domain is zero padded out to a factor (N/f_res) times the original size.
                                        Recommended value is around N/40.

        Returns:
                E 			array-like
                                        Energy in eV

                P 			array-like
                                        conversion probability as a function of energy
        '''
        if domain is None:
            domain = self.domain

        domain.omega_p = alpro.models.omega_p(domain.ne)

        if self.mass > np.max(domain.omega_p):
            regime = "massive"
        elif self.mass < np.min(domain.omega_p):
            regime = "massless"
        else:
            regime = "resonant"

        options = ["auto", "massive", "massless"]
        assert (mode in ["auto", "massive", "massless"]
                ), ("mode must be in", options)

        if mode == "auto":
            if regime == "resonant":
                print("Warning! resonant crossing point found.")
                print(
                    "Adopting massive ALP regime for calculation with inaccurate results")
                print(
                    "m_a: {:.2e} eV, plasma frequency range: {.2e} to {.2e} eV".format(
                        self.mass, np.min(
                            domain.omega_p), np.max(
                            domain.omega_p)))
                regime_to_use = "massive"
            else:
                regime_to_use = regime

        elif mode != regime:
            print(
                "Warning! forced {} ALP calculation in real {} ALP regime!".format(
                    mode, regime))
            print(
                "m_a: {:.2e} eV, plasma frequency range: {:.2e} to {:.2e} eV".format(
                    self.mass, np.min(
                        domain.omega_p), np.max(
                        domain.omega_p)))
            regime_to_use = mode

        else:
            regime_to_use = mode

        if regime_to_use == "massless":
            E, P = fourier.pga_massless(self, f_res=f_res, Nsamples=N, pol=pol)
        elif regime_to_use == "massive":
            E, P = fourier.pga_massive(self, f_res=f_res, Nsamples=N, pol=pol)

        self.P = P
        self.energies = E

        return (E, P)


def pure_pol_vectors_like(arr):
    size = len(arr)
    x = np.zeros((size, 3, 3), dtype=np.complex128)
    y = np.zeros((size, 3, 3), dtype=np.complex128)
    a = np.zeros((size, 3, 3), dtype=np.complex128)
    x[:, 0, 0] = 1.0 + 0.0j
    y[:, 1, 1] = 1.0 + 0.0j
    a[:, 2, 2] = 1.0 + 0.0j
    return (x, y, a)
