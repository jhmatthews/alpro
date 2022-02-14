import alpro.models as models 
import numpy as np 
import alpro
import alpro.util as util 
import matplotlib.pyplot as plt

class Survival:
	'''
	High-level class which interfaces with actual ALP
	calculation as well as cluster models to compute 
	survival probability curves. 
	'''
	def __init__(self, ModelType = None, implementation="numba", pol_matrix=False, xmin=3.5, xmax=10.0):
		self.model = ModelType
		self.coherence_func = None

		self.available_models = ["1821", "1275a", "1275b", "custom", "file", "uniform", "libanov", None]

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
			raise ValueError("Model type not recognised, must be one of:", self.available_models)

		if implementation == "c":
			self.get_P = alpro.core.get_P 
		elif implementation == "python":
			self.get_P = alpro.pure.get_P
		elif implementation == "numba":
			self.get_P = alpro.get_P

		self.pol_matrix_bool = pol_matrix

	def set_coherence_pl(self, n=1.2, xmin=3.5, xmax=10.0, r0=None):
		pl = util.my_powerlaw(n=n, xmin=xmin, xmax=xmax)
		self.coherence_func = pl.rvs
		self.set_coherence_r0(r0)

	def show_available_models(self):
	 	print (self.available_models)

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
		initialise the ALP model.  
		'''
		if self.coherence_func == None:
			self.coherence_func = lcorr

		if self.model == "libanov" or self.model == "uniform":
			self.domain = models.FieldModel(profile=None)
		else:
			self.domain = models.FieldModel(profile=self.cluster.profile, coherence_r0=self.coherence_r0)

	def setup_regular_model(self, B = None, L = None, ne = None, N = 1, phi = 0.0):
		self.init_model()

		if self.model == "libanov":
			self.domain.create_libanov_field()
		elif self.model == "uniform":
			self.domain.single_box(phi, B, L, ne, N=1)

	def initialise_domain(self, profile = None):
		self.domain = models.FieldModel(profile = profile)


	# def setup_cell_model(self):

	def get_curve(self, energies, random_seed, L, r0=10.0, radial_profile = False, rm_reject = np.inf, cell_centered=True):
		
		if self.model == "libanov":
			self.domain.create_libanov_field()
		else: 
			self.domain.create_box_array(L, random_seed, self.coherence_func, r0=r0, cell_centered=cell_centered) 

		propagation_func = self.propagate

		# only compute curve if RM Acceptable 
		#if self.domain.rm < rm_reject:
		P, P_radial = propagation_func(self.domain, energies)
		#else: 
		#	P, P_radial = None, None 

		if radial_profile:
			return (P, P_radial)
		else:
			return (P)

	def set_params(self, g = 1e-13 * 1e-9, mass = 1e-13, invgev = False):
		'''
		Set the ALP coupling constant and mass

		Parameters:
			g_a 		float 
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

	def propagate_with_pruning(self, domain, energies, pol="both", threshold = 0.1, refine = 10, required_res = 3, return_pruned_domain=False):
		#resonance_prune(self, mass, threshold = 0.1, refine = 50, required_res = 3)
		domain_pruned = self.domain.resonance_prune(self.mass, 
			            threshold = threshold, refine = refine, required_res = required_res)
		P, P_radial = self.propagate(domain_pruned, energies, pol=pol)
		if return_pruned_domain:
			return (P, P_radial, domain_pruned)
		else:
			return (P, P_radial)


	def propagate(self, domain = None, energies = None, pol="both", pol_matrix = None):
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

		if pol_matrix == None:
			pol_matrix = self.pol_matrix_bool

		# set up the initial state vectors depending on the mode
		if pol_matrix:
			init_x, init_y, _ = pure_pol_vectors_like(energies)
			calculate_P = alpro.get_P_matrix
		else:
			if ypol:
				init_y = np.zeros( (len(energies),3), dtype=np.complex128)
				init_y[:,1] = 1.0

			if xpol:	
				init_x = np.zeros( (len(energies),3), dtype=np.complex128)
				init_x[:,0] = 1.0

			calculate_P = self.get_P

		
		self.P_radial = np.zeros( (len(domain.deltaL),len(energies)) )
		
		# loop over the domain 
		for i in range(len(domain.deltaL)):
			L = domain.deltaL[i]
			B = domain.B[i]
			phi = domain.phi[i] * np.ones_like(energies)
			ne = domain.ne[i]

			if ypol:
				P_y, new_y = calculate_P(energies, init_y, phi, B, L, self.g_a, self.mass, ne)
				init_y = new_y

			if xpol:
				P_x, new_x = calculate_P(energies, init_x, phi, B, L, self.g_a, self.mass, ne)
				init_x = new_x

			if xpol and ypol:
				Ptot = 0.5 * (P_y + P_x)
			elif xpol:
				Ptot = P_x
			elif ypol:
				Ptot = P_y

			# also store the radial profile
			self.P_radial[i,:] = Ptot

		self.P = Ptot
		self.energies = energies

		return (self.P, self.P_radial)

	def default_plot(self, plot_kwargs={}, mode="conversion", theory=None, theory_kwargs={}):
		'''
		Mke a basic plot of the survival probability
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
			plt.plot(self.energies, theory(self.energies, **theory_kwargs), **plot_kwargs, ls="--")
		plt.semilogx()
		plt.xlabel("$E$ (eV)", fontsize=18)
		return (fig)

	def analytic_pga(self, energies):
		P = util.P_gamma_analytic(energies, self.domain.B, self.domain.ne, self.mass, self.g_a, self.domain.deltaL)
		return P



def pure_pol_vectors_like(arr):
	size = len(arr)
	x = np.zeros((size,3,3), dtype=np.complex128)
	y = np.zeros((size,3,3), dtype=np.complex128)
	a = np.zeros((size,3,3), dtype=np.complex128)
	x[:,0,0] = 1.0 + 0.0j
	y[:,1,1] = 1.0 + 0.0j
	a[:,2,2] = 1.0 + 0.0j
	return (x,y,a)