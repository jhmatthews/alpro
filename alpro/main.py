import alpro.models as models 
#import powerlaw
import numpy as np 
import alpro
import matplotlib
import alpro.util as util 

class Survival:
	'''
	High-level class which interfaces with actual ALP
	calculation as well as cluster models to compute 
	survival probability curves
	'''
	def __init__(self, ModelType):
		self.model = ModelType
		self.coherence_func = None

		if self.model == "1821":
			self.cluster = models.ClusterProfile(model="russell")
			pl = util.my_powerlaw(n=1.2, xmin=3.5, xmax=10.0)
			self.coherence_func = pl.rvs

		elif self.model == "1275a":
			self.cluster = models.ClusterProfile(model="a")
			pl = util.my_powerlaw(n=1.2, xmin=3.5, xmax=10.0)
			self.coherence_func = pl.rvs

		elif self.model == "1275b":
			self.cluster = models.ClusterProfile(model="b")
			pl = util.my_powerlaw(n=1.2, xmin=3.5, xmax=10.0)
			self.coherence_func = pl.rvs

	def init_model(self, lcorr=None):
		if self.coherence_func == None:
			self.coherence_func = lcorr

		# this allows for variation of coherence length with radius in model B
		if self.model == "1275b":
			coh_r0 = 50.0
		else:
			coh_r0 = None

		if self.model == "libanov":
			self.domain = models.FieldModel(profile=None)
		else:
			self.domain = models.FieldModel(profile=self.cluster.profile, coherence_r0 = coh_r0)

	def get_curve(self, energies, random_seed, L, r0=10.0, radial_profile = False, rm_reject = np.inf):
		if self.model == "libanov":
			self.domain.create_libanov_field()
		else: 

			self.domain.create_box_array(L, random_seed, self.coherence_func, r0=r0) 

		# only compute curve if RM Acceptable 
		if self.domain.rm < rm_reject:
			P, P_radial = self.propagate(self.domain, energies)
		else: 
			P, P_radial = None, None 

		if radial_profile:
			return (P, P_radial)
		else:
			return (P)

	def set_params(self, g_a, mass):
		self.g_a = g_a 
		self.mass = mass

	def propagate(self, domain, energies):
		'''
		Propagate an unpolarised beam through a domain and 
		calculate conversion into axion-like particles 

		Parameters:
			domain 		object
						an alpro.models.FieldModel instance containing
						magnetic field components, individual cell sizes
						and densities.
			energies	array-like
						array of photon energies in electron volts 
		'''

		Ainit1 = np.zeros( (len(energies),6))
		Ainit1[:,2] = 1.0

		Ainit2 = np.zeros( (len(energies),6))
		Ainit2[:,0] = 1.0

		P_radial = np.zeros( (len(domain.r),len(energies)) )
		
		for i in range(len(domain.r)):
			L = domain.deltaL[i]
			B = domain.B[i]
			phi = domain.phi[i] * np.ones_like(energies)
			ne = domain.ne[i]

			P1, Anew1 = alpro.get_P(energies, Ainit1, phi, B, L, self.g_a, self.mass, ne)
			P2, Anew2 = alpro.get_P(energies, Ainit2, phi, B, L, self.g_a, self.mass, ne)

			Ainit1 = Anew1
			Ainit2 = Anew2

			P_radial[i,:] = 0.5 * (P1 + P2)

		P = 0.5 * (P1 + P2)

		return (P, P_radial)
