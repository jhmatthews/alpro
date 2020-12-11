import alpro.models as models 
#import powerlaw
import numpy as np 
import alpro
import alpro.util as util 


class Survival:
	'''
	High-level class which interfaces with actual ALP
	calculation as well as cluster models to compute 
	survival probability curves
	'''
	def __init__(self, ModelType, implementation="c"):
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

		if implementation == "c":
			self.get_P = alpro.core.get_P 
		elif implementation == "python":
			self.get_P = alpro.pure.get_P
		elif implementation == "numba":
			self.get_P = alpro.get_P

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

	def get_curve(self, energies, random_seed, L, r0=10.0, radial_profile = False, rm_reject = np.inf, propagation="C"):
		if self.model == "libanov":
			self.domain.create_libanov_field()
		else: 
			self.domain.create_box_array(L, random_seed, self.coherence_func, r0=r0) 

		if propagation == "pure":
			propagation_func = self.propagate_pure
		else:
			propagation_func = self.propagate
		# only compute curve if RM Acceptable 
		if self.domain.rm < rm_reject:
			P, P_radial = propagation_func(self.domain, energies)
		else: 
			P, P_radial = None, None 

		if radial_profile:
			return (P, P_radial)
		else:
			return (P)

	def set_params(self, g_a, mass):
		self.g_a = g_a 
		self.mass = mass

	def propagate_pure(self, domain, energies):
		'''
		Propagate an unpolarised beam through a domain and 
		calculate conversion into axion-like particles. Use Pure
		python routines.

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

			P1, Anew1 = alpro.pure.get_P(energies, Ainit1, phi, B, L, self.g_a, self.mass, ne)
			P2, Anew2 = alpro.pure.get_P(energies, Ainit2, phi, B, L, self.g_a, self.mass, ne)

			Ainit1 = Anew1
			Ainit2 = Anew2

			P_radial[i,:] = 0.5 * (P1 + P2)
			#for j in range(len(P1)):
			#	if (P_radial[i,j]<0) or (P_radial[i,j]>1)
			#		print (Anew1[j,:], Anew2[j,:], P1[j], P2[j]) 

		P = 0.5 * (P1 + P2)

		return (P, P_radial)

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
			#for j in range(len(P1)):
			#	if (P_radial[i,j]<0) or (P_radial[i,j]>1)
			#		print (Anew1[j,:], Anew2[j,:], P1[j], P2[j]) 

		P = 0.5 * (P1 + P2)

		return (P, P_radial)
