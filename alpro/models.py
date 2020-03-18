import numpy as np 
from scipy.interpolate import interp1d
import types, os

class units:
    def __init__(self):
        self.kpc = 3.086e21
        self.pc = 3.086e18
        self.c = 2.997925e10
        self.yr = 3.1556925e7
        self.myr = 3.1556925e13
        self.kyr = 3.1556925e10
        self.radian = 57.29577951308232
        self.msol = 1.989e33
        self.mprot = 1.672661e-24
        self.ev = 1.602192e-12
        self.kb = 1.38062e-16
        self.h = 6.6262e-27
        self.g = 6.670e-8

# class to use for units
unit = units()

def random_angle():
	'''
	compute a random isotropic angle 
	'''
	costheta = 2.0 * np.random.random() - 1.0
	phi = 2.0 * np.pi * np.random.random()
	theta = np.arccos(costheta)
	return (theta, phi)

def get_libanov_B(r):
	x = r/(93.0)
	Bx = (0.00312443*(x**18)) - (0.0319991*(x**16)) + (0.260311*(x**14)) - (1.63197*(x**12)) + (7.58002*(x**10)) - (24.721*(x**8)) + (52.3929*(x**6)) - (63.8794*(x**4)) + (35.8973*(x**2)) - 5.86899	
	By = (0.0102459*(x**17))-(0.0937683*(x**15)) + (0.671841*(x**13)) - (3.6406*(x**11)) + (14.2479*(x**9)) - (37.7455*(x**7)) + (61.3611*(x**5)) - (51.7231*(x**3)) + (16.9128*x)
	return 1e-6*Bx, 1e-6*By

# class powerlaw:
# 	def __init__

# 	# def generate_random(r):

class ClusterProfile:
	'''
	container for a magentic field and density profile for a cluster
	'''
	def __init__(self, model="a", plasma_beta = 100, t_keV = 1.0):
		self.plasma_beta = plasma_beta

		if model == "a":
			# Model B from Reynolds+ 2020
			self.get_B = self.B_modA 
			self.density = self.churasov_density
			self.n0 = self.density(0.0)
			self.n25 = self.density(25.0)

		elif model == "b":
			# Model B from Reynolds+ 2020
			self.get_B = self.B_modB 
			self.density = self.churasov_density
			self.n0 = self.density(0.0)
			self.n25 = self.density(25.0)

		elif model == "russell":
			# Russell+ 2010 paper on 1821
			# read the data, which should be in data subfolder 
			folder = os.path.dirname(__file__)
			filename = "{}/data/russell-pressure.dat".format(folder)
			r, P = np.genfromtxt(filename, unpack=True)

			filename = "{}/data/russell-density.dat".format(folder)
			r2, n = np.genfromtxt(filename, unpack=True)

			# convert from Pascals_to_dynes
			P *= 10.0

			# get magnetic field from pressure 
			B = np.sqrt(8.0 * np.pi * P / plasma_beta)

			#kT = t_keV * unit.ev * 1000.0
			#n = P / kT

			# create interpolation functions
			self.get_B = interp1d(r, B, kind="quadratic", fill_value="extrapolate")

			# create interpolation functions
			self.density = interp1d(r2, n, kind="quadratic", fill_value="extrapolate")

		else:
			raise ValueError("ClusterProfile did not understand model type {}".format(model))

	def B_modA(self, r):
		'''
		Model A from Reynolds et al. 2020

		Parameters:
			r 	float
				distance from cluster centre in kpc
		'''
		return (2.5e-5 * (self.density(r) / self.n0)**0.7)

	def B_modB(self, r):
		'''
		Model B from Reynolds et al. 2020

		Parameters:
			r 	float
				distance from cluster centre in kpc
		'''
		B25 = 7.5e-6
		B = B25 * np.sqrt(self.density(r) / self.n25 * 100.0 / self.plasma_beta)
		return (B)

	def churasov_density(self, r):
		'''
		Density function from Churasov et al 2003
		Given as equation (2) of Reynolds et al. 2020

		Parameters:
			r 	float
				distance from cluster centre in kpc
		'''
		term1 = 3.9e-2 / (1.0 + (r/80.0))**1.8
		term2 = 4.05e-3 / (1.0 + (r/280.0))**0.87
		return (term1 + term2)

	def profile(self, r):
		return (self.density(r), self.get_B(r))



class FieldModel:
	def __init__(self, profile, plasma_beta=100):
		self.profile = profile
		self.beta = plasma_beta

	def create_libanov_field(self, deltaL=1.0, Lmax=93.0):
		self.r = np.arange(0, Lmax-deltaL, deltaL)
		self.deltaL = np.ones_like(self.r) * deltaL
		self.rcen = self.r + (0.5 * self.deltaL)
		self.Bx, self.By = get_libanov_B(self.rcen)
		self.B = np.sqrt(self.Bx**2 + self.By**2) 
		self.phi = np.arctan(self.Bx/self.By) 
		self.ne = 1e-20 * np.ones_like(self.r)	# vanishing

	def create_box_array(self, L, random_seed, coherence, r0=10.0):
		'''
		create an array of random magnetic field boxes by drawing
		random angles and box sizes from coherence_func. 

		Parameters:
			L				float
							size of domain in kiloparsecs 

			random_seed 	int
							random number seed 

			coherence_func	function or float
							function that computes coherence length at distance r,
							or a single-value floating point number if the coherence
							length is constant.

			r0 				float
							inner radius of the calculation (used to excide an inner region)
		'''

		if isinstance(coherence, float) == False and callable(coherence) == False:
			raise TypeError("kwarg coherence must be callable or a float.")

		# set random number seed 
		np.random.seed(random_seed)

		# initialise arrays and counters 
		r = r0
		rcen = r0
		Bx_array, By_array = [], []
		rcen_array, r_array = [], []
		B_array, ne_array = [], []
		deltaL_array = []

		# IMPROVE?
		# wonder if there's a better way to do this?
		# vectorise it, you mug.
		# print (r, L)
		while r < L:
			# get a coherence length which will be the size of the box 
			# this can be a function or a float 
			if callable(coherence):
				lc = coherence()
			else:
				lc = coherence

			if rcen == r0:
				rcen += lc / 2.0
			else:
				rcen += lc

			# draw a random isotropic angle 
			theta, phi = random_angle()

			# get density and magnetic field strength at centre of box
			density, B = self.profile(rcen)

			# get the x and y components and increment r
			B_array.append(B)
			ne_array.append(density)
			Bx_array.append(B * np.sin(theta) * np.cos(phi))
			By_array.append(B * np.sin(theta) * np.sin(phi))
			rcen_array.append(rcen)
			r_array.append(r)
			deltaL_array.append(lc)
			r += lc

		# copy to class 
		self.deltaL = np.array(deltaL_array)
		self.ne = np.array(ne_array)
		self.B = np.array(B_array)
		self.Bx = np.array(Bx_array)
		self.By = np.array(By_array)
		self.phi = np.arctan(self.Bx/self.By) 
		self.r = np.array(r_array)
		self.rcen = np.array(rcen_array)







