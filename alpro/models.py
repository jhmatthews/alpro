import numpy as np 
from scipy.interpolate import interp1d
import os

def random_angle():
	costheta = 2.0 * np.random.random() - 1.0
	phi = 2.0 * np.pi * np.random.random
	theta = np.arccos(costheta)
	return (theta, phi)

class RussellProfile:
	'''
	container for a pressure profile for the Russell+ 2010
	paper on 1821
	'''
	def __init__(self, plasma_beta = 100):

		self.plasma_beta = plasma_beta

		# read the data, which should be in data subfolder 
		folder = os.path.dirname(__file__)
		filename = "{}/data/russell-pressure.dat".format(folder)
		r, P = np.genfromtxt(filename, unpack=True)

		# convert from Pascals_to_dynes
		P *= 10.0

		# get magnetic field from pressure 
		B = np.sqrt(8.0 * np.pi * P / plasma_beta)

		# create interpolation function
		self.interp_func = interp1d(r, B, kind="linear", fill_value="extrapolate")

	def profile (self, r):
		return self.interp_func(r)

class PerseusProfile:
	def __init__(self, model="a", plasma_beta = 100):
		self.n0 = self.density(0.0)
		self.n25 = self.density(25.0)
		self.plasma_beta = plasma_beta

		if model == "a":
			self.get_B = self.B_modA 
		elif model == "b":
			self.get_B = self.B_modB 

	def B_modA(self, r):
		return (2.5e-5 * (self.density(r) / self.n0)**0.7)

	def B_modB(self, r):
		B25 = 7.5e-6
		B = B25 * np.sqrt(self.density(r) / self.n25 * 100.0 / self.plasma_beta)
		return (B)

	def density(self, r):
		term1 = 3.9e-2 / (1.0 + (r/80.0))**1.8
		term2 = 4.05e-3 / (1.0 + (r/280.0))**0.87
		return (term1 + term2)



class BoxModel:
	def __init__(self, profile, plasma_beta=100):
		self.profile = profile
		self.beta = plasma_beta

	# def create_box_array(self, L, B0, random_seed, coherence_func):
	# 	'''
	# 	create an array of random magnetic field boxes by drawing
	# 	random angles and box sizes from coherence_func. 

	# 	Parameters:
	# 		L				float
	# 						size of domain in kiloparsecs 

	# 		random_seed 	int
	# 						random number seed 

	# 		coherence_func	function
	# 						function that computes coherence length at distance r
	# 	'''


	# 	# r = 0
	# 	w







