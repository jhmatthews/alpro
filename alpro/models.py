import numpy as np 
from scipy.interpolate import interp1d
import types, os
from scipy.integrate import simps

class units:
	'''
	class containing some units. Should probably use astropy units 
	but I find them a bit annoying.
	'''
	def __init__(self):
		self.kpc = 3.0857e21
		self.pc = 3.0857e18
		self.c = 2.997925e10
		self.yr = 3.1556925e7
		self.myr = 3.1556925e13
		self.kyr = 3.1556925e10
		self.radian = 57.29577951308232
		self.msol = 1.989e33
		self.mprot = 1.672661e-24
		self.melec = 9.10938356e-28
		self.melec_csq = self.melec * self.c * self.c
		self.mprot_csq = self.mprot * self.c * self.c
		self.e = 4.8032045057134676e-10     # fundamental charge 
		self.kev = 1.602192e-9  # kilo electron volts in CGS
		self.ev = 1.602192e-12  # electron volts in CGS
		self.kb = 1.38062e-16   # boltzmann 
		self.h = 6.6262e-27     # plank 
		self.hbar = self.h / np.pi / np.pi      
		self.g = 6.670e-8       # gravitational 
		self.hbar_c = self.hbar * self.c
		self.alpha = self.e * self.e / self.hbar_c
		self.thomson = 0.66524e-24

# class to use for units
unit = units()

def random_angle(size=None):
	'''
	compute a random isotropic angle 
	'''
	costheta = (2.0 * np.random.random(size=size)) - 1.0
	phi = 2.0 * np.pi * np.random.random(size=size)
	theta = np.arccos(costheta)
	return (theta, phi)

def fterm(r, C, alpha):
	term1 = -alpha * np.cos(alpha * r)
	term2 = np.sin(alpha * r) / r
	F0 = C * alpha * alpha * (alpha * np.cos(alpha) - np.sin(alpha))
	term3 = F0 * r * r / alpha / alpha
	f = C * (term1 + term2) + term3
	return f

def fprime(r, C, alpha):
	term1 = alpha * alpha * np.sin(alpha * r)
	term2 = (alpha * np.cos(alpha * r) / r) - (np.sin(alpha * r) / r / r)
	F0 = C * alpha * alpha * (alpha * np.cos(alpha) - np.sin(alpha))
	term3 = 2.0 * F0 * r / alpha / alpha
	f = C * (term1 + term2) + term3
	return f
    

def libanov_Br(r, alpha=5.76, theta=np.pi/4.0, C=6e-8):
	f =fterm(r, C, alpha)
	Br = 2.0 * np.cos(theta) * f / r / r
	return (-Br)

def get_libanov_B(r, theta=np.pi/4, Rcavity=93.0, alpha=5.76, C=6e-8):
    rnorm = r/Rcavity
    fr = fterm(rnorm, C, alpha)
    Br = 2.0 * np.cos(theta) * fr / rnorm / rnorm
    Btheta = -np.sin(theta) * fprime(rnorm, C, alpha) / rnorm
    Bphi = alpha * np.sin(theta) * fr / rnorm
    return (Br, Btheta, Bphi)

def get_libanov_B_old(r, include_radial=True):
	x = r/(93.0)
	Bx = (0.00312443*(x**18)) - (0.0319991*(x**16)) + (0.260311*(x**14)) - (1.63197*(x**12)) + (7.58002*(x**10)) - (24.721*(x**8)) + (52.3929*(x**6)) - (63.8794*(x**4)) + (35.8973*(x**2)) - 5.86899	
	By = (0.0102459*(x**17))-(0.0937683*(x**15)) + (0.671841*(x**13)) - (3.6406*(x**11)) + (14.2479*(x**9)) - (37.7455*(x**7)) + (61.3611*(x**5)) - (51.7231*(x**3)) + (16.9128*x)
	
	if include_radial:
		Bz = libanov_Br(x)
		return 1e-6*Bx, 1e-6*By, Bz
	else:
		return 1e-6*Bx, 1e-6*By

def churasov_density(r):
	'''
	Density function from Churasov et al 2003
	Given as equation (2) of Reynolds et al. 2020

	Parameters:
		r 	float
			distance from cluster centre in kpc
	'''
	term1 = 3.9e-2 / (( 1.0 + (r/80.0)**2)**1.8)
	term2 = 4.05e-3 / ((1.0 + (r/280.0)**2)**0.87)
	return (term1 + term2)

class ClusterProfile:
	'''
	container for a magentic field and density profile for a cluster
	'''
	def __init__(self, model="a", plasma_beta = 100, B_rms=None, n=None):
		
		self.plasma_beta = plasma_beta

		if model == "a":
			# Model B from Reynolds+ 2020
			self.get_B = self.B_modA 
			self.density = churasov_density
			self.n0 = self.density(0.0)
			self.B0 = 2.5e-5
			self.B_exponent = 0.7

		elif model == "b":
			# Model B from Reynolds+ 2020
			self.get_B = self.B_modB 
			self.density = churasov_density
			self.n25 = self.density(25.0)

		elif model == "flat":
			'''
			allow to just have a uniform field
			'''
			self.B_rms = B_rms
			self.n = n
			self.get_B = self.Bflat
			self.density = churasov_density

		elif model == "murgia":
			self.n0 = 1e-3 
			self.r0 = 400.0 
			self.B_exponent = 0.5
			self.beta = 0.6
			self.B0 = 1e-6
			self.density = self.beta_density
			self.get_B = self.B_modA 


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
			self.get_B = interp1d(r, B, kind="slinear", fill_value="extrapolate")

			# create interpolation functions
			self.density = interp1d(r2, n, kind="slinear", fill_value="extrapolate")
		elif model == "custom":
			print ("Warning: Custom model specified - need to make sure get_B and density methods are populated!")

		else:
			raise ValueError("ClusterProfile did not understand model type {}".format(model))

	def beta_density(self, r):
		exponent = -3.0 * self.beta / 2.0
		n = self.n0 * (1 + (r/self.r0)**2) ** exponent
		return (n)

	def Bflat(self):
		return (self.B_rms)

	def nflat(self):
		return (self.n)

	def B_modA(self, r):
		'''
		Model A from Reynolds et al. 2020

		Parameters:
			r 	float
				distance from cluster centre in kpc
		'''
		return (self.B0 * (self.density(r) / self.n0)**self.B_exponent)

	def B_modB(self, r, B25=7.5e-6):
		'''
		Model B from Reynolds et al. 2020

		Parameters:
			r 	float
				distance from cluster centre in kpc
		'''
		B = B25 * np.sqrt(self.density(r) / self.n25 * 100.0 / self.plasma_beta)
		return (B)

	def profile(self, r):
		return (self.density(r), self.get_B(r))


class ClusterFromFile:
	def __init__(self, fname="Bfield.npy", model_type="cube"):
		# load the array
		self.Bfull = np.load(fname)
		self.N = self.Bfull.shape[0]
		self.mid = self.N//2
		self.density = churasov_density

		if model_type == "cube":
			if any(i != self.N for i in self.Bfull.shape[:-1]):
				raise ValueError("File supplied must be cube shaped but has shape {}".format(self.Bfull.shape)) 
		elif model_type == "1d":
			self.z = self.Bfull[0,:]
			self.B = np.transpose(self.Bfull[1:,:])
			interp_x_temp = interp1d(self.z, self.B[:,0], kind='slinear')
			interp_y_temp = interp1d(self.z, self.B[:,1], kind='slinear')
			interp_z_temp = interp1d(self.z, self.B[:,2], kind='slinear')
			# actual interpolation always done using 2nd order interp 
			kind = "quadratic"
			# kind='quadratic'
			self.interp_x = interp1d(self.z, self.B[:,0], kind=kind)
			self.interp_y = interp1d(self.z, self.B[:,1], kind=kind)
			self.interp_z = interp1d(self.z, self.B[:,2], kind=kind)


	def slice(self, z, L=100.0, axis=0, sign=1, degrade=1, normalise = 1.0):

		if axis == 0:
			self.B = self.Bfull[:,self.mid,self.mid,:]
		elif axis == 1:
			self.B = self.Bfull[self.mid,:,self.mid,:]
		elif axis == 2:
			self.B = self.Bfull[self.mid,self.mid,:,:]

		if sign > 0:
			self.B = self.B[self.mid:,:]
		else:
			self.B = self.B[:self.mid,:]


		self.B *= normalise
		# take a slice along the B field

		from scipy.interpolate import interp1d
		ztrue = z
		self.z = np.linspace(0,L,len(self.B[:,0])/degrade)

		if degrade > 1:
			# these functions will allow us to degrade the resolution using linear spline interp
			interp_x_temp = interp1d(ztrue, self.B[:,0], kind='slinear')
			interp_y_temp = interp1d(ztrue, self.B[:,1], kind='slinear')
			interp_z_temp = interp1d(ztrue, self.B[:,2], kind='slinear')

			self.B = np.zeros((len(self.z),3))

			self.B[:,0] = interp_x_temp (self.z)
			self.B[:,1] = interp_y_temp (self.z)
			self.B[:,2] = interp_z_temp (self.z)

		elif degrade < 1:
			raise ValueError("degrade needs to be >= 1!")

		# actual interpolation always done using 2nd order interp 
		kind = "quadratic"
		# kind='quadratic'
		self.interp_x = interp1d(self.z, self.B[:,0], kind=kind)
		self.interp_y = interp1d(self.z, self.B[:,1], kind=kind)
		self.interp_z = interp1d(self.z, self.B[:,2], kind=kind)


	def get_Bz(self, z):
		Bz = self.interp_z(z)
		return (Bz)

	def get_B(self, z):
		Bx = self.interp_x(z)
		By = self.interp_y(z)
		return (Bx, By)

	def profile(self, r):
		return (self.density(r), self.get_B(r))


class FieldModel:
	def __init__(self, profile, plasma_beta=100, coherence_r0 = None):
		self.profile = profile
		self.beta = plasma_beta
		# coherence_r0 scales the coherence lengths with radius 
		# by a factor of (1 + r/coherence_r0), in kpc
		self.coherence_r0 = coherence_r0
		self.Bz = 1.0

	def create_libanov_field(self, deltaL=1.0, Lmax=93.0, density=None, theta=np.pi/4.0):
		'''
		set arrays according to uniform field model of Libanox et al.
		'''
		self.r = np.arange(0, Lmax, deltaL)
		self.deltaL = np.ones_like(self.r) * deltaL
		self.rcen = self.r + (0.5 * self.deltaL)
		self.Bx, self.By, self.Bz = get_libanov_B(self.rcen, theta=theta)
		self.B = np.sqrt(self.Bx**2 + self.By**2) 
		self.phi = np.arctan(self.Bx/self.By) 
		
		if density == None:
			self.ne = 1e-20 * np.ones_like(self.r)	# vanishing density 
		elif density == "churasov":
			self.ne = churasov_density(self.r)
		#self.rm = self.get_rm()

	def uniform_field_z(self, deltaL=1.0, Lmax=1800.0):
		self.r = np.arange(0, Lmax-deltaL, deltaL)
		self.deltaL = np.ones_like(self.r) * deltaL
		self.rcen = self.r + (0.5 * self.deltaL)
		self.Bx, self.By = 0.0, 0.0
		self.ne, self.Bz = self.profile(self.rcen)
		self.B = np.sqrt(self.Bx**2 + self.By**2) 
		self.phi = np.zeros_like(self.Bx)


	def get_rm(self):
		prefactor = (unit.e ** 3) / 2.0 / np.pi / unit.melec_csq / unit.melec_csq
		prefactor = 812.0

		# integrate using simpson 
		# units are cm^-3, kpc and microgauss
		#integral = simps(self.rcen * unit.kpc, self.ne * self.Bz * 1e6)
		#costheta = 
		#n_dot_Bz = 
		integral = simps(self.rcen, self.ne * self.Bz * 1e6)
		# ntegral = 1.0

		#integral = simps(self.rcen * unit.kpc, self.ne * self.Bz)
		
		# convert to rad / m^2 and return 
		return (prefactor * integral)

	def domain_from_slice(self, Cluster, deltaL=1.0, Lmax=500.0):
		self.r = np.arange(0, Lmax, deltaL)
		self.deltaL = np.ones_like(self.r) * deltaL
		self.rcen = self.r + (0.5 * self.deltaL)
		self.Bx, self.By = Cluster.get_B(self.rcen)
		self.Bz = Cluster.get_Bz(self.rcen)
		self.B = np.sqrt(self.Bx**2 + self.By**2) 
		self.phi = np.arctan(self.Bx/self.By) 
		self.ne = Cluster.density(self.r) 
		#self.rm = self.get_rm()

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
		rcen_array, r_array = [], []
		deltaL_array = []

		# wonder if there's a better way to do this?
		while r < L:
			# get a coherence length which will be the size of the box 
			# this can be a function or a float 
			if callable(coherence):
				lc = coherence()
			else:
				lc = coherence

			if self.coherence_r0 != None:
				lc *= (1.0 + (r/(self.coherence_r0)))

			# ensure the simulation is truncated at distance L
			#if (r + lc) > L:
			#	lc = (L-r) + 1e-10

			if rcen == r0:
				rcen += lc / 2.0
			else:
				rcen += lc

			rcen_array.append(rcen)
			r_array.append(r)
			deltaL_array.append(lc)

			r += lc

		# now we have box sizes and radii, get the field and density in each box 
		Ncells = len(r_array)
		self.r = np.array(r_array)
		self.rcen = np.array(rcen_array)
		self.deltaL = np.array(deltaL_array)

		# draw random isotropic angles and save phi
		theta, phi = random_angle(size = Ncells)
		#phi = phi

		# get density and magnetic field strength at centre of box
		self.ne, Btot = self.profile(self.r)
  
		# get the x and y components and increment r
		#Bx_array.append(B * np.sin(theta2))
		#y_array.append(B * np.cos(theta2))
		self.Bx = Btot * np.sin(theta) * np.cos(phi)
		self.By = Btot * np.sin(theta) * np.sin(phi)

		# note B is actually Bperp
		self.B = np.sqrt(self.Bx**2  + self.By **2)
		self.phi = np.arctan(self.Bx/self.By) 
		#self.phi = phi

		self.Bz = Btot * np.cos(theta) 
		self.rm = self.get_rm()
		#print (self.rm)







