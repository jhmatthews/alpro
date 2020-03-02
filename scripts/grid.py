import matplotlib.pyplot as plt 
import numpy as np 
from constants import *
import constants as c
import alpro 
import matplotlib
import sys
import time

UNIT_TIME  = 1519252268630104.8
UNIT_LENGTH= 50676.79373667135
UNIT_MASS  = 5.6095363761802584e+32
UNIT_GAUSS = 0.06925147467360344
HBAR_EV = 6.582119569e-16

def get_density_perseus(r):
	x1 = 3.9e-2 / ((1.0 + (r/80.0))**1.8)
	x2 = 4.05e-3 / ((1.0 + (r/280.0))**0.87)
	return (x1 + x2)

class Bfield:
	def __init__(self, field_type="simple", fname="Bfield.npy", seed=0):
		
		self.type = field_type

		if field_type == "load":
			# load the array
			self.Bfull = np.load("Bfield_{}.npy".format(seed))
		elif field_type == "gen":
			#some code to generate a Bfield
		elif field_type == ""


	def slice(self, B0, size=100, slices=(100,100),degrade=1):
		# take a slice along the B field
		self.B = self.Bfull[:,slices[0],slices[1],:]
		Bstrength = np.linalg.norm(self.B, axis=1)
		Brms = np.mean(np.sqrt(Bstrength**2))
		normalisation = B0 / Brms
		self.B *= normalisation

		from scipy.interpolate import interp1d
		ztrue = np.linspace(0,size,len(self.B[:,0]))
		self.z = np.linspace(0,size,len(self.B[:,0])/degrade)

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

		# actual interpolation always done using 2nd order interpolation  
		self.interp_x = interp1d(self.z, self.B[:,0], kind='quadratic')
		self.interp_y = interp1d(self.z, self.B[:,1], kind='quadratic')


	def get_B(self, z, B0=None):
		'''
		get the magnetic field at distance z 
		'''
		if self.type == "simple":
			if B0 == None:
				raise AttributeError("No B0 kwarg specified!")
			costheta = 2.0 * np.random.random() - 1.0
			theta = np.arccos(costheta)
			phi = np.random.random() * 2.0 * np.pi
			Bx = B0 * np.sin(theta) * np.cos(phi)
			By = B0 * np.sin(theta) * np.cos(phi)
		else:	
			Bx = self.interp_x(z)
			By = self.interp_y(z)
		return (Bx, By)


def Bfield1d_plot(z, Bfield_model, i, Lmax, lc):
	# plot field
	z = z[z<Lmax]
	Bx, By = Bfield_model.get_B(z)
	plt.plot(z, Bx, lw=3)
	plt.plot(z, By, lw=3)

	# plot coherence points 
	zz = np.arange(0,Lmax,lc)
	Bx, By = Bfield_model.get_B(zz)
	BB = np.sqrt(Bx**2 + By**2)
	plt.scatter(zz, Bx, alpha=1)
	plt.scatter(zz, By, alpha=1)

	# label and save 
	plt.ylabel("r (kpc)")
	plt.savefig("1d_fields/Bfield1d_{}.png".format(i))
	plt.clf()

def propagate(Bfield, Lmax, L):
	'''
	Propagate an ALP through a given Bfield model

	Parameters:
		Bfield 		object
					Bfield class that contains the Bfield model. 
					gets queried for B at radius r.
		Lmax 		float
					maximum radius in kpc 
		L 			float
					resolution / box size in kpc 
	'''
	Lsteps = np.arange(0,Lmax,L)
	rsteps = np.arange(L/2,Lmax-(L/2),L)

	for i, r in enumerate(rsteps):

		# get magnetic field
		ne = get_density(r)
		Bx, By = Bfield_model.get_B(r)
		B = np.sqrt(Bx**2 + By**2)

		Bx, By = Bfield_model.get_B(r)
		B = np.sqrt(Bx**2 + By**2)

		# get angle of rotation 
		phi = (np.arctan(Bx/By) * np.ones_like(energys)) 
		phi2 = (np.arctan(Bx/By) * np.ones_like(energys)) 
		
		# get the axion probability for both polarisations 
		P1, Anew = alpro.get_P(energys, Ainit, phi, B, L, g_a * 1e-9, mass, 1e-20)
		P2, Anew2 = alpro.get_P(energys, Ainit2, phi2, B, L, g_a * 1e-9, mass, 1e-20)

		# record new initial vectors 
		Ainit = Anew
		Ainit2 = Anew2

	# average the probabilities
	P = 0.5 * (P1 + P2)

	return (P)



np.random.seed(12)

# initialise B field model 
Bfield_model = Bfield(1e-5)

DEGRADE_FACTOR = 1

# define global arrays and variables 
z = np.linspace(0,100,1024)

# sample on coherence length and grid length
Ls = np.array([0.78/2.0 * 0.2, 12.471])

Lmax = 100.0
energy2 = np.logspace(3,7,1000)

# different slices
slices_to_take = ((0,0),(25,25),(50,50),(75,75),(100,100))
#slices_to_take = ((0,0),(25,25))

# parameters - 1 value at the moment 
masses = np.logspace(-11,-8,1)
g_as = np.logspace(-12,-11,1)

# shape of arrays to store 
NSEEDS = 10
random_seeds = np.arange(NSEEDS)
shape = (len(slices_to_take),NSEEDS,len(energy2))
residuals = np.zeros(shape)

for iseed,seed in enumerate(random_seeds):
	# initialise B field model 
	Bfield_model = Bfield(1e-5, seed=seed)
	print ("Bfield seed:", seed)


	for islice,slices in enumerate(slices_to_take):
		print ("Slice", islice)

		# take a slice along the B field
		Bfield_model.slice(1e-5, slices=slices, degrade=DEGRADE_FACTOR)

		Bx, By = Bfield_model.get_B(z)
		#Bfield1d_plot(z, Bfield_model, slices[0], Lmax, lc)

		theta = np.zeros_like(energy2)

		Probabilities = []

		my_cmap = matplotlib.cm.get_cmap('viridis')
		colors = my_cmap(np.linspace(0,1,num=len(Ls)+2))

		EPSILON = 1e-6

		g_as = np.logspace(-12,-11,10)
		g_as = [1e-11]

		for g_a in g_as:
			for mass in masses:
				for il, L in enumerate(Ls):

					# create the initial state vector 
					Ainit = np.zeros( (len(energy2),6))
					Ainit2 = np.zeros((len(energy2),6))
					energys = energy2
					# set the X and Y polarisations 
					Ainit[:,2] = 1.0
					Ainit2[:,0] = 1.0

					# these are the steps we are going to take in L
					P = propagate(Bfield, Lmax, L)

					if il == 0:
						one_minus_P0 = 1-P
						fine__probability
					else:
						residuals[islice,iseed,:] = ((1-P)-one_minus_P0)/one_minus_P0

				# now do random boxes approach



np.save("residual_array.npy", residuals)


residuals = np.load("residual_array.npy")

# for iseed,seed in enumerate(random_seeds):
# 	for islice,slices in enumerate(slices_to_take):
# 		plt.plot(energy2/1e3,residuals[islice,iseed,:], alpha=0.5)

for iseed,seed in enumerate(random_seeds):
	x = np.mean(residuals[:,iseed,:], axis=0)
	plt.plot(energy2/1e3,x, alpha=1)


#plt.ylim(-0.3,0.3)
plt.xlim(1,1e4)
plt.legend()
plt.ylabel(r"$P_{\gamma -> \gamma}$ Residual")
plt.semilogx()
plt.xlabel("Energy (keV)")
#plt.semilogy()
plt.show()
#plt.savefig("prop_curves/propagation_s{}_deg{}_{}.png".format(slices[0],DEGRADE_FACTOR,mode), dpi=300)

plt.clf()

