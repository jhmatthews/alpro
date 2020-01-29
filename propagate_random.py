import matplotlib.pyplot as plt 
import numpy as np 
from constants import *
import constants as c
import alpro 
import matplotlib

UNIT_TIME  = 1519252268630104.8
UNIT_LENGTH= 50676.79373667135
UNIT_MASS  = 5.6095363761802584e+32
UNIT_GAUSS = 0.06925147467360344
HBAR_EV = 6.582119569e-16

def get_density(r):
	x1 = 3.9e-2 / ((1.0 + (r/80.0))**1.8)
	x2 = 4.05e-3 / ((1.0 + (r/280.0))**0.87)
	return (x1 + x2)

class Bfield:
	def __init__(self, fname="Bfield.npy"):
		# load the array
		self.Bfull = np.load("Bfield.npy")

	def slice(self, B0, size=100, slices=(100,100)):
		# take a slice along the B field
		self.B = self.Bfull[:,slices[0],slices[1],:]
		Bstrength = np.linalg.norm(self.B, axis=1)
		Brms = np.mean(np.sqrt(Bstrength**2))
		normalisation = B0 / Brms
		self.B *= normalisation
		print (self.B.shape)
		#print (Brms, np.sqrt(self.B[0,:]**2))

		from scipy.interpolate import interp1d
		self.z = np.linspace(0,size,len(self.B[:,0]))
		#rint (self.B[:,0].shape)
		#import sys
		#sys.exit()
		print (self.B[:,0], self.B[:,1])

		self.interp_x = interp1d(self.z, self.B[:,0])
		self.interp_y = interp1d(self.z, self.B[:,1])

	def get_B(self, z):
		Bx = self.interp_x(z)
		By = self.interp_y(z)
		return (Bx, By)

def generate_field_arrays(domain_size, distribution):
	x = []
	L = 0.0
	while L < domain_size: 
		size = 0.1 * PARSEC * 1000.0
		x.append(size)
		L += size

	x = np.array(x)
	r = np.cumsum(x)-(x/2.0)
	ne = get_density(r)
	B = get_B(r,ne)
	phi = np.random.random(size=len(B)) * np.pi / 2.0
	Bx = np.sin(phi)
	By = np.cos(phi)
	return (Bx, By, ne, r, x)

def Bfield1d_plot(z, Bfield_model, i):
	# plot field
	Bx, By = Bfield_model.get_B(z)
	plt.plot(z, Bx, lw=3)
	plt.plot(z, By, lw=3)

	# plot coherence points 
	zz = np.arange(0,100,lc)
	Bx, By = Bfield_model.get_B(zz)
	BB = np.sqrt(Bx**2 + By**2)
	plt.scatter(zz, Bx, alpha=1)
	plt.scatter(zz, By, alpha=1)

	# label and save 
	plt.ylabel("r (kpc)")
	plt.savefig("1d_fields/Bfield1d_{}.png".format(i))
	plt.clf()

np.random.seed(12)

# initialise B field model 
Bfield_model = Bfield(1e-5)

# define global arrays and variables 
z = np.linspace(0,100,256)
lc = 22.98
Ls = np.array([0.01,0.1,0.2,0.5,1,2]) * lc
Lmax = 100.0 - np.max(Ls)
energy2 = np.logspace(3,7,1000)

# different slices
slices_to_take = ((0,0),(25,25),(50,50),(75,75),(100,100))

for islice,slices in enumerate(slices_to_take):

	# take a slice along the B field
	Bfield_model.slice(1e-5, slices=slices)

	Bx, By = Bfield_model.get_B(z)
	Bfield1d_plot(z, Bfield_model, slices[0])

	theta = np.zeros_like(energy2)

	Probabilities = []

	my_cmap = matplotlib.cm.get_cmap('viridis')
	colors = my_cmap(np.linspace(0,1,num=len(Ls)))

	for il, L in enumerate(Ls):
		r = -L / 2.0
		myL = 0.0
		EPSILON = 1e-10

		# g_a = 1e-11
		density = 0.0
		mass = 1e-9

		fudge_factor = 1.0

		g_as = np.logspace(-12,-11,10)
		g_as = [1e-11]

		for g_a in g_as:
			for mass in np.logspace(-11,-8,1):
				Ainit = np.zeros( (len(energy2),6))
				Ainit2 = np.zeros((len(energy2),6))
				energys = energy2
				Ainit[:,2] = 1.0
				Ainit2[:,0] = 1.0
				myL = 0.0
				r = -L / 2.0
				while myL < (Lmax-EPSILON):

					myL += L 
					r += L 

					print (myL, r)

					ne = get_density(r)
					Bx, By = Bfield_model.get_B(r)
					B = np.sqrt(Bx**2 + By**2)

					#print (B, Bx, By)
					# phi = np.arctan(By/Bx)
					# B = np.sqrt(Bx**2 + By**2)
					# theta = np.arctan(Ainit[:,1]/Ainit[:,2])
					phi = (np.arctan(Bx/By) * np.ones_like(energys)) 
					phi2 = (np.arctan(Bx/By) * np.ones_like(energys)) 
					
					#phi = np.random.random(size=len(energy2)) * np.pi/2.0
					#print (Anew[0])
					#B *= 1e-6 / fudge_factor
					P1, Anew = alpro.get_P(energys, Ainit, phi, B, L, g_a * 1e-9, mass, 1e-20)
					P2, Anew2 = alpro.get_P(energys, Ainit2, phi2, B, L, g_a * 1e-9, mass, 1e-20)
					Ainit = Anew
					Ainit2 = Anew2
					# theta = -np.arctan(Anew[:,0]/Anew[:,2])

				P = 0.5 * (P1 + P2)

				# P = 0.5 * (P1 + P2)
				print (myL)
				plt.plot(energys/1e3, 1-P, lw=3, label=str(L/lc), ls="-", alpha=0.7, color=colors[il])

	# plt.plot(energy2, Anew[:,4]**2, lw=3, label="Code Discretised", c="C3", ls=":")
	#print (Anew)
	# x,y = np.loadtxt("marsh_data.txt", unpack=True)
	# plt.plot(x,y, label="Marsh Data", c="k", lw=3, ls="--")
	plt.ylim(0,1.1)
	plt.xlim(1,1e4)
	plt.legend()
	plt.ylabel(r"$P_{\gamma -> \gamma}$")
	plt.semilogx()
	plt.xlabel("Energy (keV)")
	#plt.semilogy()
	plt.savefig("prop_curves/propagation_{}.png".format(slices[0]), dpi=300)

	plt.clf()
# x = np.arange(0,93,0.1)
# plt.plot(x, get_B(x)[0])
# plt.plot(x, get_B(x)[1])
# plt.show()
# 