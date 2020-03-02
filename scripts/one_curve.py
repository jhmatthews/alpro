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
plt.rcParams["text.usetex"] = "True"
plt.rcParams['xtick.labelsize']=14
plt.rcParams['ytick.labelsize']=14
plt.rcParams['font.serif']=['cm']
plt.rcParams['font.family']='serif'	
plt.rcParams["text.usetex"] = "True"
plt.rcParams["lines.linewidth"] = 3
plt.rcParams["axes.linewidth"] = 1.5
plt.rcParams["xtick.major.width"] = 1.5
plt.rcParams["ytick.major.width"] = 1.5

def get_density(r):
	x1 = 3.9e-2 / ((1.0 + (r/80.0))**1.8)
	x2 = 4.05e-3 / ((1.0 + (r/280.0))**0.87)
	return (x1 + x2)

class Bfield:
	def __init__(self, fname="Bfield.npy"):
		# load the array
		self.Bfull = np.load("Bfield.npy")

	def slice(self, B0, size=100, slices=(100,100),degrade=1):
		# take a slice along the B field
		self.B = self.Bfull[:,slices[0],slices[1],:]
		Bstrength = np.linalg.norm(self.B, axis=1)
		Brms = np.mean(np.sqrt(Bstrength**2))
		normalisation = B0 / Brms
		self.B *= normalisation
		#print (self.B.shape)
		#print (Brms, np.sqrt(self.B[0,:]**2))

		from scipy.interpolate import interp1d
		ztrue = np.linspace(0,size,len(self.B[:,0]))
		self.z = np.linspace(0,size,len(self.B[:,0])/degrade)
		#rint (self.B[:,0].shape)
		#import sys
		#sys.exit()
		#print (self.B[:,0], self.B[:,1])
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
		self.interp_x = interp1d(self.z, self.B[:,0], kind='quadratic')
		self.interp_y = interp1d(self.z, self.B[:,1], kind='quadratic')


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

np.random.seed(12)

# initialise B field model 
Bfield_model = Bfield(1e-5)

DEGRADE_FACTOR = 1

# define global arrays and variables 
mode = "lgrid"
z = np.linspace(0,100,1024)


lc = 0.78/2.0 * DEGRADE_FACTOR
Ls = np.array([0.05]) * lc

Lmax = 10.0
energy2 = np.logspace(3,7,1000)

# different slices
slices_to_take = ((50,50),)

for islice,slices in enumerate(slices_to_take):
	print ("Slice", islice)

	# take a slice along the B field
	Bfield_model.slice(1e-5, slices=slices, degrade=DEGRADE_FACTOR)

	Bx, By = Bfield_model.get_B(z)
	Bfield1d_plot(z, Bfield_model, slices[0], Lmax, lc)

	theta = np.zeros_like(energy2)

	Probabilities = []

	my_cmap = matplotlib.cm.get_cmap('viridis')
	colors = my_cmap(np.linspace(0,1,num=len(Ls)+2))

	plt.figure(figsize=(6,5))

	for il, L in enumerate(Ls):
		EPSILON = 1e-6

		# g_a = 1e-11
		#density = 0.0
		#mass = 1e-9

		g_as = np.logspace(-12,-11,3)
		g_as = [1e-12]

		for ig, g_a in enumerate(g_as):
			for mass in np.logspace(-11,-8,1):

				# create the initial state vector 
				Ainit = np.zeros( (len(energy2),6))
				Ainit2 = np.zeros((len(energy2),6))
				energys = energy2
				Ainit[:,2] = 1.0
				Ainit2[:,0] = 1.0

				# these are the steps we are going to take in L
				Lsteps = np.arange(0,Lmax,L)
				rsteps = np.arange(L/2,Lmax-(L/2),L)

				for i, r in enumerate(rsteps):

					# get magnetic field
					ne = get_density(r)
					Bx, By = Bfield_model.get_B(r)
					B = np.sqrt(Bx**2 + By**2)

					# get angle of rotation 
					phi = (np.arctan(Bx/By) * np.ones_like(energys)) 
					phi2 = (np.arctan(Bx/By) * np.ones_like(energys)) 
					
					# get the axion probability for both polarisations 
					t1 = time.time()
					P1, Anew = alpro.get_P(energys, Ainit, phi, B, L, g_a * 1e-9, mass, 1e-2)
					P2, Anew2 = alpro.get_P(energys, Ainit2, phi2, B, L, g_a * 1e-9, mass, 1e-2)

					print (time.time() - t1)
					
					# record new initial vectors 
					Ainit = Anew
					Ainit2 = Anew2

				# average the probabilities
				P = 0.5 * (P1 + P2)


				plt.plot(energys/1e3, (1-P), lw=3, label=str(L/lc), ls="-", alpha=0.7, color=colors[ig])


	#plt.title(mode)
	plt.text(2,0.9984,r"$g_a = 10^{-12}~\mathrm{GeV}^{-1}, m_a = 10^{-11}~\mathrm{eV}$", fontsize=16)
	#plt.ylim(-0.3,0.3)
	plt.xlim(1,1e2)
	#plt.legend()
	plt.ylabel(r"$P_{\gamma -> \gamma}$", fontsize=16)
	plt.semilogx()
	plt.xlabel("Energy (keV)", fontsize=16)
	#plt.semilogy()
	plt.savefig("onecurve.png".format(slices[0],DEGRADE_FACTOR,mode), dpi=300)

	plt.clf()

