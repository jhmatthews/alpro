import matplotlib.pyplot as plt 
import numpy as np 
from constants import *
import constants as c
import alpro 
import matplotlib
import sys

UNIT_TIME  = 1519252268630104.8
UNIT_LENGTH= 50676.79373667135
UNIT_MASS  = 5.6095363761802584e+32
UNIT_GAUSS = 0.06925147467360344
HBAR_EV = 6.582119569e-16
alpro.util.set_default_plot_params() 

def make_3panel_plot(E, P, Pnorm, Bfield_model, L, Lmax, savename):

	# plot field
	z = np.linspace(0, Lmax, 256)

	plt.figure(figsize=(6,8))
	plt.subplot(311)
	Bx, By = Bfield_model.get_B(z)
	plt.plot(z, Bx, lw=3)
	plt.plot(z, By, lw=3)

	# plot coherence points 
	zz = np.arange(0,Lmax,L)
	Bx, By = Bfield_model.get_B(zz)
	BB = np.sqrt(Bx**2 + By**2)
	plt.scatter(zz, Bx, alpha=1)
	plt.scatter(zz, By, alpha=1)
	plt.xlabel("$z$")

	# plot survival curve 
	plt.subplot(312)
	plt.plot(E, P, lw=2, label=str(L/lc), ls="--", alpha=1, color="k")
	plt.plot(E, Pnorm, lw=3, label=str(L/lc), ls="-", alpha=0.5, color="k")
	plt.ylim(0.5,1)
	plt.semilogx()
	plt.xlabel("$E$")
	plt.ylabel("Probability")

	# plot survival curve 
	plt.subplot(313)
	plt.plot(E, (P-Pnorm)/Pnorm, lw=2, label=str(L/lc), ls="-", alpha=1, color="k")
	plt.ylim(-0.3,0.3)
	plt.semilogx()
	plt.xlabel("$E$")
	plt.ylabel("Relative Error")

	plt.savefig("prop_curves/3panel_{}.png".format(savename), dpi=300)
	plt.clf()



def get_density(r):
	x1 = 3.9e-2 / ((1.0 + (r/80.0))**1.8)
	x2 = 4.05e-3 / ((1.0 + (r/280.0))**0.87)
	return (x1 + x2)

class Bfield:
	def __init__(self, fname="Bfield.npy"):
		# load the array
		self.Bfull = np.load(fname)

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
		kind = "quadratic"
		# kind='quadratic'
		self.interp_x = interp1d(self.z, self.B[:,0], kind=kind)
		self.interp_y = interp1d(self.z, self.B[:,1], kind=kind)


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
Bfield_model = Bfield(fname="Bfield_0.npy")

DEGRADE_FACTOR = 1

# define global arrays and variables 
mode = sys.argv[1]
z = np.linspace(0,100,1024)

if mode == "lc":
	lc = 12.471
	Ls = np.array([0.01,0.04,0.1,0.2,0.5,1,2]) * lc
elif mode == "lgrid":
	lc = 0.78/2.0 * DEGRADE_FACTOR
	Ls = np.array([0.05,0.1,0.5,1,2]) * lc
	Ls = np.linspace(0.1,2,20) * lc

#Ls = np.linspace(0.039,25.0,100.0)
Ls = np.logspace(-1.5,1.4,100)

Lmax = 100.0 - np.max(Ls)
Lmax = 70.0
energy2 = np.logspace(3,7,1000)

# different slices
slices_to_take = ((0,0),(25,25),(50,50),(75,75),(100,100))
slices_to_take = ((25,25),)


for islice,slices in enumerate(slices_to_take):
	print ("Slice", islice)

	# take a slice along the B field
	Bfield_model.slice(1e-5, slices=slices, degrade=DEGRADE_FACTOR)

	Bx, By = Bfield_model.get_B(z)
	Bfield1d_plot(z, Bfield_model, slices[0], Lmax, lc)

	theta = np.zeros_like(energy2)

	probabilities = np.zeros((len(Ls),len(energy2)))

	my_cmap = matplotlib.cm.get_cmap('viridis')
	colors = my_cmap(np.linspace(0,1,num=len(Ls)+2))

	for il, L in enumerate(Ls):
		print (il, len(Ls))
		EPSILON = 1e-6

		# g_a = 1e-11
		#density = 0.0
		#mass = 1e-9

		g_as = np.logspace(-13,-11,10)
		g_as = [5e-13]

		for g_a in g_as:
			for mass in np.logspace(-13,-8,1):

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
					P1, Anew = alpro.get_P(energys, Ainit, phi, B, L, g_a * 1e-9, mass, 1e-20)
					P2, Anew2 = alpro.get_P(energys, Ainit2, phi2, B, L, g_a * 1e-9, mass, 1e-20)
					
					# record new initial vectors 
					Ainit = Anew
					Ainit2 = Anew2

				# average the probabilities
				P = 0.5 * (P1 + P2)

				if il == 0:
					one_minus_P0 = 1-P

				# plot the result 
				#print (myL, myL-L, Lmax-EPSILON, Lmax, 100 - np.max(Ls))
				else:
					#savename = "smode+"_"+str(il)
					savename = "s{}_Lmax{}_{}{:03d}".format(slices[0], Lmax, mode, il)
					#savename = "s{}_deg{}_Lmax{}_{}{:.1f}".format(slices[0], DEGRADE_FACTOR, Lmax, mode, L/lc)
					make_3panel_plot(energys/1e3, 1-P, one_minus_P0, Bfield_model, L, Lmax, savename)
					#make_3panel_plot(E, P, Pnorm, Bfield_model, L, Lmax, savename)
					plt.close("all")

					probabilities[il] = 1-P

plt.close("all")
def L_residual_plot(lengths, resid, fname="residual_v_length", ylabel="Residual", 
	                lims=(-0.2,0.2)):
	plt.figure(figsize=(7,5))

	plt.plot(lengths, resid, ls="steps")
	plt.xlabel("Resolution (kpc)", fontsize=16)
	plt.ylabel(ylabel, fontsize=16)
	plt.xlim(0.039,25.0)
	plt.ylim(lims[0], lims[1])
	plt.gca().vlines([0.39,12.471], lims[0],lims[1], color="k", ls="--", lw=2)
	plt.text(14,-0.1,"$l_c$", fontsize=14)
	plt.text(0.43,-0.1,r"$l_{{ \rm grid}}$", fontsize=14)
	plt.semilogx()
	plt.savefig(fname+"png", dpi=200)


from scipy.interpolate import interp1d
interped_probs = np.zeros_like(probabilities)
E_linear = np.linspace(1e3,1e4,100)
E_linear = np.logspace(3,4,100)
interped_probs = np.zeros( (len(Ls),len(E_linear)) )
for i in range(len(Ls)):
	interp_func = interp1d(energys, probabilities[i])
	interped_probs[i,:] = interp_func(E_linear)

probabilities_use = interped_probs

one_minus_P0 = probabilities_use[1] 
resid = (probabilities_use - one_minus_P0) / one_minus_P0
mean_resid = np.mean(resid, axis=1) 
median_resid = np.median(resid, axis=1) 
max_resid = np.max(resid, axis=1) 
mean_P = np.mean(probabilities, axis=1)

L_residual_plot(Ls, mean_resid, fname="residual_l_mean", ylabel="Mean Residual")
L_residual_plot(Ls, median_resid, fname="residual_l_median", ylabel="Median Residual")
#L_residual_plot(Ls, max_resid, fname="residual_l_max", ylabel="Max Residual")
L_residual_plot(Ls, mean_P, fname="prob_l_mean", ylabel=r"Mean $P_{\gamma \rightarrow \gamma}$", lims=(0.5,1))

plt.figure(figsize=(7,5))
plt.pcolormesh(np.log10(E_linear), Ls, probabilities_use, vmin=0.7)
cbar = plt.colorbar()
plt.semilogy()
plt.ylim(0.05,25.0)
plt.gca().hlines([0.39,12.471], 3,4, color="C3", ls="--", lw=2)
plt.text(14,-0.1,"$l_c$", fontsize=14, color="C3")
plt.text(0.43,-0.1,r"$l_{{ \rm grid}}$", fontsize=14, color="C3")
cbar.set_label(r"$P_{\gamma \rightarrow \gamma}$")
plt.xlabel(r"$\log E$ (eV)", fontsize=16)
plt.ylabel("Resolution (kpc)", fontsize=16)
plt.savefig("heatmap.png")

