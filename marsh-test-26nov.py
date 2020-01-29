import matplotlib.pyplot as plt 
import numpy as np 
from constants import *
import constants as c

UNIT_TIME  = 1519252268630104.8
UNIT_LENGTH= 50676.79373667135
UNIT_MASS  = 5.6095363761802584e+32
UNIT_GAUSS = 0.06925147467360344
HBAR_EV = 6.582119569e-16

def get_density(r):
	x1 = 3.9e-2 / ((1.0 + (r/80.0))**1.8)
	x2 = 4.05e-3 / ((1.0 + (r/280.0))**0.87)
	return (x1 + x2)

def get_B(r):
	x = r/(93.0)
	Bx = (0.00312443*(x**18)) - (0.0319991*(x**16)) + (0.260311*(x**14)) - (1.63197*(x**12)) + (7.58002*(x**10)) - (24.721*(x**8)) + (52.3929*(x**6)) - (63.8794*(x**4)) + (35.8973*(x**2)) - 5.86899	
	By = (0.0102459*(x**17))-(0.0937683*(x**15)) + (0.671841*(x**13)) - (3.6406*(x**11)) + (14.2479*(x**9)) - (37.7455*(x**7)) + (61.3611*(x**5)) - (51.7231*(x**3)) + (16.9128*x)
	return Bx, By


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


np.random.seed(12)

### Physical Parameters
mass = 10.0**(-12.7)
M = 1.0 / (1e-12) * 1e9
g = 1e-11 * 1e-9
ne = 1e-40
omega_pl = np.sqrt(4.0 * PI * E * E * ne / MELEC) * HBAR_EV
omega_pl = 1e-100
B = 1e-5 * UNIT_GAUSS
phi = 0.0
L = 1.0
distance = L * 1000.0 * PARSEC * UNIT_LENGTH

energy2 = np.logspace(8,11,1000)
import powerlaw, alpro 
# powerlaw PDF of coherence lengths 
distribution = powerlaw.Power_Law(xmin=3.5,xmax=10.0,parameters=[1.2])

# get Ainit 
Ainit = np.zeros( (len(energy2),6))
Ainit2 = np.zeros( (len(energy2),6))
# orientated in the y direction 
direction = 2

# Ainit[:,0] = 1.0
Ainit[:,direction] = 1.0
Ainit2[:,0] = 1.0
r = 0.0

Lmax = 93.0
Anew = Ainit
theta = np.zeros_like(energy2)

L  = 1
r = -L / 2.0
myL = 0.0
EPSILON = 1e-10



Probabilities = []

L  = 1
r = -L / 2.0
myL = 0.0
EPSILON = 1e-10

# g_a = 1e-11
density = 0.0
mass = 1e-9



g_as = np.logspace(-12,-11,10)
g_as = [1e-11]
for g_a in g_as:
	for mass in np.logspace(-9,-8,1):
		Ainit = np.zeros( (len(energy2),6))
		Ainit2 = np.zeros((len(energy2),6))
		energys = energy2
		Ainit[:,direction] = 1.0
		Ainit2[:,0] = 1.0
		myL = 0.0
		r = -L / 2.0
		while myL < (Lmax-EPSILON):
		# Bx, By, ne, r, x = generate_field_arrays(Lmax, distribution)	

		# for i in range(len(Bx)):
			# random length 
			myL += L 
			r += L 
			# L = x[i]
			#Anew = Ainit


			# ne = get_density(r)
			Bx,By = get_B(r)
			B = np.sqrt(Bx**2 + By**2)

			print (B)
			# phi = np.arctan(By/Bx)
			# B = np.sqrt(Bx**2 + By**2)
			# theta = np.arctan(Ainit[:,1]/Ainit[:,2])
			phi = (np.arctan(Bx/By) * np.ones_like(energys)) 
			phi2 = (np.arctan(Bx/By) * np.ones_like(energys)) 
			
			#phi = np.random.random(size=len(energy2)) * np.pi/2.0
			#print (Anew[0])
			B *= 1e-6 
			P1, Anew = alpro.get_P(energys, Ainit, phi, B, L, g_a * 1e-9, mass, 1e-20)
			P2, Anew2 = alpro.get_P(energys, Ainit2, phi2, B, L, g_a * 1e-9, mass, 1e-20)
			Ainit = Anew
			Ainit2 = Anew2
			# theta = -np.arctan(Anew[:,0]/Anew[:,2])

		P = 0.5 * (P1 + P2)

		# P = 0.5 * (P1 + P2)
		print (myL)
		plt.plot(energys/1e9, 1-P, lw=3, label="My Code", ls="-")

# plt.plot(energy2, Anew[:,4]**2, lw=3, label="Code Discretised", c="C3", ls=":")
#print (Anew)
x,y = np.loadtxt("marsh_data.txt", unpack=True)
plt.plot(x,y, label="Marsh Data", c="k", lw=3, ls="--")
x,y = np.loadtxt("libanov.dat", unpack=True)
plt.plot(x,y, label="Libanov Data", c="r", lw=2, alpha=0.8, ls="-")
plt.ylim(0,1.0)
plt.xlim(0.1,100)
plt.legend()
plt.ylabel(r"$P_{\gamma -> \gamma}$")
plt.semilogx()
plt.xlabel("Energy (GeV)")
#plt.semilogy()
plt.savefig("marsh_test.png", dpi=100)

plt.clf()
# x = np.arange(0,93,0.1)
# plt.plot(x, get_B(x)[0])
# plt.plot(x, get_B(x)[1])
# plt.show()
# 