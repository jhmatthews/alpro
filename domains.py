import matplotlib.pyplot as plt 
import numpy as np 
from constants import *
import constants as c

def get_density(r):
	x1 = 3.9e-2 / ((1.0 + (r/80.0))**1.8)
	x2 = 4.05e-3 / ((1.0 + (r/280.0))**0.87)
	return (x1 + x2)

def get_B(r, ne):
	B = 25.0 * 1e-6 * (ne / get_density(0.0))
	return B

def generate_field_arrays(domain_size, distribution):
	x = []
	L = 0.0
	while L < domain_size: 
		size = distribution.generate_random()
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
UNIT_TIME  = 1519252268630104.8
UNIT_LENGTH= 50676.79373667135
UNIT_MASS  = 5.6095363761802584e+32
UNIT_GAUSS = 0.06925147467360344
HBAR_EV = 6.582119569e-16

### Physical Parameters
mass = 10.0**(-12.7)
M = 1.0 / (1e-12) * 1e9
g = 1e-12 * 1e-9
ne = 0.1
omega_pl = np.sqrt(4.0 * PI * E * E * ne / MELEC) * HBAR_EV
B = 1e-5 * UNIT_GAUSS
phi = 0.0
L = 1.0
distance = L * 1000.0 * PARSEC * UNIT_LENGTH

energy2 = np.logspace(3,4,100)
import powerlaw, alpro 
# powerlaw PDF of coherence lengths 
distribution = powerlaw.Power_Law(xmin=3.5,xmax=10.0,parameters=[1.2])

# get Ainit 
Ainit = np.zeros( (len(energy2),6))
# orientated in the y direction 
direction = 2

# Ainit[:,0] = 1.0
Ainit[:,direction] = 1.0
r = 0.0

Lmax = 1000.0

while r < Lmax:
# Bx, By, ne, r, x = generate_field_arrays(Lmax, distribution)	

# for i in range(len(Bx)):
	# random length 
	L = distribution.generate_random()
	# L = x[i]

	r += L / 2.0
	ne = get_density(r)
	B = get_B(r, ne)
	# B = np.sqrt(Bx**2 + By**2)
	# phi = np.arctan(Bx/By)


	omega_pl = np.sqrt(4.0 * PI * E * E * ne / MELEC) * HBAR_EV
	phi = np.random.random(size=len(energy2)) * np.pi/2.0

	P, Anew = alpro.get_P(energy2, Ainit, phi, B, L, g, mass, ne)
	Ainit = Anew


plt.plot(energy2, 1-P, lw=3, label="Code Discretised", c="C0", ls="-")
# plt.plot(energy2, Anew[:,4]**2, lw=3, label="Code Discretised", c="C3", ls=":")
print (Anew)
plt.ylim(0.75,1.0)
plt.legend()
plt.ylabel(r"$P_{\gamma -> a}$")
plt.semilogx()
plt.xlabel("Energy (eV)")
# plt.semilogy()
plt.savefig("domain_test.png", dpi=100)
# 