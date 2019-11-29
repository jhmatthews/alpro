import matplotlib.pyplot as plt 
import numpy as np 
from constants import *
import constants as c

UNIT_TIME  = 1519252268630104.8
UNIT_LENGTH= 50676.79373667135
UNIT_MASS  = 5.6095363761802584e+32
UNIT_GAUSS = 0.06925147467360344
HBAR_EV = 6.582119569e-16

### Physical Parameters
mass = 1e-12
M = 1.0 / (1e-12) * 1e9
g = 1e-12 * 1e-9
ne = 0.1
omega_pl = np.sqrt(4.0 * PI * E * E * ne / MELEC) * HBAR_EV
B = 1e-5 * UNIT_GAUSS
phi = 0.0
L = 1.0
distance = L * 1000.0 * PARSEC * UNIT_LENGTH

energy2 = np.logspace(3,4,1000)

import alpro
Ainit = np.zeros( (len(energy2),6))

# orientated in the y direction 
direction = 2

# Ainit[:,0] = 1.0
Ainit[:,direction] = 1.0
phi = np.ones_like(energy2) * np.pi/2 * 0.0
P, Anew = alpro.get_P(energy2, Ainit, phi, B / UNIT_GAUSS, L, g, mass, ne)

plt.plot(energy2, P, c="C0", lw=3, label="Code")
# print (Ainit, Anew)
PP = 1.0
Ainit = np.zeros( (len(energy2),6))
Ainit[:,direction] = 1.0
Anew = Ainit
Nshells = 10


for i in range(Nshells):
	#Anew = Ainit
	# phi = np.arctan(Anew[:,0]/Anew[:,1])
	# phi = np.pi/2
	# Ainit[:,2] = Anew[:,2]
	Ainit = Anew
	
	P2, Anew = alpro.get_P(energy2, Ainit, phi, B / UNIT_GAUSS, L/Nshells, g, mass, ne)
	print (Anew)
	# phi = np.pi/2.0 * np.random.random(len(phi))
	# PP *= (1.0 - P2)

plt.plot(energy2, 1-P2, lw=3, label="Code Discretised", c="C2", ls="--")
# plt.plot(energy2, Anew[:,4]**2, lw=3, label="Code Discretised", c="C3", ls=":")
print (Anew)

plt.legend()
plt.ylabel(r"$P_{\gamma -> a}$")
plt.semilogx()
plt.xlabel("Energy (eV)")
plt.semilogy()
plt.savefig("discretised_test.png", dpi=100)

plt.clf()


Nshells = 2
Ainit = np.zeros( (len(energy2),6))
Ainit[:,direction] = 1.0
Anew = Ainit
for i in range(Nshells):

	Ainit = Anew
	
	P2, Anew = alpro.get_P(energy2, Ainit, phi, B / UNIT_GAUSS, L/Nshells, g, mass, ne)
	# phi += np.pi/2.0

plt.plot(energy2, P2, lw=3, label="With B Flip", c="C2", ls="-")
print (Anew)
plt.legend()
plt.ylabel(r"$P_{\gamma -> a}$")
plt.semilogx()
plt.xlabel("Energy (eV)")
# plt.semilogy()
plt.savefig("flip_test.png", dpi=100)

