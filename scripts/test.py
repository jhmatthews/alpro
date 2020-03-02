import matplotlib.pyplot as plt 
import numpy as np 
from constants import *
import constants as c

def mod_squared(vector, i):
	'''
	get mod squared, ith element, for 3 complex vector
	''' 
	# if vector.shape == (6,):
	z = np.complex(vector[i], vector[i+3])
	# else:
		# z = [mod_squared(Anew[j], 2) for i in range(len(Anew))]
	return np.abs(z)**2

def get_prediction(energy, mass, M, omega_pl, B, distance):
	Delta_PL = -(omega_pl * omega_pl) / 2.0 / energy
	Delta_AA = -(mass * mass) / 2.0 / energy
	Delta_AG = B / 2.0 / M
	x = (Delta_PL - Delta_AA) ** 2
	Delta_OSC = np.sqrt (x + (4.0 * Delta_AG * Delta_AG))
	#Delta_OSC = np.sqrt (x )
	#print (x, Delta_AG)
	alpha = 0.5 * np.arctan(2.0 * Delta_AG / (Delta_PL - Delta_AA))

	term1 = np.sin(2.0*alpha)**2
	term2 = np.sin(Delta_OSC * distance / 2.0)**2
	prediction = term1 * term2

	return (prediction, term2)

def get_prediction2(energy, mass, M, omega_pl, B, distance):
	Delta_PL = -(omega_pl * omega_pl) / 2.0 / energy
	Delta_AA = -(mass * mass) / 2.0 / energy
	Delta_AG = B / 2.0 / M
	x = (Delta_PL - Delta_AA) ** 2
	Delta_OSC = np.sqrt (x + (4.0 * Delta_AG * Delta_AG))
	#Delta_OSC = np.sqrt (x )
	#print (x, Delta_AG)
	alpha = 0.5 * np.arctan(2.0 * Delta_AG / (Delta_PL - Delta_AA))

	term1 = (B / M / Delta_OSC)**2
	term2 = np.sin(Delta_OSC * distance / 2.0)**2
	prediction = term1 * term2

	return (prediction, term1)


energy2 = np.logspace(3,4,1000)


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


# prediction, alpha = get_prediction2 (energy2, mass, M, omega_pl, B, distance)
# prediction, alpha = get_prediction (energy2, mass, M, omega_pl, B, distance)

import alpro
Ainit = np.zeros( (len(energy2),6))
Ainit[:,0] = 1.0
phi = np.zeros_like(energy2) * np.pi/2.0
P, Anew = alpro.get_P(energy2, Ainit, phi, B / UNIT_GAUSS, L, g, mass, ne)

plt.plot(energy2, P, c="C0", lw=3, label="Code")
# print (Ainit, Anew)
PP = 1.0
Ainit = np.zeros( (len(energy2),6))
Ainit[:,0] = 1.0
Anew = Ainit
Nshells = 100
for i in range(Nshells):
	#Anew = Ainit
	# phi = np.arctan(Anew[:,0]/Anew[:,1])
	# phi = np.pi/2
	Ainit = Anew
	
	P2, Anew = alpro.get_P(energy2, Ainit, phi, B / UNIT_GAUSS, L/Nshells, g, mass, ne)

	PP *= (1.0 - P2)
plt.plot(energy2, P2, lw=3, label="Code Discretised", c="C2", ls="--")

print (P)

pp = [mod_squared(Anew[i], 2) for i in range(len(Anew))]

print (pp)
print (Anew)
# plt.plot(energy2, pp, c="r", lw=5, label = "From vectors", alpha=0.5)


theta = 2.0 * B * energy2 / M / (mass**2 - omega_pl**2)
Delta = (omega_pl**2 - mass**2) / 4.0 / energy2
prob = theta**2 / (theta**2+1) * (np.sin( Delta * distance * np.sqrt(1.0 + theta**2)))**2
plt.plot(energy2, prob, ls="-.", label="Theory", c="k", lw=2, alpha=0.5)

plt.legend()
plt.ylabel(r"$P_{\gamma -> a}$")
plt.semilogx()
plt.xlabel("Energy (keV)")
plt.semilogy()
plt.show()