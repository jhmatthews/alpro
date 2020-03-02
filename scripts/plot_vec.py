import matplotlib.pyplot as plt 
import numpy as np 
from constants import *
import constants as c

def mod_squared(vector, i):
	'''
	get mod squared, ith element, for 3 complex vector
	''' 
	# if vector.shape == (6,):
	print (vector[0])
	z = np.complex(vector[i], vector[i+3])
	# else:
	# 	z = np.complex(vector[:,i], vector[:,i+3])
	return np.abs(z)**2


print (c.__file__)

energy, Ax, Ay, a = np.loadtxt("vectors.dat", usecols=(1,2,3,4), unpack=True)

# amp = (a+ia) * (a-ia)

#plt.plot(energy, a, label="Code", lw=3)

print (energy)

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
# energy2 = energy 
# UNIT_TIME = 1.0 / H * 2.0 * PI * EV2ERGS
# UNIT_LENGTH = 1.0 / H * 2.0 * PI / C * EV2ERGS
# UNIT_MASS = C * C / EV2ERGS
# UNIT_MAXWELL = (UNIT_LENGTH** 1.5) * (UNIT_MASS ** 0.5) / UNIT_TIME

# UNIT_GAUSS = ((C * H / 2.0 / PI)**1.5) / EV2ERGS / EV2ERGS

# UNIT_GAUSS2 = (UNIT_LENGTH**-0.5) * (UNIT_MASS ** 0.5) / UNIT_TIME
# print (UNIT_LENGTH, UNIT_GAUSS2, UNIT_GAUSS, UNIT_TIME, UNIT_MASS)

UNIT_TIME  = 1519252268630104.8
UNIT_LENGTH= 50676.79373667135
UNIT_MASS  = 5.6095363761802584e+32
UNIT_GAUSS = 0.06925147467360344
HBAR_EV = 6.582119569e-16

mass = 1e-12
M = 1.0 / (1e-12) * 1e9
g = 1e-12 * 1e-9
ne = 1e-3
omega_pl = np.sqrt(4.0 * PI * E * E * ne / MELEC) * HBAR_EV

print (omega_pl)
# / EV2ERGS
# omega_pl = 1e20
#omega_pl = 5.64e4 * np.sqrt(1e-2) * HEV / 2.0 / PI 
# B = (1e-5 ** 2) / 8.0 / PI * 1e7
B = 1e-6 * UNIT_GAUSS
phi = 0.0
L = 1.0
distance = L * 1000.0 * PARSEC * UNIT_LENGTH


# for logM in range(-20,20,1):
# 	M = 10.0 ** logM
# 	print (logM)
prediction, alpha = get_prediction2 (energy2, mass, M, omega_pl, B, distance)
# plt.plot(energy2, prediction)

prediction, alpha = get_prediction (energy2, mass, M, omega_pl, B, distance)
# plt.plot(energy2, prediction, c="r")

import alpro
Ainit = np.zeros( (len(energy2),6))
print (Ainit.shape)
Ainit[:,0] = 1.0

P, Anew = alpro.get_P(energy2, Ainit, 0, B / UNIT_GAUSS, L, g, mass, ne)
plt.plot(energy2, P, c="C0", lw=3, label="Code")

print (Anew.shape)

pp = [mod_squared(Anew[i], 2) for i in range(len(Anew))]
plt.plot(energy2, pp, c="r", lw=3)
# print (alpro.__doc__)

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