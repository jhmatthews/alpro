import alpro
import numpy as np
import matplotlib.pyplot as plt
import time 

def mod_squared(vector, i):
	'''
	get mod squared, ith element, for 3 complex vector
	''' 
	z = np.complex(vector[i], vector[i+3])
	return abs(z)**2





plt.clf()
energies = np.logspace(3,4,1000)

N = 100
Ls = np.linspace(1,10,3) 
gs = np.logspace(-14,-10,N)
ms = np.logspace(-14,-10,N)
gs = [1e-12]
ms = [1e-12]
ne = 1e-2
Ainit = np.array([0.0, 1.0, 0.0, 0.0, 0.0, 0.0])

t1 = time.time()
i = 0
for L in Ls:
	for g in gs:
		for m in ms:
			P = alpro.get_P(energies, Ainit, 0.0, 1e-6, L, g * 1e-9, m, ne)
			# print (i)
			plt.plot(energies, P)
			i += 1

print ( (time.time()-t1) / N**3)
plt.show()
