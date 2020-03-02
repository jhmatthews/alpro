import alpro
import numpy as np
import matplotlib.pyplot as plt

plt.clf()
energies = np.logspace(3,4,1000)

Ls = [1.0,5.0,10.0]
for L in Ls:
	plt.plot(energies, alpro.get_P(energies, 0.0, 1e-3, L))

plt.show()
