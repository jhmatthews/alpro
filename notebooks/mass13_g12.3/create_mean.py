import numpy as np 
import matplotlib.pyplot as plt 

Ptot = np.zeros( (100,7997) )
for seed in range(1,101):
	file = "PE_seed_{}_ma_13_g_12.3_mod.dat".format(seed)
	E1, E2, P = np.genfromtxt(file, unpack=True)
	Ptot[seed-1] = P
	E = 0.5 * (E1 + E2)
	plt.plot(E, P, alpha=0.3, lw=1, c="k")

Pmean = np.mean(Ptot, axis=0)
Pmedian = np.median(Ptot, axis=0)
Ps = np.std(Ptot, axis=0)
arr_to_save = np.column_stack( (E1, E2, Pmean, Ps) )
fname = "PE_mean_ma_13_g_12.3_mod.dat"

np.savetxt(fname, arr_to_save)

plt.plot(E, Pmean)
plt.plot(E, Pmedian)
plt.semilogx()
plt.semilogy()
plt.xlim(1,10)
plt.ylim(0.8,1.0)
plt.show()
