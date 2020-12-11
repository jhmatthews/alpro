import numpy as np 
import matplotlib.pyplot as plt 

m = 13 
g = 12.4
alpro_PE_ma_13.0_g_12.4.dat



for seed in range(0,500):
	file = "Seed{:03d}/alpro_PE_ma_{:.1f}_g_{:.1f}.dat".format(seed, m, g)
	E1, E2, P = np.genfromtxt(file, unpack=True)
	if seed == 0:
		Ptot = np.zeros((500,len(P)))
	Ptot[seed,:] = P
	#E = 0.5 * (E1 + E2)
	#plt.plot(E, P, alpha=0.3, lw=1, c="k")

Pmean = np.mean(Ptot, axis=0)
Pmedian = np.median(Ptot, axis=0)
Ps = np.std(Ptot, axis=0)
arr_to_save = np.column_stack( (E1, E2, Pmean, Ps) )
fname = "meanPE_ma_{:.1f}_g_{:.1f}.dat".format(seed, m, g)
np.savetxt(fname, arr_to_save)
