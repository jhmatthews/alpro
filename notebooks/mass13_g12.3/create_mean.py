import numpy as np 

Ptot = np.zeros( (100,7997) )
for seed in range(1,101):
	file = "PE_seed_{}_ma_13_g_12.3_mod.dat".format(seed)
	E1, E2, P = np.genfromtxt(file, unpack=True)
	Ptot[seed-1] = P

Pmean = np.mean(Ptot, axis=0)
Ps = np.std(Ptot, axis=0)
arr_to_save = np.column_stack( (E1, E2, Pmean, Ps) )
fname = "PE_mean_ma_13_g_12.3_mod.dat"

np.savetxt(fname, arr_to_save)

