import matplotlib.pyplot as plt 
import numpy as np 
import alpro 
alpro.util.set_default_plot_params(tex=True)


def custom_energy_grid(obs_frame_cut = 1e4 , z=0.299, 
	                   dE_fine = 1.0, dE_coarse = 10.0):
	# 10 keV cutoff in observer frame 
	rest_frame_cut = obs_frame_cut * (1.0 + z)

	range_fine = rest_frame_cut - 1.0
	nbins_fine = range_fine / dE_fine
	E1 = np.linspace(1e3, rest_frame_cut,nbins_fine)

	range_coarse = 1e5 - rest_frame_cut
	nbins_coarse = range_coarse / dE_coarse
	E2 = np.linspace(rest_frame_cut, 1e5, nbins_coarse)
	Eall = np.concatenate( (E1[:-1], E2))
	return (Eall)



# initialise model 
s = alpro.Survival("1821")
s.init_model()

# set axion parameters
logm = 13
logg = 12.3
g = (10.0**(-1.0 * logg)) * 1e-9   # GeV^-1
mass = 10.0 ** (-1.0 * logm)      # 1e-12 eV

s.set_params(g, mass)
delta = 9.0 / 1000.0


Emin = np.log10(0.3)
Emax = np.log10(13)
Eall = custom_energy_grid(obs_frame_cut = 1e4 , z=0.299, 
	                      dE_fine = 1.0, dE_coarse = 10.0)

print (len(Eall))
E1 = Eall[:-1]
E2 = Eall[1:]

#E1 = np.arange(1,10,)
#E1, E2, Ptest = np.genfromtxt("PE_seed_1_ma_13_g_12.3_mod.dat", unpack=True)
energies = 0.5 * (E1 + E2) * 1e3


# delta_e = E2 - E1 
# plt.plot(delta_e)
# plt.semilogy()
# plt.show()


N = 5
plt.figure(figsize=(8,6))
import time
for seed in np.arange(N):
	t1 = time.time()
	P = 1.0 - s.get_curve(energies, seed, 583.0, r0=1)
	print (time.time() - t1)

	plt.plot(energies, P, label="Seed {}".format(seed))

	arr_to_save = np.column_stack( (E1, E2, P) )
	fname = "alpro_PE_seed_{}_ma_{}_g_{}_1821.dat".format(seed, logm, logg)
	np.savetxt(fname, arr_to_save)

plt.semilogx()
plt.xlabel("$E$ (keV)", fontsize=16)
plt.ylabel("$P_{\gamma->\gamma}$", fontsize=16)
plt.legend()
plt.savefig("alpro_PE_seeds_ma_{}_g_{}_1821.png".format(seed, logm, logg))


# r1 = s1.domain.r
# r2 = s1.domain.r + s1.domain.deltaL
# Bx = s1.domain.Bx
# By = s1.domain.By



