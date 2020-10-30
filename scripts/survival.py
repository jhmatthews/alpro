import matplotlib.pyplot as plt 
import numpy as np 
import alpro 
alpro.util.set_default_plot_params(tex=True)

def format_plot():
	plt.xlabel("$E$ (keV)", fontsize=16)
	plt.ylabel("$P_{\gamma->\gamma}$", fontsize=16)
	plt.subplots_adjust(top=0.97, bottom=0.12, right=0.95, left=0.1)
	ax1 = plt.gca()
	ax1.set_xscale('log')
	ax1.set_xlim(0.5,9.9)
	from matplotlib.ticker import ScalarFormatter, NullFormatter
	ax1.get_xaxis().set_minor_formatter(NullFormatter())
	xticks = [1,2,5]
	xtick_str = [str(x) for x in xticks]
	ax1.set_xticks(xticks)
	ax1.set_xticklabels(xtick_str)
	#ax1.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
	plt.ylim(0.8,1.0)
	plt.legend(fontsize=12)
	plt.yticks([0.8,0.9,1])

# initialise cluster models 
s1 = alpro.Survival("1275a")
s1.init_model()
s2 = alpro.Survival("1821")
s2.init_model()

# set axion parameters
logm = 13
logg = 12.3
g = (10.0**(-1.0 * logg)) * 1e-9   # GeV^-1
mass = 10.0 ** (-1.0 * logm)      # 1e-12 eV
s1.set_params(g, mass)
s2.set_params(g, mass)

# compute curves for different model instances
energies = np.logspace(2.5,5,1000)
N = 100
Ptot = np.zeros( (2,N,len(energies)) )
for seed in np.arange(N):
    Ptot[0,seed,:] = 1.0 - s1.get_curve(energies, seed, 500.0)
    Ptot[1,seed,:] = 1.0 - s2.get_curve(energies, seed, 500.0)

# make plots 
plt.clf()
labels = ["Perseus, Model A", "H1821$+$643 mean and std. dev."]
for i in range(2):
    P = Ptot[i,:,:]
    mean_P = np.mean(P, axis=0)
    plt.plot(energies/1e3, mean_P, label=labels[i], lw=4)
    sigma = np.std(P, axis=0)
    plt.fill_between(energies/1e3, mean_P-sigma, mean_P+sigma, alpha=0.5)


plt.plot(energies/1e3, Ptot[1,10,:], lw=2, c="k", ls=":", label="H1821$+$643, Realisation 1")
plt.plot(energies/1e3, Ptot[1,20,:], lw=2, c="k", ls="--", label="H1821$+$643, Realisation 2")
format_plot()
plt.savefig("survival-curves-rest.png", dpi=300)
plt.clf()

# now make observer frame plot 
redshifts = [0.017559, 0.299]
labels = ["Perseus, Model A, $z=0.018$", "H1821$+$643 mean and std. dev., $z=0.299$"]
# energies = energies / (1.0 + redshift)
for i in range(2):
	energies_obs = energies / (1.0 + redshifts[i])
	P = Ptot[i,:,:]
	mean_P = np.mean(P, axis=0)
	plt.plot(energies_obs/1e3, mean_P, label=labels[i], lw=4)
	sigma = np.std(P, axis=0)
	plt.fill_between(energies_obs/1e3, mean_P-sigma, mean_P+sigma, alpha=0.5)


plt.plot(energies_obs/1e3, Ptot[1,10,:], lw=2, c="k", ls=":", label="H1821$+$643, Realisation 1")
plt.plot(energies_obs/1e3, Ptot[1,20,:], lw=2, c="k", ls="--", label="H1821$+$643, Realisation 2")
format_plot()
plt.savefig("survival-curves-observed.png", dpi=300)



# plt.semilogy()