from alpro import models, util
import numpy as np 
import matplotlib.pyplot as plt 
util.set_default_plot_params(tex=True)
r = np.logspace(0,3,num=1000)

russell = models.ClusterProfile(model="russell")
modA = models.ClusterProfile(model="a")
modB = models.ClusterProfile(model="b")

plt.plot(r, russell.get_B(r)/1e-6, label="Model for H1821$+$643", lw=3, c="k")
plt.plot(r, modA.get_B(r)/1e-6, label="Reynolds+ 2020 Model A, Perseus", lw=3, c="k", ls="--")
#plt.plot(r, modB.get_B(r), label="Model B, Perseus")
plt.loglog()
plt.ylabel(r"$B (\mu\mathrm{G})$", fontsize=16)
plt.xlabel("$r$ (kpc)", fontsize=16)
plt.ylim(5,55)
#plt.xlim(1,500)
yticks = [5,10,20,30,40,50]
plt.yticks(yticks)
plt.gca().set_yticklabels([str(x) for x in yticks])

xticks = [1,10,100]
plt.xticks(xticks)
plt.gca().set_xticklabels([str(x) for x in xticks])
# plt.ylim(5,50)
plt.xlim(1,500)
plt.legend(fontsize=12)
plt.subplots_adjust(top=0.98, bottom=0.12, right=0.95)
plt.savefig("cluster-profiles.png", dpi=300)