#!/usr/bin/env python
import numpy as np 
import constants as c
import alpro 
import unittest

alpro.util.set_default_plot_params()
UNIT_LENGTH= 50676.79373667135
UNIT_GAUSS = 0.01953548032
HBAR_EV = 6.582119569e-16
c.MELEC = 9.10938356e-28
c.PARSEC = 3.0857e18
c.E = 4.8032045057134676e-10

#this our energy array in units of eV
energy = np.logspace(5,7,100000)
B = 1e-5            # 10 micro G
g = 1e-13 * 1e-9    # 1e-9 GeV^-1
mass = 1e-13        # 1e-12 eV
L = 10.0            # 10 kpc
ne = 100            # 1 particle cm^-3 (sets plasma frequency)


# function to get analytic prediction
def get_prediction(energy, mass, M, omega_pl, B, distance, norm=1):
    '''
    Get analytic prediction for idealised conversion probability.
    This function is actually written according to the notation
    of de Angelis et al. 2010 rather than the above equation -
    but is equivalent.
    '''
    SIM_UNITS = norm
    energy2 = energy / SIM_UNITS
    mass2 = mass / SIM_UNITS
    M2 = M / SIM_UNITS
    omega_pl2 = omega_pl / SIM_UNITS
    distance2 = distance * SIM_UNITS
    B2 = B / SIM_UNITS / SIM_UNITS

    # calculate deltas 
    Delta_PL = -(omega_pl2 * omega_pl2) / 2.0 / energy2
    Delta_AA = -(mass2 * mass2) / 2.0 / energy2
    Delta_AG = B2 / 2.0 / M2
    x = (Delta_PL - Delta_AA) ** 2
    Delta_OSC = np.sqrt (x + (4.0 * Delta_AG * Delta_AG))

    # put it all together 
    term1 = (B2 / M2 / Delta_OSC)**2
    term2 = np.sin(Delta_OSC * distance2 / 2.0)**2
    prediction = term1 * term2
    return (prediction, term1)


class RunTest(unittest.TestCase):

    def analytic_test(self):
        # Initial state {0,1,0}
        Ainit = np.zeros( (len(energy),6))
        Ainit[:,2] = 1.0
        phi = np.zeros_like(energy) #aligned field

        # get code result and plot 
        P, Anew = alpro.get_P(energy, Ainit, phi, B, L, g, mass, ne)

        # get theoretical prediction and plot
        omega_pl = np.sqrt(4.0 * np.pi * c.E * c.E * ne / c.MELEC) * HBAR_EV
        prediction, alpha = get_prediction (energy, mass, 1.0/g, omega_pl, B * UNIT_GAUSS, L * 1000.0 * c.PARSEC * UNIT_LENGTH)

        frac_error = (prediction - P) / P
        self.assertTrue (frac_error.all() < 1e-6)

    def discretise_test(self):
        # re-initialise 
        Ainit = np.zeros( (len(energy),6))
        Ainit[:,2] = 1.0
        phi = np.zeros_like(energy)

        # get single probability and plot
        P1, Anew = alpro.get_P(energy, Ainit, phi, B, L, g, mass, ne)

        # now do the discretised version
        Ainit = np.zeros( (len(energy),6))
        Ainit[:,2] = 1.0
        Anew = Ainit
        phi = np.zeros_like(energy)
        Nshells = 10
        for i in range(Nshells):
            Ainit = Anew
            P2, Anew = alpro.get_P(energy, Ainit, phi, B, L/Nshells, g, mass, ne)
        frac_error = (P2 - P1) / P1
        self.assertTrue (frac_error.all() < 1e-6)

    def marsh_test(self):

        # load data from Marsh and Libanov codes 
        import os
        folder = os.path.dirname(alpro.__file__)
        energies_m, P_m = np.genfromtxt("{}/data/marsh_test.dat".format(folder), unpack=True)

        energies = np.logspace(8,11,1000)
        s = alpro.Survival("libanov")
        s.init_model()
        mass = 1e-9
        g = 1e-11 * 1e-9
        s.set_params(g, mass)
        s.domain.rm = 0.0
        P = 1.0 - s.get_curve(energies_m, 1.0, 93.0, r0=0.0)

        frac_error = (P - P_m) / P
        self.assertTrue (frac_error.all() < 1e-6)

        

       







# import matplotlib.pyplot as plt 
# import numpy as np 
# import constants as c
# import alpro 
# alpro.util.set_default_plot_params()
# # now we're going to test against Jamie Davies' code.
# s1 = s = alpro.Survival("davies")
# s1.domain = alpro.models.FieldModel(None)
# mass = 1e-8
# g = 1e-11 * 1e-9
# s1.set_params(g, mass)

# r, B, phi, ne = np.genfromtxt("ClusterField.txt", unpack=True)
# last_index =-10
# s1.domain.ne = ne
# s1.domain.phi = phi
# r = np.append(r, 500)
# # s1.domain.r=r[:-1]
# s1.domain.r = 0.5 * (r[1:] + r[:-1])
# s1.domain.B = B * 1e-6
# #s1.domain.deltaL = 0.111111111111 * np.ones_like(s1.domain.r)
# #r = np.linspace(0,500,num=4500)
# s1.domain.deltaL = r[1:] - r[:-1] 

# energies, P_jamie = np.genfromtxt('noqedorcmb_james_test_m10_g1e11_pggs.txt', unpack=True)
# energies2 = np.logspace(9,14,1000)
# P, Prad = s1.propagate(s1.domain, energies2)

# plt.plot(energies, P_jamie, label="JD code")
# plt.plot(energies2, 1.0 - P, label="JM code", ls="--")
# plt.semilogx()
# plt.semilogy()
# plt.xlim(1e9,1e10)
# plt.legend()
# plt.xlabel("$E$ (eV)", fontsize=16)
# plt.ylabel("$P_{\gamma->\gamma}$", fontsize=16)
# plt.savefig("test.png", dpi=200)


# # In[8]:


# energies, P_jamie = np.genfromtxt('james_field_mm13_gm12p3_pggs.txt', unpack=True)
# plt.plot(energies, P_jamie, label="JD code", ls="-")

# energies = np.logspace(3,4,10000)
# # now we're going to test against Jamie Davies' code.
# s1 = s = alpro.Survival("davies")
# s1.domain = alpro.models.FieldModel(None)

# logg = 12.3
# logm = 13
# mass = (10.0**-logm)
# g = (10.0**-logg) * 1e-9
# s1.set_params(g, mass)

# r1, r2, Bx, By, ne = np.genfromtxt("domain_seed_100_1275b.dat", unpack=True)
# s1.domain.ne = ne 
# s1.domain.phi = np.arctan(By/Bx) 
# B = np.sqrt(Bx * Bx + By * By) 
# s1.domain.r = r2
# s1.domain.B = B 
# s1.domain.deltaL = (r2 - r1) 

# P, Prad = s1.propagate(s1.domain, energies)
# p_to_plot = 1.0 - P

# plt.plot(energies, p_to_plot, label="JM code", ls="--")
# plt.semilogx()
# plt.semilogy()
# plt.xlim(1e3,1e4)
# plt.xlabel("$E$ (eV)", fontsize=16)
# plt.ylabel("$P_{\gamma->\gamma}$", fontsize=16)
# plt.savefig("test_1275.png", dpi=200)


# # In[ ]:

if __name__ == '__main__':
    unittest.main()




