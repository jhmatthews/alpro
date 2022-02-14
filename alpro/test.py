#!/usr/bin/env python
import numpy as np 
import alpro 
import unittest
import os
import matplotlib.pyplot as plt 
from datetime import date

#this our energy array in units of eV
energy_gamma = np.logspace(5,7,1000)
energy_x = np.logspace(3,4,1000)
B = 1e-5            # 10 micro G
g = 1e-12 * 1e-9    # 1e-9 GeV^-1
mass = 1e-11        # 1e-12 eV
L = 10.0            # 10 kpc
ne = 0.01            # 1 particle cm^-3 (sets plasma frequency)



class RunTest(unittest.TestCase):

    def make_plot(self, energy, P1, P2, subplot=221, logy=False):
        plt.subplot(subplot)
        plt.plot(energy, P2, ls="-")
        plt.plot(energy, P1, ls="--")
        plt.semilogx()
        if logy:
            plt.semilogy()

    def test_analytic(self):

        # set up uniform modek 
        s = alpro.Survival("uniform")
        s.setup_regular_model(B = B, L = L, ne = ne, N = 1, phi = 0.0)


        s.set_params(g = g, mass = mass)   # set up ALP parameters in eV^-1 and eV 
        P, _ = s.propagate(energies=energy_x, pol="y")        # propagate the photon-ALP beam

        # get theoretical prediction
        prediction = s.analytic_pga(energy_x)

        # compare results
        self.assertTrue (np.allclose(P, prediction))
        self.make_plot(energy_x, prediction, P, subplot=221, logy=True)

    def test_discretise(self):

        # set up uniform modek 
        s1 = alpro.Survival("uniform")
        s1.setup_regular_model(B = B, L = L, ne = ne, N = 10, phi = 0.0)

        s2 = alpro.Survival("uniform")
        s2.setup_regular_model(B = B, L = L, ne = ne, N = 1, phi = 0.0)

        # set up ALP parameters in eV^-1 and eV
        s1.set_params(g = g, mass = mass)    
        s2.set_params(g = g, mass = mass) 

        P1, _ = s1.propagate(energies=energy_x, pol="y")        # propagate the photon-ALP beam
        P2, _ = s2.propagate(energies=energy_x, pol="y")        # propagate the photon-ALP beam

        # compare results
        self.assertTrue (np.allclose(P1, P2))
        self.make_plot(energy_x, P1, P2, subplot=222, logy=True)

    def test_marsh(self):

        # load data from Marsh code
        folder = os.path.dirname(alpro.__file__)
        energies_m, P_m = np.genfromtxt("{}/data/marsh_test.dat".format(folder), unpack=True)

        s = alpro.Survival("libanov")
        s.init_model()
        mass = 1e-9
        g = 1e-11 * 1e-9
        s.set_params(g, mass)
        s.domain.rm = 0.0

        # Marsh energies are in GeV
        P = 1.0 - s.get_curve(energies_m*1e9, 1.0, 93.0, r0=0.0)
        self.assertTrue (np.allclose(P, P_m, rtol=4e-3))
        self.make_plot(energies_m, P_m, P, subplot=223)

    def test_save(self):
        plt.savefig("Test_{}.png".format(date.today().strftime("%d-%m-%Y")))


if __name__ == '__main__':
    print ("Testing ALPRO")
    print ("Location: {}".format(alpro.__file__))
    unittest.main()



