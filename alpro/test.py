#!/usr/bin/env python
import numpy as np 
import constants as c
import alpro 
import unittest
import os
import matplotlib.pyplot as plt 

#this our energy array in units of eV
energy = np.logspace(5,7,100000)
B = 1e-5            # 10 micro G
g = 1e-13 * 1e-9    # 1e-9 GeV^-1
mass = 1e-13        # 1e-12 eV
L = 10.0            # 10 kpc
ne = 100            # 1 particle cm^-3 (sets plasma frequency)

class RunTest(unittest.TestCase):

    def test_analytic(self):
        # Initial state {0,1,0}
        Ainit = np.zeros( (len(energy),6))
        Ainit[:,2] = 1.0
        phi = np.zeros_like(energy) #aligned field

        # get code result
        P, Anew = alpro.get_P(energy, Ainit, phi, B, L, g, mass, ne)

        # get theoretical prediction
        prediction = alpro.util.P_gamma_analytic(energy, B, ne, mass, g, L)

        # compare 
        frac_error = (prediction - P) / P
        self.assertTrue (np.all(frac_error < 1e-6))

    def test_discretise(self):
        # re-initialise 
        Ainit = np.zeros( (len(energy),6))
        Ainit[:,2] = 1.0
        phi = np.zeros_like(energy)

        # get single probability and plot
        P1, Anew = alpro.get_P(energy, Ainit, phi, B, L, g, mass, ne)

        # now do the discretised version
        Ainit = np.zeros( (len(energy),6))
        Ainit[:,2] = 1.0
        phi = np.zeros_like(energy)
        Nshells = 10
        for i in range(Nshells):
            P2, Anew = alpro.get_P(energy, Ainit, phi, B, L/Nshells, g, mass, ne)
            Ainit = Anew

        frac_error = np.fabs((P2 - P1) / P1)
        self.assertTrue (np.all(frac_error < 1e-6))

    def test_marsh(self):

        # load data from Marsh and Libanov codes 
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
        frac_error = np.fabs((P - P_m) / P)
        print ("Marsh test, max frac error:", np.max(frac_error))
        self.assertTrue (np.all(frac_error< 0.01))


if __name__ == '__main__':
    print ("Testing ALPRO")
    print ("Location: {}".format(alpro.__file__))
    unittest.main()




