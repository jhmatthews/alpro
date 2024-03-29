#!/usr/bin/env python
import numpy as np
import alpro
import unittest
import os
import matplotlib.pyplot as plt
from datetime import date, datetime

# this our energy array in units of eV
energy_gamma = np.logspace(5, 7, 1000)
energy_x = np.logspace(3, 4, 1000)
B = 1e-5            # 10 micro G
g = 1e-12 * 1e-9    # 1e-9 GeV^-1
mass = 1e-11        # 1e-12 eV
L = 10.0            # 10 kpc
ne = 0.01            # 1 particle cm^-3 (sets plasma frequency)


class RunTest(unittest.TestCase):

    def make_plot(self, energy, P1, P2, subplot=111, logy=False, title=None):
        ax = plt.subplot(subplot)
        ax.plot(energy, P2, ls="-", lw=2.5)
        ax.plot(energy, P1, ls="--", lw=2.5)
        ax.set_xscale("log")
        if logy:
            ax.set_yscale("log")
        plt.title(title)

    def test_analytic(self, subplot=321):

        # set up uniform modek
        s = alpro.Survival("uniform")
        s.setup_regular_model(B=B, L=L, ne=ne, N=1, phi=0.0)

        s.set_params(g=g, mass=mass)   # set up ALP parameters in eV^-1 and eV
        # propagate the photon-ALP beam
        P, _ = s.propagate(energies=energy_x, pol="y")

        # get theoretical prediction
        prediction = s.analytic_pga(energy_x)

        # compare results
        self.assertTrue(np.allclose(P, prediction))
        self.make_plot(
            energy_x,
            prediction,
            P,
            subplot=subplot,
            logy=True,
            title="Analytic")

    def test_discretise(self, subplot=322):

        # set up uniform modek
        s1 = alpro.Survival("uniform")
        s1.setup_regular_model(B=B, L=L, ne=ne, N=10, phi=0.0)

        s2 = alpro.Survival("uniform")
        s2.setup_regular_model(B=B, L=L, ne=ne, N=1, phi=0.0)

        # set up ALP parameters in eV^-1 and eV
        s1.set_params(g=g, mass=mass)
        s2.set_params(g=g, mass=mass)

        # propagate the photon-ALP beam
        P1, _ = s1.propagate(energies=energy_x, pol="y")
        # propagate the photon-ALP beam
        P2, _ = s2.propagate(energies=energy_x, pol="y")

        # compare results
        self.assertTrue(np.allclose(P1, P2))
        self.make_plot(
            energy_x,
            P1,
            P2,
            subplot=subplot,
            logy=True,
            title="Discretisation")

    def test_marsh(self, subplot=323):

        # load data from Marsh code
        folder = os.path.dirname(alpro.__file__)
        energies_m, P_m = np.genfromtxt(
            "{}/data/marsh_test.dat".format(folder), unpack=True)
        P_m = 1.0 - P_m

        s = alpro.Survival("libanov")
        s.init_model()
        mass = 1e-9
        g = 1e-11 * 1e-9
        s.set_params(g, mass)
        s.domain.rm = 0.0

        # Marsh energies are in GeV
        P = s.get_curve(energies_m * 1e9, 1.0, 93.0, r0=0.0)
        self.assertTrue(np.allclose(1 - P, 1 - P_m, rtol=4e-3))
        self.make_plot(
            energies_m,
            P_m,
            P,
            subplot=subplot,
            title="Libanov Field")

    def test_fourier_single(self, subplot=324):

        # set up uniform modek
        s1 = alpro.Survival("uniform")
        s1.setup_regular_model(B=B, L=L, ne=0.0, N=1, phi=0.0)
        s2 = alpro.Survival("uniform")
        s2.setup_regular_model(B=B, L=L, ne=0.0, N=1, phi=0.0)

        # set up ALP parameters in eV^-1 and eV
        s1.set_params(g=g, mass=mass)
        s2.set_params(g=g, mass=mass)

        E, P1 = s1.propagate_fourier(
            s1.domain, pol="y", mode="massive", N=10000, f_res=200)
        prediction = s2.analytic_pga(E)
        #P1, _ = s1.propagate(s1.domain, E, pol="y")

        #self.assertTrue (np.allclose(1-P1[E>1e4], prediction[E>1e4], rtol=2))
        self.make_plot(
            E,
            prediction,
            P1,
            subplot=subplot,
            title="Fourier Single",
            logy=True)

    def test_fourier_massive(self, subplot=325):

        # set up cell model with uniform cell sizes
        s1 = alpro.Survival("1275b")
        s1.init_model()
        s1.set_params(1e-14 * 1e-9, 1e-9)
        s1.set_coherence_r0(None)
        s1.domain.create_box_array(100.0, 0, 10.0, r0=0)
        #print (s1.domain.deltaL)

        s2 = alpro.Survival("1275b")
        s2.init_model()
        s2.set_params(1e-14 * 1e-9, 1e-9)
        s2.set_coherence_r0(None)
        s2.domain.create_box_array(100.0, 0, 10.0, r0=0)
        #print (s2.domain.deltaL)

        E, P2 = s2.propagate_fourier(
            s2.domain, pol="both", mode="auto", N=100000, f_res=2000)
        P1, _ = s1.propagate(s1.domain, E, pol="both")

        select = (P1 > 1e-7)
        self.assertTrue(np.allclose(P1[select], P2[select], rtol=1e-2))
        self.make_plot(
            E,
            P1,
            P2,
            subplot=subplot,
            title="Fourier Massive",
            logy=True)
        # plt.ylim(1e-9,2)

    def test_fourier_massless(self, subplot=326):

        # set up cell model with uniform cell sizes
        s1 = alpro.Survival("1275b")
        s1.init_model()
        s1.set_params(1e-12 * 1e-9, 1e-13)
        s1.set_coherence_r0(None)
        s1.domain.create_box_array(500.0, 0, 10.0, r0=0)

        s2 = alpro.Survival("1275b")
        s2.init_model()
        s2.set_params(1e-12 * 1e-9, 1e-13)
        s2.set_coherence_r0(None)
        s2.domain.create_box_array(500.0, 0, 10.0, r0=0)

        N = 100000
        f_res = N // 40
        E, P2 = s2.propagate_fourier(
            s2.domain, pol="both", mode="massless", N=N, f_res=f_res)

        select = (E > 1e2) * (E < 1e4)
        P1, _ = s1.propagate(s1.domain, E[select], pol="both")

        #self.assertTrue (np.allclose(P1, P2[select], rtol=0.1))
        self.make_plot(
            E[select],
            P1,
            P2[select],
            subplot=subplot,
            title="Fourier Massless",
            logy=True)

    def title_string(self):
        t = """ALPro version {}
        Basic test on {}
        """.format(alpro.__version__, datetime.now().strftime("%d/%m/%Y at %H:%M:%S"))
        return (t)

    def test_save(self):
        t = self.title_string()
        plt.suptitle(t)
        plt.gcf().set_size_inches(8, 11)
        plt.savefig("Test_{}.png".format(date.today().strftime("%d-%m-%Y")))


if __name__ == '__main__':
    print("Testing ALPRO")
    print("Location: {}".format(alpro.__file__))
    unittest.main()
