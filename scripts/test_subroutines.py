import matplotlib.pyplot as plt 
import numpy as np 
from constants import *
import constants as c
import alpro

UNIT_TIME  = 1519252268630104.8
UNIT_LENGTH= 50676.79373667135
UNIT_MASS  = 5.6095363761802584e+32
UNIT_GAUSS = 0.06925147467360344
HBAR_EV = 6.582119569e-16

def get_density(r):
	x1 = 3.9e-2 / ((1.0 + (r/80.0))**1.8)
	x2 = 4.05e-3 / ((1.0 + (r/280.0))**0.87)
	return (x1 + x2)

def get_B(r):
	x = r/(93.0)
	Bx = (0.00312443*(x**18)) - (0.0319991*(x**16)) + (0.260311*(x**14)) - (1.63197*(x**12)) + (7.58002*(x**10)) - (24.721*(x**8)) + (52.3929*(x**6)) - (63.8794*(x**4)) + (35.8973*(x**2)) - 5.86899	
	By = (0.0102459*(x**17))-(0.0937683*(x**15)) + (0.671841*(x**13)) - (3.6406*(x**11)) + (14.2479*(x**9)) - (37.7455*(x**7)) + (61.3611*(x**5)) - (51.7231*(x**3)) + (16.9128*x)
	return Bx, By


def run_marsh_test(energy, g_a, mass):
	Ainit = np.zeros( (len(energy),6))
	Ainit2 = np.zeros( (len(energy),6))
	Ainit[:,2] = 1.0
	Ainit2[:,0] = 1.0
	myL = 0.0
	L = 1
	EPSILON = 1e-10
	r = -L / 2.0
	Lmax = 93.0

	while myL < (Lmax-EPSILON):
		myL += L 
		r += L 

		Bx,By = get_B(r)
		B = np.sqrt(Bx**2 + By**2)

		phi = (np.arctan(Bx/By) * np.ones_like(energy)) 
		phi2 = (np.arctan(Bx/By) * np.ones_like(energy)) 
		
		B *= 1e-6 
		P1, Anew = alpro.get_P(energy, Ainit, phi, B, L, g_a, mass, 1e-20)
		P2, Anew2 = alpro.get_P(energy, Ainit2, phi2, B, L, g_a, mass, 1e-20)
		Ainit = Anew
		Ainit2 = Anew2

	P = 0.5 * (P1 + P2)
	return (P)