import numpy as np 
import matplotlib.pyplot as plt 
from alpro.models import unit

class my_powerlaw:
	'''
	container for generating random numbers from a powerlaw with slope < 1. 
	For a power of form x**-n, with b>1, from xmin to xmax, the resulting
	CDF can be shown to take the analytic form

		CDF(x) = (xmin^alpha - x^alpha) / (xmin^alpha - xmax^alpha)

	where alpha = 1-n. Thus, if z is a random number in range (0,1),
	then a random variable rv can be found by inverting this expression,
	i.e. rv = [z * (b^alpha - a^alpha) + a^alpha] ^ (1/alpha).

	This is embedded in a class so an individual function can be passed without
	additional arguments, the function rvs() generates the actual random numbers
	similarly to scipy.stats distribution syntax.
	'''
	def __init__(self, n=1.2, xmin=3.5, xmax=10.0):
		'''
		initialise the powerlaw parameters
		'''
		self.n = n
		self.alpha = 1.0 - self.n
		self.xmin = xmin 
		self.xmax = xmax 

	def rvs(self, size=None):
		'''
		generate (size) random variables. if size is None, 
		generate a single float.
		'''
		if size == None:
			z = np.random.random()
		else:
			z = np.random.random(size=size)

		term1 = z * (self.xmax ** self.alpha)
		term2 = (1. - z) * (self.xmin ** self.alpha)
		rv = (term1 + term2) ** (1.0 / self.alpha)
		return (rv)

class modelb_powerlaw:
	'''
	container for generating random numbers from a powerlaw with slope < 1. 
	For a power of form x**-n, with b>1, from xmin to xmax, the resulting
	CDF can be shown to take the analytic form

		:math:`CDF(x) = (xmin^alpha - x^alpha) / (xmin^alpha - xmax^alpha)`

	where alpha = 1-n. Thus, if z is a random number in range (0,1),
	then a random variable rv can be found by inverting this expression,
	i.e. rv = [z * (b^alpha - a^alpha) + a^alpha] ^ (1/alpha).

	This is embedded in a class so an individual function can be passed without
	additional arguments, the function rvs() generates the actual random numbers
	similarly to scipy.stats distribution syntax.
	'''
	def __init__(self, n=1.2, xmin=3.5, xmax=10.0):
		'''
		initialise the powerlaw parameters
		'''
		self.n = n
		self.alpha = 1.0 - self.n
		self.xmin = xmin 
		self.xmax = xmax 

	def rvs(self, r, size=None):
		'''
		generate (size) random variables. if size is None, 
		generate a single float.
		'''
		if size == None:
			z = np.random.random()
		else:
			z = np.random.random(size=size)

		term1 = z * (self.xmax ** self.alpha)
		term2 = (1. - z) * (self.xmin ** self.alpha)
		rv = (term1 + term2) ** (1.0 / self.alpha)
		return (rv)

def delta_eff(energy, ne, mass):
	r'''
	Calculate delta_eff for a B field in Gauss, a mass in eV and a photon energy in eV
	:math:`\Delta_{\rm eff} = \frac{m_A^2 - \omega_p^2}{4 \omega}`
	'''
	omega_p = np.sqrt (4.0 * np.pi * unit.e * unit.e * ne / unit.melec) * unit.hbar_ev;
	delta = ((mass * mass) - (omega_p * omega_p)) / 4.0 / energy
	return (delta)

def theta(energy, B, ne, mass, g_a):
	r'''
	Calculate capital Theta for a B field in Gauss, a mass in eV and a photon energy in eV
	:math:`\Theta = \frac{2 B_\perp \omega g_a}{m_A^2 - \omega_p^2}`.
	'''
	B *= unit.unit_gauss_natural
	omega_p = np.sqrt (4.0 * np.pi * unit.e * unit.e * ne /
	  unit.melec) * unit.hbar_ev;
	theta = 2.0 * B * energy * g_a 
	m_effsq = (mass * mass) - (omega_p * omega_p)
	#print (mass, omega_p)
	theta /= m_effsq
	return (theta)

def Losc(energy, B, ne, mass, g_a):
	'''
	calculate the ALP oscillation length in cm
	'''
	delta = delta_eff(energy, ne, mass)
	th = theta(energy, B, ne, mass, g_a)
	L = 2.0 * np.pi / (delta * np.sqrt(1.0 + th**2)) / unit.unit_length_natural / unit.kpc
	return (L)

def P_gamma_analytic(energy, B, ne, mass, g_a, L):
	'''
	L should be in kpc!
	ne, B in CGS
	everything else in eV or 1/eV
	'''
	delta = delta_eff(energy, ne, mass)
	th = theta(energy, B, ne, mass, g_a)
	phase = delta * np.sqrt(1.0 + th**2) * L * unit.kpc * unit.unit_length_natural
	amplitude = th**2 / (1.0 + th**2)
	sin_term = np.sin(phase)
	return (amplitude * sin_term * sin_term)



def set_default_plot_params(tex=False, dpi=300):
	if tex:
		plt.rcParams["text.usetex"] = "True"
		plt.rcParams['font.family']='serif'

	plt.rcParams['axes.labelsize'] = 16
	plt.rcParams['figure.figsize'] = [6.4, 4.8]
	plt.rcParams['savefig.dpi'] = dpi
	plt.rcParams['xtick.labelsize']=14
	plt.rcParams['ytick.labelsize']=14
	plt.rcParams['font.serif']=['cm']
	#plt.rcParams['font.family']='serif'	
	plt.rcParams["lines.linewidth"] = 3
	plt.rcParams["axes.linewidth"] = 1.5
	plt.rcParams["xtick.major.width"] = 1.5
	plt.rcParams["ytick.major.width"] = 1.5




