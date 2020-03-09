import numpy as np 
import matplotlib.pyplot as plt 

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

def set_default_plot_params(tex=False, dpi=300):
	if tex:
		plt.rcParams["text.usetex"] = "True"
		plt.rcParams['font.family']='serif'

	plt.rcParams['savefig.dpi'] = dpi
	plt.rcParams['xtick.labelsize']=14
	plt.rcParams['ytick.labelsize']=14
	plt.rcParams['font.serif']=['cm']
	#plt.rcParams['font.family']='serif'	
	plt.rcParams["lines.linewidth"] = 3
	plt.rcParams["axes.linewidth"] = 1.5
	plt.rcParams["xtick.major.width"] = 1.5
	plt.rcParams["ytick.major.width"] = 1.5




