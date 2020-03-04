import numpy as np 
import matplotlib.pyplot as plt 

class my_powerlaw:
	'''
	container for generating random numbers from a powerlaw with slope < 1
	'''
	def __init__(self, n=1.2, xmin=3.5, xmax=10.0):
		self.n = n
		self.alpha = 1.0 - self.n
		self.xmin = xmin 
		self.xmax = xmax 

	def rvs(self, size=None):
		'''
		generate random variables
		'''
		if size == None:
			z = np.random.random()
		else:
			z = np.random.random(size=size)

		term1 = z * (self.xmax ** self.alpha)
		term2 = (1. - z) * (self.xmin ** self.alpha)
		x = (term1 + term2) ** (1.0 / self.alpha)
		return (x)

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




