#import matplotlib.pyplot as plt 
import numpy as np 
import alpro 
#import matplotlib
import sys, os
from mpi4py import MPI
import time 

global FOLDER
FOLDER = "."

def custom_energy_grid(obs_frame_cut = 1e4 , z=0.299, 
	                   dE_fine = 1.0, dE_coarse = 10.0):
	# 10 keV cutoff in observer frame 
	rest_frame_cut = obs_frame_cut * (1.0 + z)

	range_fine = rest_frame_cut - 1.0
	nbins_fine = range_fine / dE_fine
	E1 = np.linspace(1e3, rest_frame_cut,nbins_fine)

	range_coarse = 1e5 - rest_frame_cut
	nbins_coarse = range_coarse / dE_coarse
	E2 = np.linspace(rest_frame_cut, 1e5, nbins_coarse)
	Eall = np.concatenate( (E1[:-1], E2))
	return (Eall)

def calculate_individual_curve(energies, i, m_a, g_a, dx=1.75, Lmax=1800.0, field_type="file"):
	
	# initialise class
	s1 = alpro.Survival(field_type)

	if field_type == "file":
		# read in B field model
		filename = "{}/saveBfield_{:03d}.npy".format(FOLDER, i)
		Bfield_model = alpro.models.ClusterFromFile(fname=filename, model_type="1d")

		# set up domain and g m pair
		s1.domain = alpro.models.FieldModel(None)
		s1.set_params(g_a, m_a)
		s1.domain.domain_from_slice(Bfield_model, deltaL=dx, Lmax=Lmax)

		# do the propagation
		P, Pradial = s1.propagate(s1.domain, energies)
	else:
		s1.init_model()
		s1.set_params(g_a, m_a)
		P = s1.get_curve(energies, i, 583.0)

	# return survival probability
	return (1.0 - P)


def save_curve(Eall, pgg, j, m_a, g_a):
	'''
	save a survival probability curve to file
	'''
	from scipy.interpolate import interp1d 
	E1 = Eall[:-1] 
	E2 = Eall[1:]
	Ecen = 0.5 * (E1 + E2)
	fint = interp1d(Ecen, pgg, kind="slinear")

	Efine = np.linspace(1e3,1.5e4,14001)
	E1 = Efine[:-1]
	E2 = Efine[1:]
	Ecen = 0.5 * (E1 + E2)
	pgg_fine = fint(Ecen)
	arr_to_save = np.column_stack( (E1/1e3, E2/1e3, pgg_fine))
	logm = -1.0 * np.log10(m_a)
	logg = -1.0 * np.log10(g_a)
	fname = "Seed{:03d}/alpro_PE_ma_{:.1f}_g_{:.1f}_1821.dat".format(int(j), int(j), logm, logg)
	np.savetxt(fname, arr_to_save, fmt="%8.3f %8.3f     %8.4e")

def set_param_arrays():
	g1 = 10.0**np.arange(-13,-11.8,0.1)
	g1 = 10.0**np.arange(-11.9,-11.8,0.1)
	g2 = 10.0**np.arange(-11.8,-10.4,0.2)
	g_a = np.concatenate((g1,g2))

	m_a = 10.0**np.arange(-13.7,-10,0.1)
	return (m_a, g_a)



def run_it(N_fields = 500, Lmax = 1800.0):
	# let's do this in parallel 
	nproc = MPI.COMM_WORLD.Get_size()       # number of processes
	my_rank = MPI.COMM_WORLD.Get_rank()     # The number/rank of this process
	my_node = MPI.Get_processor_name()      # Node where this MPI process runs

	if my_rank == 0:
		for j in range(N_fields):
			os.system("mkdir -p Seed{:03d}".format(j))

	MPI.COMM_WORLD.Barrier()


	# specify g and m arrays 
	#g_a = 10.0**np.arange(-13,-10.6,0.1)
	#m_a = 10.0**np.arange(-13.7,-10,0.1)
	g_a = [10.0**-13]
	m_a = [10.0**-13.7]
	m_a, g_a = set_param_arrays()

	n_g = len(g_a)
	n_m = len(m_a)
	N_MODELS_TOTAL = n_m * n_g

	params = np.zeros((2, N_MODELS_TOTAL))
	itot = 0

	for i, g in enumerate(g_a):
		for j, m in enumerate(m_a):
			params[0,itot] = m
			params[1,itot] = g
			itot += 1

	# this MUST be an integer division so we get the remainder right 
	n_models = N_MODELS_TOTAL // nproc       # number of models for each thread
	remainder = N_MODELS_TOTAL - ( n_models * nproc )   # the remainder. e.g. your number of models may 

	# little trick to spread remainder out among threads. If say you had 19 total models, and 4 threads
	# then n_models = 4, and you have 3 remainder. This little loop would redistribute these three 
	if remainder < my_rank + 1:
	    my_extra = 0
	    extra_below = remainder
	else:
	    my_extra = 1
	    extra_below = my_rank

	# where to start and end your loops for each thread
	my_nmin = int((my_rank * n_models) + extra_below)
	my_nmax = int(my_nmin + n_models + my_extra)

	# total number you actually do
	ndo = my_nmax - my_nmin

	print ("This is thread {} calculating pairs {} to {}".format(my_rank, my_nmin, my_nmax))
	print ("total pairs {} models {}".format(N_MODELS_TOTAL, N_MODELS_TOTAL * N_fields))

	# set barrier so print output doesn't look muddled
	# just waits for other thread
	MPI.COMM_WORLD.Barrier()
	dx = Lmax / 2048.0

	# set up energies array in eV, high resolution for Athena
	#Eall = np.arange(1e3, 1.5e4 + 1.0, 1.0)
	#energies = np.arange(1e3, 1.4e4, 1.0)
	Eall = custom_energy_grid(obs_frame_cut = 1e4 , z=0.299, 
	                      dE_fine = 1.0, dE_coarse = 10.0)
	E1 = Eall[:-1] 
	E2 = Eall[1:]
	energies = 0.5 * (E1 + E2) 

	# start a timer for each thread
	time_init = time.time()
	tinit = time_init

	for i in range(my_nmin, my_nmax):
		m_a = params[0,i]
		g_a = params[1,i]
		print ("g m pair: {}/{}, {} {}".format(i, ndo, np.log10(m_a), np.log10(g_a)))
		for j in range(N_fields):
			print ("Seed:", j)
			pgg = calculate_individual_curve(energies, j, m_a, g_a, Lmax=Lmax, field_type="1821")
			print ("Curve Done")
			save_curve(Eall, pgg, j, m_a, g_a)
			print ("Curve saved")

	# get the time taken
	time2 = time.time() 
	time_tot = time2 - time_init

	# set barrier so print output doesn't look muddled
	if my_rank == 0: print ('Waiting for other threads to finish...')

	# another barrier, wait for all to finish
	MPI.COMM_WORLD.Barrier()

	print ("Thread {} took {} seconds to calculate {} models".format(my_rank, time_tot, ndo))

if __name__ == "__main__":
	print (alpro.__file__)
	#FOLDER = sys.argv[1]
	run_it(N_fields=5)

