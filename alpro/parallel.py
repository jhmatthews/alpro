#!/usr/bin/env python
import numpy as np 
import sys, os
from mpi4py import MPI
import time
import unittest
import alpro

def dummy(i):
    return i

def get_nmin_nmax_mpi(my_rank, nproc, Nmodels, Nstart=0, mode="chunks"):
    '''
    For a parallel task involving nproc processors and Nmodels models, work out which models
    the process my_rank is given. This splits up Nmodels between the processors equally, and
    then also distributes any remainder if Nmodels is not divisible by nproc.

    Parameters:
        my_rank     int
                    processor rank, number between 0 and nproc 

        nproc       int 
                    number of parallel processors

        Nmodels     int 
                    total number of models / calculations we want to do

        Nstart      int 
                    Starting point (in case you want to start higher than 0) 

        mode        str
                    how to split up the parallel processes.
    '''

    # this MUST be an integer division so we get the remainder right 
    n_models_per_thread = int(Nmodels // nproc)       # number of models for each thread
    remainder = Nmodels - ( n_models_per_thread * nproc )   # the remainder, since your number of models may not be divisible 

    if mode == "sequence":
        my_tasks = np.arange(my_rank,(n_models_per_thread * nproc), nproc, dtype=int)
        ndo = n_models_per_thread
        if remainder > my_rank:
            my_tasks = np.append(my_tasks, (n_models_per_thread * nproc) + my_rank)
            ndo += 1

    elif mode == "chunks":
        # little trick to spread remainder out among threads. If say you had 19 total models, and 4 threads
        # then n_models = 4, and you have 3 remainder. This little loop bit of code redistributes these 3
        if remainder < my_rank + 1:
            my_extra = 0
            extra_below = remainder
        else:
            my_extra = 1
            extra_below = my_rank

        # where to start and end your loops for each thread
        my_nmin = Nstart + int((my_rank * n_models_per_thread) + extra_below)
        my_nmax = Nstart + int(my_nmin + n_models_per_thread + my_extra)

        # total number you actually do
        ndo = my_nmax - my_nmin

        # tasks to iterate over 
        my_tasks = np.arange(my_nmin, my_nmax, 1, dtype=int)

    else:
        raise ValueError("mode must be chunks or sequence")

    assert (ndo == len(my_tasks))

    return (my_tasks, ndo)


def run(function, iterable, kwargs={}, split="chunks"):
    '''
    Run a function in parallel using an iterable as arguments.

    Parameters:
        function    callable
                    function to compute 

        iterable    iterable
                    iterable containing first arguments 

        kwargs      dictionary
                    dictionary of keyword arguments to pass to function 

        split       str 
                    how to split up the parallel processes

    Returns:
        time_tot    float 
                    time taken to compute 

    '''

    # let's get some of the information about processors running in parallel 
    nproc = MPI.COMM_WORLD.Get_size()       # number of processes
    my_rank = MPI.COMM_WORLD.Get_rank()     # The number/rank of this process
    Nmodels = len(iterable)

    # where to start and end your loops for each thread
    my_tasks, ndo = get_nmin_nmax_mpi(my_rank, nproc, Nmodels, Nstart=0, mode=split)


    print ("This is thread {} calculating models {} to {}, total {}".format(my_rank, my_tasks[0], my_tasks[-1], ndo))
    

    # set barrier so print output doesn't look muddled
    # just waits for other thread
    MPI.COMM_WORLD.Barrier()
    t1 = time.time() # set a timer

    for i in my_tasks:
        # call the function with the iterable as the argument 
        function(iterable[i], **kwargs)
        
    time_tot = time.time() - t1
    
    # set barrier so print output doesn't look muddled
    if my_rank == 0: print ('Waiting for other threads to finish...')

    # another barrier, wait for all to finish
    MPI.COMM_WORLD.Barrier()

    print ("Thread {} took {} seconds to calculate {} models".format(my_rank, time_tot, ndo))
    return (time_tot)

class RunTest(unittest.TestCase):
    def test_parallel(self):
        # function to use for parallelisation
        def function_for_para(i):
            s1 = alpro.Survival("1275b")
            s1.init_model()
            s1.set_params(1e-14 * 1e-9, 1e-13)
            s1.domain.create_box_array(1800.0, i, s1.coherence_func, r0=0)

        n = np.arange(0,100)
        run(function_for_para, n)

if __name__ == "__main__":
    print ("Testing ALPRO parallelisation routine")
    print ("Location: {}".format(alpro.__file__))
    unittest.main()

