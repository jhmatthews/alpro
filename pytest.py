
from ctypes import *
so_file = "alpro.so"
my_functions = CDLL(so_file)
my_functions.main()