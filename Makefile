CC = gcc

# point this to where GSL is installed 
GSL_LOC = /Users/matthews/winds/python/
LIB = $(GSL_LOC)lib
INCLUDE = $(GSL_LOC)include
#BIN = ../bin

ifeq (D,$(firstword $(MAKECMDGOALS)))
# use pg when you want to use gprof the profiler
# to use profiler make with arguments "make D python"
# this can be altered to whatever is best
	CFLAGS = -g -pg -Wall $(EXTRA_FLAGS) -I$(INCLUDE) $(MPI_FLAG)
	FFLAGS = -g -pg
	PRINT_VAR = DEBUGGING, -g -pg -Wall flags
else
# Use this for large runs
	CFLAGS = -O3 -Wall $(EXTRA_FLAGS) -I$(INCLUDE) $(MPI_FLAG)
	FFLAGS =
	PRINT_VAR = LARGE RUNS, -03 -Wall flags
endif

LDFLAGS+= -L$(LIB) -lm -lgsl -lgslcblas

startup:
	@echo $(COMPILER_PRINT_STRING)			# prints out compiler information
	@echo 'YOU ARE COMPILING FOR' $(PRINT_VAR)	# tells user if compiling for optimized or debug

prototypes:
	cproto -I$(INCLUDE) $(alpro_source) > prototypes.h

alpro_objects = alpro.o
alpro_source = alpro.c

alpro: startup $(alpro_objects)
	$(CC) ${CFLAGS} $(alpro_objects) $(LDFLAGS) -o alpro
	#cp $@ alpro


all: clean alpro

clean :
	rm -f *.o  *~

# lslt /Library/Frameworks/Python.framework/Versions/3.7/include/python3.7m/
# gcc -fPIC -shared -I/Users/matthews/winds/python/include -I/Library/Frameworks/Python.framework/Versions/3.7/include/python3.7m/ -L/Users/matthews/winds/python/lib -lm -lgsl -lgslcblas -lpython3.7 -o alpro.so alpro.c

