CC = g++
CCOMP = /Users/krishna/Programs/build/openmpi/bin/mpic++

OPTIMIZE = -O3

TWOFAST_INC = -I/Users/krishna/Research/github/twofast/include
TWOFAST_LIB = -L/Users/krishna/Research/github/twofast/lib -ltwofast

INC = $(TWOFAST_INC)
LIBS = -lm $(TWOFAST_LIB)
CFLAGS = $(OPTIMIZE)

all: TWOPOINT TWOPOINT_MPI

TWOPOINT:
	$(CC) $(CFLAGS) $(INC) two_point.cpp $(LIBS) -o TWOPOINT
	mv TWOPOINT ../TWOPOINT

TWOPOINT_MPI:
	$(CCOMP) $(CFLAGS) $(INC) two_point_mpi.cpp $(LIBS) -o TWOPOINT_MPI
	mv TWOPOINT_MPI ../TWOPOINT_MPI

clean:
	rm -r ../TWOPOINT
	rm -r ../TWOPOINT_MPI
