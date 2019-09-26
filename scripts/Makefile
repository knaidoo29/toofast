MACHINE = KRISHMACBOOK

ifeq ($(MACHINE),KRISHMACBOOK)
	CC = g++
	CCOMP = /Users/krishna/Programs/build/openmpi/bin/mpic++
	TWOFAST_INC = -I/Users/krishna/Research/github/twofast/include
	TWOFAST_LIB = -L/Users/krishna/Research/github/twofast/lib -ltwofast
endif

ifeq ($(MACHINE),SPLINTER)
	CC = g++
	CCOMP = mpic++
	TWOFAST_INC = -I/share/data1/knaidoo/splinter_libraries/src/twofast/include
	TWOFAST_LIB = -L/share/data1/knaidoo/splinter_libraries/src/twofast/lib -ltwofast
endif

OPTIMIZE = -O3

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
