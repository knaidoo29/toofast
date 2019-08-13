CC = g++
CCOMP = /Users/krishna/Programs/build/openmpi/bin/mpic++

OPTIMIZE = -O3

TWOFAST_INC = -I/Users/krishna/Research/github/twofast/include
TWOFAST_LIB = -L/Users/krishna/Research/github/twofast/lib -ltwofast

INC = $(TWOFAST_INC)
LIBS = -lm $(TWOFAST_LIB)
CFLAGS = $(OPTIMIZE)

SRCS = script_two_point.cpp
SRCSOMP = script_two_point_omp.cpp
OBJS = TWOPOINT #script_two_point.o
OBJSOMP = TWOPOINT_MPI #script_two_point_omp.o

all: $(OBJS) $(OBJSOMP)

#.cpp.o:
#	$(CC) $(CFLAGS) $(INC) $(SRCS) $(LIBS) -o $(OBJS)
#	$(CCOMP) $(CFLAGS) $(INC) $(SRCSOMP) $(LIBS) -o $(OBJSOMP)

TWOPOINT:
	$(CC) $(CFLAGS) $(INC) $(SRCS) $(LIBS) -o $(OBJS)
	mv $(OBJS) ../$(OBJS)

TWOPOINT_MPI:
	$(CCOMP) $(CFLAGS) $(INC) $(SRCSOMP) $(LIBS) -o $(OBJSOMP)
	mv $(OBJSOMP) ../$(OBJSOMP)

clean:
	rm -r ../TWOPOINT
	rm -r ../TWOPOINT_MPI
