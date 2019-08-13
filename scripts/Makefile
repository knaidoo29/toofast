CC = g++
CCOMP = /Users/krishna/Programs/build/openmpi/bin/mpic++

OPTIMIZE = -O3

TWOFAST_INC = -I/Users/krishna/Research/github/twofast/include
TWOFAST_LIB = -L/Users/krishna/Research/github/twofast/lib -ltwofast

INC = $(TWOFAST_INC)
LIBS = -lm $(TWOFAST_LIB)
CFLAGS = $(OPTIMIZE)

SRCS = script_two_point_3d.cpp
SRCSOMP = script_two_point_3d_omp.cpp
OBJS = script_two_point_3d.o
OBJSOMP = script_two_point_3d_omp.o

all: $(OBJS) $(OBJSOMP)

.cpp.o:
	$(CC) $(CFLAGS) $(INC) $(SRCS) $(LIBS) -o $(OBJS)
	$(CCOMP) $(CFLAGS) $(INC) $(SRCSOMP) $(LIBS) -o $(OBJSOMP)

clean:
	rm -r *.o*
	rm .DS_Store
