CC = /Users/krishna/Programs/build/openmpi/bin/mpic++

OPTIMIZE = -O3

TWOFAST_INC = -I/Users/krishna/Research/github/twofast/include
TWOFAST_LIB = -L/Users/krishna/Research/github/twofast/lib -ltwofast

INC = $(TWOFAST_INC)
LIBS = -lm $(TWOFAST_LIB)
CFLAGS = $(OPTIMIZE)

OBJS = two_point_omp.o

all: $(OBJS)

.cpp.o:
	$(CC) $(CFLAGS) $(INC) $< $(LIBS) -o $@

clean:
	rm -r *.o*
	rm .DS_Store
