CC = g++
#CC = /Users/krishna/Programs/build/openmpi/bin/mpic++
CFLAGS = -c -O3

INC =
LIBS =

OBJS = src/two_point.o src/util_omp.o src/util_param.o src/util_progress.o src/util_read.o src/util_write.o

all: lib/libtwofast.a

lib/libtwofast.a: $(OBJS)
	 ar cr lib/libtwofast.a $(OBJS)

.cpp.o:
	 $(CC) $(CFLAGS) $(INC) $< $(LIBS) -o $@

clean:
	 rm src/*.o*
	 rm lib/*.a*
	 rm .DS_Store
