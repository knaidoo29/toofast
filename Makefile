CC = g++
CFLAGS = -c

INC =
LIBS =

OBJS = src/progress.o src/read.o src/write.o

all: lib/libtoofast.a

lib/libtoofast.a: $(OBJS)
	 ar cr lib/libtoofast.a $(OBJS)

.cpp.o:
	 $(CC) $(CFLAGS) $(INC) $< $(LIBS) -o $@

clean:
	 rm src/*.o*
	 rm lib/*.a*
	 rm .DS_Store
