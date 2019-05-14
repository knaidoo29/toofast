CC = g++

OPTIMIZE =

TOOFAST_INC = -I/Users/krishna/Research/github/toofast/include
TOOFAST_LIB = -L/Users/krishna/Research/github/toofast/lib -ltoofast

INC = $(TOOFAST_INC)
LIBS = -lm $(TOOFAST_LIB)
CFLAGS = $(OPTIMIZE)

OBJS = test1.o test2.o test3.o

all: $(OBJS)

.cpp.o:
	$(CC) $(CFLAGS) $(INC) $< $(LIBS) -o $@

clean:
	rm -r *.o*
	rm .DS_Store
