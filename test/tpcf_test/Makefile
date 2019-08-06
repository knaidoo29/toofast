CC = g++

OPTIMIZE =

TOOFAST_INC = -I/Users/krishna/Research/github/toofast/include
TOOFAST_LIB = -L/Users/krishna/Research/github/toofast/lib -ltoofast

INC = $(TOOFAST_INC)
LIBS = -lm $(TOOFAST_LIB)
CFLAGS = $(OPTIMIZE)

OBJS = two_point_estimator.o weighted_two_point_estimator.o three_point_estimator.o

all: $(OBJS)

.cpp.o:
	$(CC) $(CFLAGS) $(INC) $< $(LIBS) -o $@

clean:
	rm -r *.o*
	rm .DS_Store
