#===General Variables===#
CC=gcc
CFLAGS=-Wall -Wextra -g3

all: makeAll

makeAll: makeHelper makeRing makePoly makeAlgebra makeDSP makeFilter makeMain
	$(CC) $(CFLAGS) helper.o ring.o polynomial.o algebra.o dsp.o filter.o main.o -o dsp -ldl -lm -lblas -llapack -lpthread 

makeMain: main.c
	$(CC) $(CFLAGS) -c main.c -o main.o 

makeFilter: filter.c filter.h
	$(CC) $(CFLAGS) -c filter.c -o filter.o

makeDSP: dsp.c dsp.h
	$(CC) $(CFLAGS) -c dsp.c -o dsp.o

makeAlgebra: linear_algebra.c linear_algebra.h
	$(CC) $(CFLAGS) -c linear_algebra.c -o algebra.o

makePoly: polynomial.c polynomial.h
	$(CC) $(CFLAGS) -c polynomial.c -o polynomial.o

makeRing: ring.c ring.h
	$(CC) $(CFLAGS) -c ring.c -o ring.o

makeHelper: helper.c helper.h macros.h
	$(CC) $(CFLAGS) -c helper.c -o helper.o

.PHONY: clean

clean:
	rm -f *~ $(ODIR)/*.o $(IDIR)/*~
