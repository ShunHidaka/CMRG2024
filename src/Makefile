CC = gcc
CFLAGS    = -lm -g -Wall -Wfloat-equal
FASTFLAGS = -O3
LAFLAGS   = -lgfortran -lblas -llapack
#LAFLAGS   = -I/opt/OpenBLAS/include/ -L/opt/OpenBLAS/lib -lopenblas

all: gsminres.out gscocg.out inner.out

gsminres.out: gsminres.c function_util.c
	$(CC) $^ -o $@ $(CFLAGS) $(FASTFLAGS) $(LAFLAGS)
gscocg.out: gscocg.c function_util.c
	$(CC) $^ -o $@ $(CFLAGS) $(FASTFLAGS) $(LAFLAGS)
inner.out: inner.c function_util.c
	$(CC) $^ -o $@ $(CFLAGS) $(FASTFLAGS) $(LAFLAGS)

run1: all
	./gsminres.out 1> stdout_gsminres.dat 2> stderr_gsminres.dat
	./gscocg.out   1> stdout_gscocg.dat   2> stderr_gscocg.dat
	gnuplot graph1.gp
run2: inner.out
	./inner.out 0 1> stdout_inner0.dat
	./inner.out 1 1> stdout_inner1.dat
	./inner.out 2 1> stdout_inner2.dat
	./inner.out 3 1> stdout_inner3.dat
	gnuplot graph2.gp
clean:
	rm -f *.out *.dat *.pdf
