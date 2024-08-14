
# コンパイラの設定
ifdef fugaku
	CC     = fccpx
	CFLAGS = -lm -Kfast -Kopenmp
	LAFLAG = -SSL2
else
	CC     = gcc
	CFLAGS = -lm -O3 -fopenmp -Wall
	LAFLAG = -lgfortran -lopenblas
#	LAFLAG = -lgfortran -lblas -llapack
#	LAFLAG = -I/opt/OpenBLAS/include/ -L/opt/OpenBLAS/lib -lopenblas
endif

# デバッグ用コンパイルオプション
ifdef debug
	CFLAGS = -lm -O0 -Wall -g3
endif

.PHONY: all clean
all: gsminres.out gscocg.out inner.out

gsminres.out: gsminres.c function_util.c
	$(CC) $^ -o $@ $(CFLAGS) $(LAFLAG)
gscocg.out: gscocg.c function_util.c
	$(CC) $^ -o $@ $(CFLAGS) $(LAFLAG)
inner.out: inner.c function_util.c
	$(CC) $^ -o $@ $(CFLAGS) $(LAFLAG)

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
