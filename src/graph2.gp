fname0 = "stdout_inner0.dat"
fname1 = "stdout_inner1.dat"
fname2 = "stdout_inner2.dat"
fname3 = "stdout_inner3.dat"

set terminal pdfcairo
set out "compare-inner.pdf"
set grid
set colorsequence default
set xlabel "index of shifts"

set ylabel "number of iterations"
set yrange [0:]
plot   fname0 u 1:4 w l t "1e-13",\
       fname1 u 1:4 w l t "1e-12",\
       fname2 u 1:4 w l t "1e-10",\
       fname3 u 1:4 w l t "1e-8"

set yrange [1e-13:1]
set logscale y
set ylabel "true relative residual"
plot   fname0 u 1:6 w l t "1e-13",\
       fname1 u 1:6 w l t "1e-12",\
       fname2 u 1:6 w l t "1e-10",\
       fname3 u 1:6 w l t "1e-8"