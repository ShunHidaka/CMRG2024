fname0 = "run2_inner0.sh.33974282.out.txt"
fname1 = "run2_inner1.sh.33974283.out.txt"
fname2 = "run2_inner2.sh.33974284.out.txt"
fname3 = "run2_inner3.sh.33974285.out.txt"

set terminal pdfcairo
set out "./../fig/compare-inner.pdf"
set grid
set colorsequence default
set xlabel "index of shifts"

set ylabel "number of iterations"
set yrange [0:]
plot   fname0 u 1:4 w l t "1e-12",\
       fname1 u 1:4 w l t "1e-10",\
       fname2 u 1:4 w l t "1e-8",\
       fname3 u 1:4 w l t "1e-6"

set yrange [1e-13:1]
set logscale y
set ylabel "true relative residual"
plot   fname0 u 1:6 w l t "1e-12",\
       fname1 u 1:6 w l t "1e-10",\
       fname2 u 1:6 w l t "1e-8",\
       fname3 u 1:6 w l t "1e-6"