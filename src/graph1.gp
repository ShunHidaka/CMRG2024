fname1 = "stderr_gsminres.dat"
fname2 = "stderr_gscocg.dat"

set terminal pdfcairo
set out "compare-residual.pdf"
set grid
set colorsequence default
set xlabel "number of iterations"
set ylabel "relative residual"
set logscale y

plot   fname1 u 1:2 w l t "real(gsminres)", fname1 u 1:3 w l t "true(gsminres)",\
       fname2 u 1:2 w l t "real(gscocg)",   fname2 u 1:3 w l t "true(gscocg)"
