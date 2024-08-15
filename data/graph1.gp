fname1 = "run1_gsminres.sh.33974280.err.txt"
fname2 = "run1_gscocg.sh.33974281.err.txt"

set terminal pdfcairo
set out "./../fig/compare-residual.pdf"
set grid
set colorsequence default
set xlabel "number of iterations"
set ylabel "relative residual"
set logscale y


plot   fname1 every 5 u 1:2 w l t "real(gsminres)", fname1 every 5 u 1:3 w l t "true(gsminres)",\
       fname2 every 5 u 1:2 w l t "real(gscocg)",   fname2 every 5 u 1:3 w l t "true(gscocg)"
#plot   fname1 u 1:2 w l t "real(gsminres)", fname1 u 1:3 w l t "true(gsminres)",\
#       fname2 u 1:2 w l t "real(gscocg)",   fname2 u 1:3 w l t "true(gscocg)"
