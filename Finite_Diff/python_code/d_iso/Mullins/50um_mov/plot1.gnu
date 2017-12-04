set terminal postsc eps enhanced color
set output "tmp.eps"

set xlabel "x(cm)"
set ylabel "z(cm)"
set xrange [0:2.5e-4]

plot '0.out'  using 1:2 w l smooth csplines title "0.1s",\
