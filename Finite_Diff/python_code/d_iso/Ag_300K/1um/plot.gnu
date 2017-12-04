set terminal postsc eps enhanced color font ",25"
set output "0.eps"

set xlabel "x(nm)"
set ylabel "z(nm)"
set yrange [-1.0:0.2]
set xrange [0:100]

plot '0.out'  using ($1):($2) w l lw 4 smooth csplines title "0.0s",\
