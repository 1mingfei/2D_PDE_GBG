set terminal postsc eps enhanced color font ",25"
set output "200.eps"

set xlabel "x(nm)"
set ylabel "z(nm)"
set yrange [-50:50]
set xrange [0:3000]

plot '200.out'  using ($1*1e7):($2*1e7) w l lw 4 smooth csplines title "200ms",\
