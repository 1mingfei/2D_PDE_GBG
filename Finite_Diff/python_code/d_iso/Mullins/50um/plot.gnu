set terminal postsc eps enhanced color font ",25"
set output "tmp.eps"

set xlabel "x(nm)"
set ylabel "z(nm)"
set yrange [-50:20]
set xrange [0:3000]


plot '100.out'  using ($1*1e7):($2*1e7)  w l lw 4 smooth csplines title "0.1s",\
     '1000.out'  using ($1*1e7):($2*1e7) w l lw 4 smooth csplines title "1s",\

