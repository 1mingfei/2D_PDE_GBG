set terminal postsc eps enhanced color
set output "tmp.eps"

set xlabel "x(nm)"
set ylabel "z(nm)"

plot '100.out'  using 1:2 w l title "100",\
     '500.out'  using 1:2 w l title "500",\

