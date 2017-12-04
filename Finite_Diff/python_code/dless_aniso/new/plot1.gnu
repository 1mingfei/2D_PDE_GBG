set terminal postsc eps enhanced color font ",25"
set output "tmp1.eps"

set xlabel "x (dimensionless)"
set ylabel "y (dimensionless)"
set yrange [-0.1:0.02]
unset xtics
unset ytics

plot '100.iso.out'  using 1:2 w l lw 4 smooth csplines notitle,\
