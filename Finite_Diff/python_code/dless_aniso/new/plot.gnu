set terminal postsc eps enhanced color font ",25"
set output "tmp1.eps"

set xlabel "x"
set ylabel "y"

plot '100.iso.out'  using 1:2 w l lw 4 smooth csplines title "iso",\
     '100.ani.out'  using 1:2 w l lw 4 smooth csplines title "ani",\

