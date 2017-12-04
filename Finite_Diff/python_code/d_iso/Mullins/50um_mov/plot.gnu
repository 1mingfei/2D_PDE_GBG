set terminal postsc eps enhanced color font ",25"
set output "tmp.eps"

set xlabel "x(cm)"
set ylabel "z(cm)"

#plot '100.out'  using 1:2 w l smooth csplines title "0.1s",\
#     '1000.out'  using 1:2 w l smooth csplines title "1s",\

plot for [i=0:200:20] ''.i.'.out' using 1:2 w l smooth csplines  title ''.i.'ms',\
