#!~/usr/bash
for i in `seq 0 1 20`
do
    echo $i
    cp plot1.gnu.ref plot.gnu
    sed -i '.bak' s/lala_i/$i/g plot.gnu
    gnuplot plot.gnu
done
for i in `seq 0 20 200`
do
    echo $i
    cp plot1.gnu.ref plot.gnu
    sed -i '.bak' s/lala_i/$i/g plot.gnu
    gnuplot plot.gnu
done
