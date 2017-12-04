#!~/usr/bash
for i in `seq 200 200 5000`
do
    echo $i
    s=`echo "scale=1; $i/1000" |bc`
    cp plot1.gnu.ref plot.gnu
    sed -i '.bak' s/lala_i/$i/g plot.gnu
    sed -i '.bak' s/lala_s/$s/g plot.gnu
    gnuplot plot.gnu
done
#for i in `seq 2000 2000 10000`
#do
#    echo $i
#    s=`echo "scale=1; $i/1000" |bc`
#    cp plot1.gnu.ref plot.gnu
#    sed -i '.bak' s/lala_i/$i/g plot.gnu
#    sed -i '.bak' s/lala_s/$s/g plot.gnu
#    gnuplot plot.gnu
#done
