g++ ../src/Disco.cpp
./a.out 1 mag.gif > cm.g
./a.out 0 nmag.gif > sm.g

gnuplot cm.g
gnuplot sm.g

mv mag.gif ../outp/
mv nmag.gif ../outp/

rm cm.g sm.g a.out
