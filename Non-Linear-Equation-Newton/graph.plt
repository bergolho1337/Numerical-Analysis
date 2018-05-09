set title "Método Euler Implícito"
set grid
set xlabel "t"
set ylabel "y"
set terminal png
set output 'sistema.png'
plot 'dados.dat' using 1:2 title "y1" w l, 'dados.dat' using 1:3 title "y2" w l
