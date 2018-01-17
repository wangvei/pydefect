set terminal postscript eps color enhanced "Arial" 20
set encoding iso 

unset key
set xzeroaxis

step(x,from,y)=( (from <= x) ) ? y : 0/0

set xlabel "Distance from a defect (\305)"
set ylabel "Potential (eV)"

set yrange [-1:1]

set output "potential.eps"
plot step(x, 4.24689, 0.1694742),\
     "potential_Mg" using 6:7 w p,\
     "potential_O" using 6:7 w p,\
     "potential.txt" using 6:8, "" using 6:($7-$8)
