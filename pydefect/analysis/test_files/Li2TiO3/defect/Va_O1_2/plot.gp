set terminal postscript eps color enhanced "Arial" 20
set encoding iso 

unset key
set xzeroaxis

step(x,from,y)=( (from <= x) ) ? y : 0/0

set xlabel "interval from a defect (\305)"
set ylabel "potential (eV)"

set yrange [-0.5:0.5]

set output "potential.eps"
plot step(x, distance, potential_shift),\
     "potential_Li" using 6:7 w p,\
     "potential_Ti" using 6:7 w p,\
     "potential_O"  using 6:7 w p,\
     "potential.txt" using 6:8, "" using 6:($7-$8)
