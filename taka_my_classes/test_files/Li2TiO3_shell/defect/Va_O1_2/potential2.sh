#!/bin/zsh

defect_position="0.81928147980973      0.79704247333310      0.63373526153701" #Position of the defect in fractional coordinate or atomic number.
charge_state=2 # [|e|]
epsilon="42.396476 0 -1.057006 0 42.396476 0 -1.057006 0 15.812507"
accuracy=25
omega=0.01 # its original value
#EXTRINSIC=--remove
OUTCAR=OUTCAR-final
POSCAR=POSCAR-final
script_dir="../../defect_correction/"
plot_dir=../
species=`echo Mg O`

step=`echo {$1..$2}`

for i in `echo $step`; do
echo step$i
case $i in
0)
python $script_dir/potential.py -p $POSCAR --accuracy $accuracy -e $epsilon $omega --check  --omega $omega
;;
1) 
python $script_dir/get_potential_vasp.py -o $OUTCAR -r reference_potential.txt -d $defect_position >| vasp_potential.txt
;;
2) 
python $script_dir/defect_structure.py -p $POSCAR -d $defect_position >| defect_structure.txt
;;
3) 
python $script_dir/potential.py -p $POSCAR -c $charge_state -d $defect_position --accuracy $accuracy -e $epsilon --omega $omega >| model_potential.txt
;;
4) 
TMPFILE=`mktemp`; TMPFILE2=`mktemp`
line=`{wc -l defect_structure.txt | awk '{print $1}'}`; tail -n $((line-3)) defect_structure.txt >| $TMPFILE
line=`{wc -l model_potential.txt | awk '{print $1}'}`; head -n $((line-1)) model_potential.txt >| $TMPFILE2
paste $TMPFILE vasp_potential.txt $TMPFILE2 > potential.txt
;;
5)
distance=`{python $script_dir/parse_poscar.py -p $POSCAR -r | awk '{print $11}'}`
python $script_dir/parse_potential.py --potential potential.txt  -d $distance -c $charge_state >| alignment.txt
potential_shift=`{head -n 1 alignment.txt | awk '{print $3}'}`
sed s/distance/$distance/ $plot_dir/plot.gp |sed s/potential_shift/$potential_shift/ >| plot.gp
;;
6)
for i in `echo $species`
do 
  grep $i potential.txt >| potential_$i
  echo "     \"potential_$i\" using 6:7 w p,\\" >> plot.gp
done

echo '     "potential.txt" using 6:8, "" using 6:($7-$8)' >> plot.gp
gnuplot plot.gp
;;
esac
done
