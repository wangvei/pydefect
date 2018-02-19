#!/bin/zsh

defect_position="0.81928147980973      0.79704247333310      0.63373526153701" #Position of the defect in fractional coordinate or atomic number.
charge_state=2 # [|e|]
epsilon="42.396476 0 -1.057006 0 42.396476 0 -1.057006 0 15.812507"
accuracy=25
omega=0.01 # its original value
#omega=0.000874491269411 # debug
#check="--check"
#EXTRINSIC=--remove
OUTCAR=OUTCAR-final
POSCAR=POSCAR-final
script_dir="../../defect_correction/"

step=`echo {$1..$2}`

for i in `echo $step`; do
echo step$i
case $i in
1) 
python $script_dir/get_potential_vasp.py -o $OUTCAR -r reference_potential.txt -d $defect_position >| vasp_potential.txt
;;
2) 
python $script_dir/defect_structure.py -p $POSCAR -d $defect_position >| defect_structure.txt
;;
3) 
python $script_dir/potential.py -p $POSCAR -c $charge_state -d $defect_position --accuracy $accuracy -e $epsilon --omega $omega $check >| model_potential.txt
;;
4) 
TMPFILE=`mktemp`; TMPFILE2=`mktemp`
line=`{wc -l defect_structure.txt | awk '{print $1}'}`; tail -n $((line-3)) defect_structure.txt >| $TMPFILE
line=`{wc -l model_potential.txt | awk '{print $1}'}`; head -n $((line-1)) model_potential.txt >| $TMPFILE2
paste $TMPFILE vasp_potential.txt $TMPFILE2 > potential.txt
;;
5)
distance=`{python ~/defect_correction/parse_poscar.py -p POSCAR-finish -r | awk '{print $11}'}`
python $script_dir/parse_potential.py --potential potential.txt  -d $distance -c $charge_state >| alignment.txt
potential_shift=`{head -n 1 alignment.txt | awk '{print $3}'}`
sed s/distance/$distance/ ..//plot.gp |sed s/potential_shift/$potential_shift/ > plot.gp
;;
6)
grep Li potential.txt >| potential_Li
grep Ti potential.txt >| potential_Ti
grep O potential.txt >| potential_O
;;
7)
gnuplot plot.gp
;;
esac
done
