#!/bin/zsh

defect_position="0 0 0" #Position of the defect in fractional coordinate or number of defect site.
charge_state=-4 # [|e|]
epsilon="14.308715 14.308715 13.846087"
accuracy=25
#omega=
#EXTRINSIC=--remove
OUTCAR=OUTCAR
POSCAR=CONTCAR
script_dir=/home/b/b32282/defect_correction/PyDef-0.1 # path to the directory containing PyDef scripts.
plot_dir=../ # path to the directory containing plot.gp
species=`echo Zn Sn P`

step=`echo {$1..$2}`

for i in `echo $step`; do
echo step$i
case $i in
0)
python $script_dir/potential.py -p $POSCAR --accuracy $accuracy -e $epsilon $omega --check 
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
