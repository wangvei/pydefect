#!/bin/zsh
#$ -S /bin/zsh
#$ -cwd
#$ -V
#$ -j y
#$ -N pydefect
#$ -o std.log
#$ -pe all_pe* 16
#============ Shell Script ============

zsh /home/kuma/my_bin/repeatGeomOpt5.sh /home/kuma/bin/openmpi-1.8.1-intel16.0.2/bin/mpirun -np $NSLOTS /home/kuma/bin/vasp.5.4.4/bin/vasp_std
