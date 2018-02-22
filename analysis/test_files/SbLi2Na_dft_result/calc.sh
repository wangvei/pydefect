#!/bin/zsh
#$ -S /bin/zsh
#$ -cwd
#$ -V
#$ -j y
#$ -p -1000
#$ -o std.log
#$ -pe all_pe* 16
#$ -l h_rt=72:00:00
#$ -N dielectric

repeat_script_dir="/home/taka/dielectric_tensor/icsd/pbesol_dielectric_constant/"
mpi_exec="/home/common/bin/openmpi-1.8.8_intel-16.0.2/bin/mpirun"
vasp_exec="/home/common/bin/vasp_std"

zsh $repeat_script_dir/repeatGeomOpt5.sh $mpi_exec $vasp_exec 

