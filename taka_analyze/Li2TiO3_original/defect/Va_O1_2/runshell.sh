#!/bin/bash 
#============ LSF Options =============
#QSUB -q gr10261b 
#QSUB -W 336:00 
#QSUB -A p=16:t=1:c=1:m=3840M 
#QSUB -rn
#============ Shell Script ============

zsh ~/script/shell/repeatGeomOpt5.sh mpiexec.hydra  /home/b/b32282/bin/vasp/vasp533GGA
