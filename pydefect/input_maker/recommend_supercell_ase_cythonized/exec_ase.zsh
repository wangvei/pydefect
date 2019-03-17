#!/usr/bin/zsh

python setup.py build_ext -i && python ase_generation_supercell.py POSCAR_yk
#python3 setup.py build_ext -i && python3 ase_generation_supercell.py POSCAR_yk
