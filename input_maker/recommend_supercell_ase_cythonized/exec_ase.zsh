#!/usr/bin/zsh

python setup.py build_ext -i && python ase_generation_supercell.py POSCAR-unitcell
