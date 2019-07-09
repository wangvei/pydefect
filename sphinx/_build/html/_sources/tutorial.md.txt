# Simple tutorial of pydefect
-----------------------

This web page illustrates how to use the PyDefect code.

Usually, point-defect calculations are complicated following the flow below.

![](flowchart.png)

As shown, some tasks are performed concurrently.

Here, we suppose the following directory tree.

```
    <project_name>
     │
     ├ unitcell/ ── structure_opt/
     │            ├ band/
     │            └ dos/
     │            
     ├ competing_phases/ ── <competing_phase 1>
     │                    ├ <competing_phase 2>
     │                    ....
     │  
     └ defects/ ── perfect/  
                 ├ Va_X_0/
                 ├ Va_X_1/
                 ├ Va_X_2/
                 ...
```
Details of the process are examined step by step.


### 1. Relaxation of unit cell

The point defect calculations need to be generally performed at a fully relaxed theoretical structures, 
including both lattice constants and fractional coordinates of atomic positions. 
Therefore, one begins with fully relaxing the unit cell.

First, we need to prepare the POSCAR file of the pristine bulk unitcell.
And, create `unitcell/` directory and `structure_opt` subdirectory under `unitcell/` and move there.
PyDefect can construct the vasp input sets, namely INCAR, POTCAR, KPOINTS, for various tasks related to the point-defect calculations with simple command line tools.
For instance, vasp input set for structure optimization of the unitcell with the PBEsol functional can be generated using the following command.

```
python ~/my_bin/pydefect/pydefect/main.py vos -t structure_opt -x pbesol
```

vos, meaning an abbreviation for vasp_oba_set, is a subparse option.
VaspObaSet is a program that create vasp input files, analyze and visualize vasp output files.
It is written on the basis of [pymatgen](http://pymatgen.org), a robust, open-source Python library for materials analysis. 
Please also refer the document for VaspObaSet.

There are many options, so please refer the help of the sub command by 

```
python ~/my_bin/pydefect/pydefect/main.py vos -h
```

### 2. Calculation of band, DOS, and dielectric tensor

We then calculate the electronic structure, namely band structure, density of states, and dielectric constant.
In the defect calculations, we can determine the valence band maximum (VBM) and 
conduction band minimum (CBM) from the band structure and DOS.

First, make `band/`, `dos/` and `dielectric/` directories under `unitcell/`.
And then move to each directory and type the following command, respectively.
```
python ~/my_bin/pydefect/pydefect/main.py vos -t band -x pbesol
```
```
python ~/my_bin/pydefect/pydefect/main.py vos -t dos -x pbesol
```
```
python ~/my_bin/pydefect/pydefect/main.py vos -t dielectric -x pbesol
```

After finishing the vasp calculations, type
```
python ~/my_bin/phos_pbesol/obadb/obadb/analyzer/main.py pb -v vasprun.xml -k KPOINTS -f band.pdf
```
and
```
python ~/my_bin/phos_pbesol/obadb/obadb/analyzer/main.py pd -v vasprun.xml -f dos.pdf
```
for making the plot of the band structure and DOS.

Here, the band path is determined based upon the [seekpath code](https://www.materialscloud.org/work/tools/seekpath), 
so if one uses the plot, please cite the following paper.
- [Y. Hinuma, G. Pizzi, Y. Kumagai, F. Oba, I. Tanaka, Band structure diagram paths based on crystallography, Comp. Mat. Sci. 128, 140 (2017).](https://www.sciencedirect.com/science/article/pii/S0927025616305110?via%3Dihub) 
  DOI: 10.1016/j.commatsci.2016.10.015 (the "HPKOT" paper; arXiv version: arXiv:1602.06402).

### 3. Gathering unitcell information related to point-defect calculations

At this point, we can collect the unitcell information that is essential for analyzing point-defect calculations,
namely band edges, electronic and ionic contributions to dielectric tensor, unitcell volume, and unitcell DOS,
using the unitcell_results or ur subparse option.
```
python ~/my_bin/pydefect/pydefect/main.py ur --static_diele_dir dielectric --ionic_diele_dir dielectric --band_edge_dir band --volume_dir dielectric --total_dos_dir dos --json_file unitcell.json -o OUTCAR -p CONTCAR -v vasprun.xml
```
Then, unitcell.json is generated.
One can see the information with the command.
```
python ~/my_bin/pydefect/pydefect/main.py ur --print
```
which shows as follows.

```
vbm (eV): 3.4103
cbm (eV): 4.4278
static dielectric tensor:
[[ 4.680097 -0.       -0.      ]
 [-0.        4.680097  0.      ]
 [ 0.       -0.        4.680097]]
ionic dielectric tensor:
[[14.953044 -0.        0.      ]
 [-0.       14.953065 -0.      ]
 [ 0.       -0.       14.953038]]
total dielectric tensor:
[[19.633141000000002, -0.0, 0.0], [-0.0, 19.633162, 0.0], [0.0, -0.0, 19.633135]]
volume (A^3): 68.98647886527607
Total DOS: Exists
```

### 4. Calculation of competing phases
We then calculate the competing phases. 

First, we prepare their calculation directories including VASP input sets.
By using the chempotdiag program developed and managed by Akira Takahashi, we can obtain structures of the stable or slightly unstable competing phases.
```
python ~/my_bin/phos_pbesol/obadb/obadb/analyzer/chempotdiag/main.py cpd -m -el B N -ch 0.05
```
In this case, the materials of which energy above hull is less than 0.5 meV/atom are collected.
We then generate the other input files.    
```
python ~/my_bin/pydefect/pydefect/main.py vos -t structure_opt -x pbesol -is ENCUT 520 --dirs *_*/
```

Note, if competing phases are gases, we need to change the ISIF to 2 to fix the lattice constants.
In such case, like N2molecule_pydefect, type as follows,
```
python ~/my_bin/pydefect/pydefect/main.py vos -t structure_opt -x pbesol -is ENCUT 520 ISIF 2 --dirs N2molecule_pydefect
```

Next, we generate the chemical potential diagram.
```
python ~/my_bin/phos_pbesol/obadb/obadb/analyzer/chempotdiag/main.py cpd -v */ -p POSCAR-finish -o OUTCAR-finish -c BN -y -s cpd.pdf
```
Note that there are some parameters for gas phases, so read help for more details.

### 5. Select of a supercell
So far, we have done the calculations of the unit cell and competing phases.
So, let's move on to defect/ and copy unitcell POSCAR file from unitcell/dos/ to defect/

Firstly, we need to determine the supercell.
pydefect recommends a supercell which usually have proper isotropy with moderate number of atoms.
For this purpose, type
```
python ~/my_bin/pydefect/pydefect/main.py rs
```

### 6. Construction of defect initial setting
We then generate the defect initial setting file by typing 
```
python ~/my_bin/pydefect/pydefect/main.py is
```

With this command, one can build the `defect.in` file, which contains the information which kind of defects one want to generate.



### 7. Decision of interstitial sites
After constructing an input for standard defect set, one may want to add the interstitial-type defects.
For this purpose, we need to determine the interstitial sites.
Most people determine the interstitial sites by seeing the host crystal structures.
On the contrary, there are a couple of procedures that recommend the interstitial sites.
Still, it is difficult to predict the most likely interstitial sites because they strongly depend on the substituted atoms.
For instance, when positively charged cations with closed shells are substituted (e.g., Mg<sup>2+</sup>, Al<sup>3+</sup>), 
the interstitial sites with the largest vacant space should be most likely. 
On the other hand, when a proton (H<sup>+</sup>) is inserted, it should prefer to locate near O<sup>2-</sup> or N<sup>3-</sup> to form the strong O-H or N-H bonding.
Therefore, we need to be careful when determining the interstitial sites.

To add the interstitial site, e.g., 0.375 0.375 0.375, we type
```
python ~/my_bin/pydefect/pydefect/main.py i -c 0.375 0.375 0.375
```
The interstitials.yaml is then generated, which show various information related to the interstitial sites.
```
i1:
  representative_coords: [0.375, 0.375, 0.375]
  wyckoff: a
  site_symmetry: -43m
  symmetry_multiplicity: 64
  coordination_distances:
    Mg: [1.84, 1.84, 1.84, 1.84]
    O: [1.84, 1.84, 1.84, 1.84, 3.52, 3.52, 3.52, 3.52, 3.52, 3.52, 3.52, 3.52, 3.52,
      3.52, 3.52, 3.52]
  method: manual
```

If you try to add the site that is very close to the constituent atoms or other interstitial sites,
you will get the error message as 

```
Inserted position is too close to another interstitial site.
The distance is 0.042 A.
```
If you really want to add the site even in such cases, add --force_add option.

Once you generate the interstitial.yaml, you also need to modify the defect.in file.
When you modify it, add the following line to the defect.in.
```
Interstitials: i1
```
Or you can again type as follows again.
```
python ~/my_bin/pydefect/pydefect/main.py is --interstitial_sites i1
```

### 8. Construction of defect calculation directories
After constructing the defect.in file, we 


### 9. Generation of supercell information

### 10. Correction of defect formation energy in finite-size supercell 

### 11. Check of defect eigenvalues

### 12. Plot of defect formation energies
