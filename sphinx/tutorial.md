# Simple tutorial of pydefect
-----------------------

Point-defect calculations are bit complicated and performed following the flow shown below.

<!--
![flow](flowchart.png)
-->

Here, we suppose the following directory tree.

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

### 1. Relaxation of unit cell

The point defect calculations need to be performed at a relaxed theoretical atomic structures.
Therefore, one needs to fully relax the atomic position of the unit cell.
PyDefect allows us to construct the vasp input file set, namely INCAR, POTCAR,
KPOINTS, with rather simple command line tool.
For instance, initial setting for structure optimization of the unit-cell with the PBEsol functional can be
generated using the following command.

```
python ~/my_bin/pydefect/pydefect/main.py vos -t structure_opt -x pbesol
```

There are many options, so please refer the help of the sub command.

### 2. Calculation of band, DOS, and dielectric tensor

We then calculate the electronic structure, namely band and density of states (DOS) plot, and dielectric constant.
In the defect calculations, we can determine the valence band maximum (VBM) and 
conduction band minimum (CBM) from the band and DOS.

First, make band/, dos/ and dielectric and type the following commands, respectively.
```
python ~/my_bin/pydefect/pydefect/main.py vos -t band -x pbesol
```
```
python ~/my_bin/pydefect/pydefect/main.py vos -t dos -x pbesol
```
```
python ~/my_bin/pydefect/pydefect/main.py vos -t dielectric -x pbesol
```

After finishing the vasp calculations, run 
```
python ~/my_bin/phos_pbesol/obadb/obadb/analyzer/main.py pb -v vasprun.xml -k KPOINTS -f band.pdf
```
and
```
python ~/my_bin/phos_pbesol/obadb/obadb/analyzer/main.py pd -v vasprun.xml -f dos.pdf
```
for making the plot of the band structure and DOS.

### 3. Collection of unitcell information related to point-defect calculations

At this point, we can collect the unitcell information related to point-defect calculations,
namely electronic and ionic contributions to dielectric tensor, unitcell volume, unitcell DOS. 
```
python ~/my_bin/pydefect/pydefect/main.py ur --static_diele_dir dielectric --ionic_diele_dir dielectric --band_edge_dir band --volume_dir dielectric --total_dos_dir dos --json_file unitcell.json -o OUTCAR-finish -p POSCAR-finish -v repeat-1/vasprun.xml
```


### 4. Calculation of competing phases
We then calculate the competing phases. 
First, we prepare their calculation directories including VASP input files.
By using the chempotdiag program, we can obtain structures of the stable or slightly unstable competing phases.
```
python ~/my_bin/phos_pbesol/obadb/obadb/analyzer/chempotdiag/main.py cpd -m -el B N -ch 0.05
```

Then, we generate the other input files.    
```
python ~/my_bin/pydefect/pydefect/main.py vos -t structure_opt -x pbesol -is ENCUT 520 --dirs *_*/
```

Note, if competing phases are gas phases, we need to change the ISIF to 2.
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

With this command, one can 


### 7. Decision of interstitial sites
After constructing an input for standard defect set, one may want to add the interstitial-type defects.
For this purpose, we need to determine the interstitial sites.
Most people determine the interstitial sites by seeing the host crystal structures.
On the contrary, there are a couple of procedures that recommend the interstitial sites.
Still, it is difficult to predict the most likely interstitial sites because they strongly depend on the substituted atoms.
For instance, when positively charged cations with closed shells are substituted (e.g., Mg2+, Al3+), 
the interstitial sites with the largest vacant space should be most likely. 
On the other hand, when a proton (H+) is inserted, it should prefer to locate near O2- or N3- to form the strong O-H or N-H bonding.
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



### 9. Generation of supercell information

### 10. Correction of defect formation energy in finite-size supercell 

### 11. Check of defect eigenvalues

### 12. Plot of defect formation energies
