0. prerequisite: 
    a. Unit-cell calculations with 
        a-1. unit-cell lattice being optimized
        a-2. band-edge positions from band structure calculation for the defect formation energies
        a-3. accurate DOS for the carrier concentration analysis
    b. Relaxed POSCAR for a perfect supercell.
    c. INCAR and KPOINTS for defect calculations.
       INCAR must not include NELECT that will be determined automatically.
    d. ~/.pydefect.yaml that contains the information of director for 
       default potcars. An example is located in this directory.
       e.g., DEFAULT_POTCAR: /home/common/default_POTCAR
    e. The DEFAULT_POTCAR directory contains POTCAR files named as "POTCAR_" + element name. 
       e.g.,
        /home/common/default_POTCAR/POTCAR_H
        /home/common/default_POTCAR/POTCAR_He
        /home/common/default_POTCAR/POTCAR_Li ...

    An example of initial file configurations
       MgO_defect/       <---- parent directory
       MgO_defect/{INCAR,KPOINTS,POSCAR}
       ~/.pydefect.yaml  <---- home directory

1. Make defect.in file from a POSCAR that is a supercell.
    % python defect_input.py  --help
------------------------------------------------------------------------------
  -h, --help            show this help message and exit
  -p POSCAR, --poscar POSCAR
                        POSCAR name. (default: POSCAR)
  -d DOPANTS [DOPANTS ...], --dopants DOPANTS [DOPANTS ...]
                        Dopant elements. Eg. Al Ga In. (default: )
  -i INTERSTITIAL_COORDS [INTERSTITIAL_COORDS ...]
                        Inetrstitials. Eg. 0 0 0 0.5 0.5 0.5. (default: False)
  -a, --antisite        Set if antisites are considered. (default: True)
  -e ELNEG_DIFF, --ElNeg_diff ELNEG_DIFF
                        Criterion of the electronegativity difference for
                        constructing antisite and impurities. (default: 1.0)
  --include INCLUDE     Exceptionally included defect type. Eg Va_O2_-1.
                        (default: )
  --exclude EXCLUDE     Exceptionally excluded defect type. Eg Va_O2_0.
                        (default: )
  -s, --symbreak        Set if symmetry is not broken. (default: False)
  --displace DISPLACE   Displacement distance. (default: 0.2)
  --symprec SYMPREC     Set the symprec [A]. (default: 0.01)
  --print_dopant PRINT_DOPANT
                        Print Dopant information. (default: None)
------------------------------------------------------------------------------
    % python default_input.py -p POSCAR  xxxxxxxx
    This command generates DPOSCAR and the following defect.in file.
    Note that DPOSCAR is a modified POSCAR file, in which sequence of atoms 
    is sorted by element and site symmetry. Example is POSCAR-xx and DPOSCAR-xx
    defect.in file contains a full information of point defect calculations.
    The user is allowed to modify the input depending on the purpose.
------------------------------------------------------------------------------
  Name: Mg1
   Rep: 1
 Equiv: 1..32
 Coord: 0.0000000 0.0000000 0.0000000
EleNeg: 1.31
Charge: 2

  Name: O1
   Rep: 33
 Equiv: 33..64
 Coord: 0.2500000 0.2500000 0.2500000
EleNeg: 3.44
Charge: -2

Antisite: Mg_O1 O_Mg1

Symbreak: 0.15
Include: Va_O1_-1 Va_O1_-2
Exclude: Va_O1_1 Va_O1_2
Symprec: 0.001

Int_site: 0.1 0.1 0.1
------------------------------------------------------------------------------

2. Make directories for the point defect calculations. 
   Each directory's name *xx_yy_zz* contains information of the defect type.
        xx : Inserted element name. 
             "Va" means vacancy or no inserted element.
        yy : Removed irreducible atom, e.g. Mg1. 
             "iN", where N is an integer,  means Nth interstitial site.
        zz : Charge state.
    The name of the perfect calculation is *perfect*. 

   Thus, the directory names for this example are as following:
    - MgO_defect/perfect
    - MgO_defect/Va_Mg1_0
    - MgO_defect/Va_Mg1_-1
    - MgO_defect/Va_Mg1_-2
    - MgO_defect/Va_O1_0
....

    Each directory contains INCAR POSCAR POTCAR KPOINTS.
     - INCAR is copied from a parent directory. Then, NELECT will be added based on the defect type and charge state.
     - POSCAR is constructed based on the representative atoms written in defect.in and POSCAR in parent directory.
       First line for comments denotes the defect position. 
     - POTCAR is constructed based on element names written at 6th line in DPOSCAR file.
       The default potcar directory path is set in .pydefect.yaml.
     - KPOINTS is simply copied from a parent directory.

    Warning:
        1. When incorrect element name(s) are used, assert IncorrectElementNameError.

------------------------------------------------------------------------------
3. Run vasp for each directory.

    After calculations, the following three files are needed for postprocess.
    - POSCAR-initial: initial structure file
    - POSCAR-final: final structure file
    - OUTCAR-final: OUTCAR file including the converged results 

------------------------------------------------------------------------------
4. Analyze calculation data.
    a. Correct defect formation energies.
        a-1. Construct a full information JSON file of the defect in each directory.
    b. Quantify local structures and make local structure files for visualization.
        b-1. Check the change of local structure via structural optimization.
        b-2. Plot the square of (sums of) one-electron wavefunction(s). 
            --> Construct vasp input files automatically
    c. Judge if the defect contains shallow defect state.
    d. Gather total and correction energies and plot energy vs Fermi level. 
        d-1. Calculate the Fermi level and carrier concentration at room temperature.
        d-2. Calculate the Fermi level and carrier concentration quenched from a growth temperature to room temperature.
    e. Defect eigenvalues with respect to the band edge positions.
 
------------------------------------------------------------------------------
Future Works: 
# A supercell should be recommended and/or automatically constructed from a 
  unit cell. 
# Interstitial site(s) should be recommended and/or automatically constructed 
  from a unit cell.
# Input of cell-size dependence should be automatically constructed. 
# A unit cell calculation may be also included if needed by option.
