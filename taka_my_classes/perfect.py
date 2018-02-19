# TODO: Make unittest
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.outputs import Outcar, Vasprun


class Perfect:

    def __init__(self, directory_name,
                 poscar_name = "POSCAR",
                 outcar_name = "OUTCAR",
                 vasprun_name = "vasprun.xml"):
        path_poscar = directory_name + poscar_name
        path_outcar = directory_name + outcar_name
        path_vasprun = directory_name + vasprun_name
        poscar = Poscar.from_file(path_poscar)
        outcar = Outcar(path_outcar)
        vasprun = Vasprun(path_vasprun)
        self._structure = poscar.structure
        self._energy = outcar.final_energy
        self._dielectric_tensor \
            = outcar.dielectric_tensor + outcar.dielectric_ionic_tensor
        # TODO: eigenvalue is right?
        # (Correct: EIGENVAL file. check if implemented and EIGENVAL is same )
        self._eigen_value = vasprun.eigenvalues
        self._electrostatic_potential = outcar.electrostatic_potential

    @property
    def structure(self):
        return self._structure

    @property
    def energy(self):
        return self._energy

    @property
    def dielectric_tensor(self):
        return self._dielectric_tensor

    @property
    def eigen_value(self):
        return self._eigen_value

    @property
    def electrostatic_potential(self):
        return self._electrostatic_potential
