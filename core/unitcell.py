import math
from itertools import product
import numpy as np
import scipy
import scipy.constants as sconst
import scipy.stats.mstats as mstats
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.outputs import Outcar, Vasprun
__author__ = "Akira Takahashi"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Akira Takahashi"
__email__ = "takahashi.akira.36m@gmail.com"
__status__ = "Development"
__date__ = "February 19, 2018"


class Unitcell:
    """
    DFT result of unitcell without any defect.
    """

    def __init__(self, directory_path,
                 poscar_name = "/POSCAR",
                 outcar_name = "/OUTCAR",
                 vasprun_name = "/vasprun.xml"):
        """
        Args:
            directory_path (str): path of directory. 
            poscar_name (str): Name of POSCAR file. Defaults to POSCAR.
            outcar_name (str): Name of OUTCAR file. Defaults to OUTCAR.
            vasprun_name (str): Name of vasprun.xml file.
                                Defaults to vasprun.xml.
        """
        path_poscar = directory_path + poscar_name
        path_outcar = directory_path + outcar_name
        path_vasprun = directory_path + vasprun_name
        poscar = Poscar.from_file(path_poscar)
        outcar = Outcar(path_outcar)
        vasprun = Vasprun(path_vasprun)
        self._structure = poscar.structure
        self._energy = outcar.final_energy
        self._dielectric_tensor \
            = np.array(outcar.dielectric_tensor) + \
            np.array(outcar.dielectric_ionic_tensor)
        self._eigen_value = vasprun.eigenvalues
        self._ewald_param = None

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

    def generate_neighbor_lattices(self,
                                   max_length,
                                   include_self,
                                   is_reciprocal = False):
        """
        Generator of a set of lattice vectors within the max length.
        Note that angles between any two axes are assumed to be between 60 and
        120 deg.
        Args:
            max_length (float): Max length to search lattice set.
            include_self (bool): Flag if (0, 0, 0) will be yield.
            is_reciprocal (bool): If true, generate reciprocal lattice set.
                                  Default is False.
        Yields (numpy.ndarray): Cartesian vector of lattice point.
        """
        if is_reciprocal:
            lattice_vectors = self._structure.lattice.reciprocal_lattice.matrix
        else:
            lattice_vectors = self._structure.lattice.matrix
        max_int = [int(max_length / np.linalg.norm(lattice_vectors[i])) + 1
                   for i in range(3)]
        for index in product(range(-max_int[0], max_int[0] + 1),
                             range(-max_int[1], max_int[1] + 1),
                             range(-max_int[2], max_int[2] + 1)):
            if (not include_self) and index == (0, 0, 0):
                continue
            cart_vector = np.dot(lattice_vectors.transpose(), np.array(index))
            if np.linalg.norm(cart_vector) < max_length:
                yield cart_vector

    def get_ewald_param(self,
                        computes_again = False,
                        initial_value = None,
                        convergence = 1.05,
                        prod_cutoff_fwhm = 25.0):
        """
        Get optimized ewald parameter.
        Once optimized parameter is calculated (usually slow),
        the value is cached and from next calculation,
        return the cached value if you don't specify
        computes_again = True.
        Args:
            initial_value (float): Initial guess of parameter.
            computes_again (bool): If you want to re-calculate parameter
             (e.g. change initial value), make this flat true.
            convergence (float):
                If 1/convergence < n_(real)/n_(reciprocal) < convergence,
                where n_(real) and n_(reciprocal) is number of real lattices
                and reciprocal lattices, finishes optimization and
                returns ewald_param.
            prod_cutoff_fwhm (float):
                product of cutoff radius of G-vector and gaussian FWHM.
        Returns (float):
            Optimized ewald_param.
        """
        root_det_dielectric = math.sqrt(np.linalg.det(self.dielectric_tensor))
        if self._ewald_param and not computes_again:
            return self._ewald_param
        real_lattice = self._structure.lattice.matrix
        reciprocal_lattice = self._structure.lattice.reciprocal_lattice.matrix
        cube_root_vol = math.pow(self._structure.lattice.volume, 1 / 3)
        if initial_value is not None:
            ewald_param = initial_value
        else:
            # determine initial ewald parameter to satisfy following:
            # max_int(Real) = max_int(Reciprocal)
            # in generate_neighbor_lattices function.
            # Left term:
            # max_int(Real) = 2 * x * Y  / l_r where x, Y, and l_r are ewald,
            # prod_cutoff_fwhm, and axis length of real lattice, respectively.
            # Right term:
            # max_int(reciprocal) = Y  / (x * l_g)
            # where l_g is axis length of reciprocal lattice, respectively.
            # Then, x = sqrt(l_g / l_r / 2)
            # gmean : geometric mean, like (a1 * a2 * a3)^(1/3)
            l_r = mstats.gmean([np.linalg.norm(v) for v in real_lattice])
            l_g = mstats.gmean([np.linalg.norm(v) for v in reciprocal_lattice])
            ewald_param \
                = np.sqrt(l_g / l_r / 2) *\
                cube_root_vol / root_det_dielectric
        while True:
            ewald = ewald_param / cube_root_vol * root_det_dielectric
            # count number of real lattice
            num_real_lattice = 0
            max_r_vector_norm = prod_cutoff_fwhm / ewald
            for _ in self.generate_neighbor_lattices(max_r_vector_norm,
                                                     include_self = True,
                                                     is_reciprocal = False):
                num_real_lattice += 1
            # count number of reciprocal lattice
            num_reciprocal_lattice = 0
            max_g_vector_norm = 2 * ewald * prod_cutoff_fwhm
            for _ in self.generate_neighbor_lattices(max_g_vector_norm,
                                                     include_self = False,
                                                     is_reciprocal = True):
                num_reciprocal_lattice += 1
            diff_real_reciprocal = num_real_lattice / num_reciprocal_lattice
            if 1 / convergence < diff_real_reciprocal < convergence:
                return ewald
            else:
                ewald_param *= diff_real_reciprocal ** 0.17


