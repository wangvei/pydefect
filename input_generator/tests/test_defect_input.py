from pydefect.input_generator.defect_input import *

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "December 4, 2017"

a = DefectIn.from_str_file(poscar="POSCAR-ScN64atoms-test", dopants=["Sb"], interstitials="0 0 0.125", symbreak=False, symprec=1e-05, is_antisite=True, ElNeg_diff=3, include="Ag_Sr1", exclude="")

#a = DefectIn.from_str_file(poscar="POSCAR-ScN64atoms-test")
#a.to()
print(a.setting.as_dict())

#b = DefectIn.from_)defect_in()
