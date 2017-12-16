import ../defect_input

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2017, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "December 4, 2017"

a = defect_input.DefectIn.from_str_file(dopants=["Sb"], interstitials="0 0 0.125", symbreak=False, is_antisite=True, ElNeg_diff=3, irregular="Ag_Sr1")

a.to()
