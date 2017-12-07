#!/usr/bin/env python

import sys
from defect_property import DefectProperty

dir_name = sys.argv[1]
dp = DefectProperty.from_directory(dir_name)
d = dp.as_dict()
print(d["energy"])
print(d["structure"])
print(d["atomic_site_pot"])
print(d)
print(dp)
print(dp.to())
dp.to(filename="hoge.dat")
dp2 = DefectProperty.from_file("hoge.dat")
