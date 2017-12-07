#!/usr/bin/env python

import sys
from defect_property import DefectProperty

dir_name = sys.argv[1]
#dp_from_dir = DefectProperty.from_directory(dir_name)
#dp_from_dir.to(filename="dp_from_dir.json")

dp_file = DefectProperty.from_file("dp_from_dir.json")
dp_file.to(filename="dp_from_file.json")
print(dp_file.__energy)
sys.exit()

d = dp_from_dir.as_dict()
dp_dict = DefectProperty.from_dict(d)
dp_dict.to(filename="dp_from_dict.json")
s = dp_from_dir.to()
dp_str = DefectProperty.from_str(s)
dp_str.to(filename="dp_from_str.json")

