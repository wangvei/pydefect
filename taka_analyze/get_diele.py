#!/bin/env python3                                                                                                                       
import sys
from pymatgen.io.vasp.outputs import Outcar
import json
from monty.json import MontyEncoder

outcar_path=sys.argv[1]
o=Outcar(outcar_path)

fw = open("correction.json", 'w')

a = {}
a["dielectric_tensor"] = o.dielectric_tensor
a["dielectric_ionic_tensor"] = o.dielectric_ionic_tensor

json.dump(a, fw, indent=2, cls=MontyEncoder)
#print(o.dielectric_tensor)
