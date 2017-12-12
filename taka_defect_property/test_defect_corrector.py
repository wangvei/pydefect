#!/usr/bin/env python

import sys
from defect_property import DefectProperty
from defect_corrector import DefectCorrector

#defect_prop = DefectProperty.from_directory("/Users/takahashi/psi/bitbucket/BaBi2S4/Ba_Bi1_0/")
#perfect_prop = DefectProperty.from_directory("/Users/takahashi/psi/bitbucket/BaBi2S4/perfect/")
#defect_prop.to("defect_property.json")
#perfect_prop.to("perfect_property.json")
defect_prop = DefectProperty.from_file("./defect_property.json")
perfect_prop = DefectProperty.from_file("./perfect_property.json")
def_corr = DefectCorrector(defect_prop, perfect_prop)
def_corr.correct()
