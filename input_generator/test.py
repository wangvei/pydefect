import defect_input

a = defect_input.DefectIn.from_str_file(dopants=["Sb"], interstitials="0 0 0.125", symbreak=False, is_antisite=True, ElNeg_diff=3, irregular="Ag_Sr1")

a.to()
