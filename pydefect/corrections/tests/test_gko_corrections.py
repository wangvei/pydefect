from pydefect.util.testing import PydefectTest
from pymatgen.core.structure import Structure
from pydefect.corrections.gko_corrections import GkoCorrection
from pydefect.corrections.efnv_corrections import ExtendedFnvCorrection
from pydefect.core.unitcell_calc_results import UnitcellCalcResults

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class GkoCorrectionTest(PydefectTest):

    def test(self):
        file = (PydefectTest.TEST_FILES_DIR
                / "core" / "MgO" / "defects" / "Va_O1_2" / "CONTCAR")
        structure = Structure.from_file(file)

        file = (PydefectTest.TEST_FILES_DIR
                / "core" / "MgO" / "unitcell" / "unitcell.json")
        unitcell = UnitcellCalcResults.load_json(file)

        file = (PydefectTest.TEST_FILES_DIR
                / "core" / "MgO" / "defects" / "Va_O1_2" / "correction.json")
        efnv = ExtendedFnvCorrection.load_json(file)

        self.mgo = GkoCorrection.from_scratch(structure=structure,
                                              before_charge=2,
                                              after_charge=1,
                                              unitcell_dft=unitcell,
                                              efnv_correction=efnv)

        print(self.mgo)




