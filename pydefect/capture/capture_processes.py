""" Module to identify the possible carrier capture processes which can be
considered from a given set of defect data
"""

import json
import yaml
import itertools

from pprint import pprint
from pathlib import Path
from typing import List
from pydefect.analysis.defect import Defect


class CaptureProcesses:

    def __init__(self, defect_instances, defect_dir=None):
        self.defect_instances = defect_instances
        self.capture_processes = {}
        self.defect_dir = defect_dir
        self.find_processes(defect_instances)

    @classmethod
    def from_directory(cls, defect_dir: str):
        """Identifies candidate carrier capture processes from a set of
        calculations in a given directory.

        Parameters:
            defect_dir (str): Search directory.

        Returns:

        """

        path_list = Path(defect_dir).rglob('defect.json')

        defect_instances = []
        for path in path_list:
            with path.open() as file:
                data = json.load(file)
                if data['@class'] == "Defect" and data['is_converged'] is True:
                    defect_instances.append(
                        Defect.load_json(str(path.resolve())))

        return cls(defect_instances, defect_dir=defect_dir)

    def find_processes(self, defect_instances: List):
        """Identifies candidate carrier capture processes from a list of
        Defect instances"""

        for d1, d2 in itertools.combinations(defect_instances, 2):
            if self.pass_checks(d1, d2):
                self.capture_processes[
                    d1.name + "_" + str(d1.charge) + "_" + str(d2.charge)] = \
                    {'name': d1.name, 'charges': [d1.charge, d2.charge]}

    @staticmethod
    def pass_checks(d1, d2):

        def localized_state(Defect):
            string_list = [v.value for v in Defect.band_edge_states.values()]
            if "Localized state" in string_list:
                return True

            return None

        if d1.name == d2.name:
            if abs(d1.charge - d2.charge) == 1:
                if localized_state(d1) or localized_state(d2):
                    return True
        else:
            return False

    def to_yaml_file(self):
        """Prints class attributes to a YAML file"""
        with open("capture_processes.yaml", "w") as f:
            f.write(yaml.dump(self.as_dict()))

    def print_to_terminal(self):
        """Prints class attributes to terminal screen"""
        pprint(self.as_dict())

    def as_dict(self) -> dict:

        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__,
             "capture_processes": self.capture_processes,
             "defect_dir": self.defect_dir}
        return d

    # TODO: option to calc delta Q and delta E (use displacements)
    # TODO: move pass checks all in one place
    # TODO: option to do unconverged
    # TODO: option to include non-local states
    # TODO: error handling and tests
