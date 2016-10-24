import sys
from os.path import join
import unittest
import logging

from chemical_equivalence.molData import MolData
from chemical_equivalence.calcChemEquivalency import getChemEquivGroups
from chemical_equivalence.build_rings import build_rings
from chemical_equivalence.test import TESTING_DIR

class ChemicalEquivalencyTest(unittest.TestCase):

    def run_unit_check(self, base_file, expected_result_list):
        data = MolData(
            open("{base_file}.pdb".format(base_file=base_file)).read(),
            open("testing/{base_file}.mtb".format(base_file=base_file)).read()
        )
        logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - [%(levelname)s] - %(message)s  -->  (%(module)s.%(funcName)s: %(lineno)d)', datefmt='%d-%m-%Y %H:%M:%S')
        getChemEquivGroups(data, log=logging.getLogger())
        build_rings(data)
        result_set = set(map(frozenset, data.equivalenceGroups.values()))
        expectedResultSet = set(map(frozenset, expected_result_list))

        self.assertEqual(result_set, expectedResultSet)

    def test23391(self):
        self.run_unit_check("22196", [])

    def test16806(self):
        self.run_unit_check("22195", [])

suite = unittest.TestLoader().loadTestsFromTestCase(ChemicalEquivalencyTest)
unittest.TextTestRunner(verbosity=4).run(suite)
