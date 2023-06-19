import sys
from os.path import join
import unittest

from chemical_equivalence.calcChemEquivalency import partial_mol_data_for_pdbstr
from chemical_equivalence.test import TESTING_DIR

class ChemicalEquivalencyTest(unittest.TestCase):

    def run_unit_check(self, base_file):
        partial_mol_data_for_pdbstr(
            open(
                join(
                    TESTING_DIR,
                    "{base_file}.pdb".format(base_file=base_file),
                ),
            ).read(),
            united_atoms=True,
            debug=True,
        )

    def test_benzene(self):
        self.run_unit_check("benzene")

    def test_butadiene(self):
        self.run_unit_check("butadiene")

suite = unittest.TestLoader().loadTestsFromTestCase(ChemicalEquivalencyTest)
unittest.TextTestRunner(verbosity=4).run(suite)
