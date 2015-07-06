import sys
from os.path import join, abspath
sys.path.append(join(abspath(__file__), "../"))
import calcChemEquivalency
import unittest

class ChemicalEquivalencyTest(unittest.TestCase):

    def run_unit_check(self, base_file):
        calcChemEquivalency.partial_mol_data_for_pdbstr(open("{base_file}.pdb".format(base_file=base_file)).read(), united_atoms=True, debug=True)

    def test_benzene(self):
        self.run_unit_check("benzene")

    def test_butadiene(self):
        self.run_unit_check("butadiene")

suite = unittest.TestLoader().loadTestsFromTestCase(ChemicalEquivalencyTest)
unittest.TextTestRunner(verbosity=4).run(suite)