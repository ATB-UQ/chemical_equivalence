import unittest
from molData import MolData
import calcChemEquivalency

class ChemicalEquivalencyTest(unittest.TestCase):

    def testGlucose(self):
        base_file = "glucose_AA"
        data = MolData(open("testing/{base_file}.pdb".format(base_file=base_file)).read(), open("testing/{base_file}.mtb".format(base_file=base_file)).read())
        calcChemEquivalency.getChemEquivGroups(data)
        result_dict = dict( [ (at["index"], at["equivalenceGroup"]) for at in data.atoms.values()] )
        expected_result_dict = {1: -1, 2: -1, 3: -1, 4: -1, 5: -1, 6: -1, 7: -1, 8: -1, 9: -1, 10: -1, 11: -1, 12: -1, 13: -1, 14: -1, 15: -1, 16: -1, 17: -1, 18: -1, 19: -1, 20: -1, 21: -1, 22: -1, 23: -1, 24: -1}
        self.assertEqual(result_dict, expected_result_dict)

    def testTrueChiral(self):
        base_file = "trueChiral"
        data = MolData(open("testing/{base_file}.pdb".format(base_file=base_file)).read(), open("testing/{base_file}.mtb".format(base_file=base_file)).read())
        calcChemEquivalency.getChemEquivGroups(data)
        result_dict = dict( [ (at["index"], at["equivalenceGroup"]) for at in data.atoms.values()] )
        expected_result_dict = {1: -1, 2: -1, 3: -1, 4: -1, 5: -1, 6: -1, 7: 0, 8: 0, 9: 0}
        self.assertEqual(result_dict, expected_result_dict)


suite = unittest.TestLoader().loadTestsFromTestCase(ChemicalEquivalencyTest)
unittest.TextTestRunner(verbosity=2).run(suite)
