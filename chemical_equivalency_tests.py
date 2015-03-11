import unittest
from molData import MolData
import calcChemEquivalency
import logging

class ChemicalEquivalencyTest(unittest.TestCase):
    
    def run_unit_check(self, base_file, expected_result_list):
        data = MolData(open("testing/{base_file}.pdb".format(base_file=base_file)).read(), open("testing/{base_file}.mtb".format(base_file=base_file)).read())
        #logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - [%(levelname)s] - %(message)s  -->  (%(module)s.%(funcName)s: %(lineno)d)', datefmt='%d-%m-%Y %H:%M:%S')
        #calcChemEquivalency.getChemEquivGroups(data, log=logging.getLogger())
        calcChemEquivalency.getChemEquivGroups(data)
        result_set = set(map(frozenset, data.equivalenceGroups.values()))
        expectedResultSet = set(map(frozenset, expected_result_list))
        
        self.assertEqual(result_set, expectedResultSet)
        

    def testGlucose(self):
        self.run_unit_check("glucose_AA", [])
        

    def testTrueChiral(self):
        self.run_unit_check("trueChiral", [[7, 8, 9]])
        
    
    def testCNT(self):
        self.run_unit_check("CNT", [[1, 18, 20, 24, 26, 30, 32, 36, 38, 42, 44, 48, 50, 54, 56, 217, 219, 223, 225, 229, 231, 235, 237, 241, 243, 247, 249, 252], [2, 17, 19, 23, 25, 29, 31, 35, 37, 41, 43, 47, 49, 53, 55, 216, 218, 222, 224, 228, 230, 234, 236, 240, 242, 246, 248, 251], [3, 16, 21, 22, 27, 28, 33, 34, 39, 40, 45, 46, 51, 52, 57, 214, 215, 220, 221, 226, 227, 232, 233, 238, 239, 244, 245, 250], [4, 15, 58, 59, 62, 63, 66, 67, 70, 71, 74, 75, 78, 79, 82, 188, 189, 192, 193, 196, 197, 200, 201, 204, 205, 208, 209, 213], [5, 14, 60, 61, 64, 65, 68, 69, 72, 73, 76, 77, 80, 81, 108, 187, 190, 191, 194, 195, 198, 199, 202, 203, 206, 207, 210, 211], [6, 13, 83, 86, 87, 90, 91, 94, 95, 98, 99, 102, 103, 106, 107, 164, 165, 168, 169, 172, 173, 176, 177, 180, 181, 184, 185, 212], [7, 12, 84, 85, 88, 89, 92, 93, 96, 97, 100, 101, 104, 105, 109, 162, 163, 166, 167, 170, 171, 174, 175, 178, 179, 182, 183, 186], [8, 11, 110, 111, 114, 115, 118, 119, 122, 123, 126, 127, 130, 131, 134, 136, 137, 140, 141, 144, 145, 148, 149, 152, 153, 156, 157, 161], [9, 10, 112, 113, 116, 117, 120, 121, 124, 125, 128, 129, 132, 133, 135, 138, 139, 142, 143, 146, 147, 150, 151, 154, 155, 158, 159, 160]])
        
    
    def testButadiene(self):
        self.run_unit_check("butadiene", [[1, 3], [9, 10], [2, 8], [4, 6], [5, 7]])
        

    def testBenzene(self):
        self.run_unit_check("benzene", [[1, 4, 6, 8, 10, 12], [2, 3, 5, 7, 9, 11]])
    
    def testChlorobenzene(self):
        self.run_unit_check("chlorobenzene", [[3, 5], [4, 6], [7, 11], [8, 12]])
    
    def testStressTest(self):
        self.run_unit_check("stressTest", [[1, 3, 5, 7, 9, 11], [2, 4, 6, 8, 10, 12]])
        
    def testChloroethene(self):
        self.run_unit_check("chloroethene", [])
    
    def testChloroBromopropane(self):
        self.run_unit_check("1-chloro-1-bromopropane", [[9, 10, 11]])
    
    def testChlorocyclohexane(self):
        self.run_unit_check("chlorocyclohexane", [[4, 7], [10, 16], [11, 17], [12, 18], [5, 9], [6, 8]])
        
    def testCyclohexane(self):
        self.run_unit_check("cyclohexane", [[1, 3, 5, 6, 8, 9, 11, 12, 14, 15, 17, 18], [2, 4, 7, 10, 13, 16]])    
    
suite = unittest.TestLoader().loadTestsFromTestCase(ChemicalEquivalencyTest)
unittest.TextTestRunner(verbosity=4).run(suite)
