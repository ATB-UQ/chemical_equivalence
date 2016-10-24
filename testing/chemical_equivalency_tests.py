import sys
from os.path import join, abspath
import unittest
import logging
from os import devnull
devnull_stream = open(devnull, 'w')

from chemical_equivalence.molData import MolData
from chemical_equivalence.calcChemEquivalency import getChemEquivGroups
from chemical_equivalence.test import TESTING_DIR

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - [%(levelname)s] - %(message)s  -->  (%(module)s.%(funcName)s: %(lineno)d)', datefmt='%d-%m-%Y %H:%M:%S')
log = logging.getLogger()
log = None

class ChemicalEquivalencyTest(unittest.TestCase):

    def run_unit_check(self, base_file, expected_result_list):
        data = MolData(
            open(
                join(
                    TESTING_DIR,
                    "{base_file}.pdb".format(base_file=base_file),
                ),
            ).read(),
            log,
        )
        results = getChemEquivGroups(data, log=log)

        self.assertEqual(
            results,
            expected_result_list,
            msg='Expected: {0}; GOT: {1}'.format(expected_result_list, results),
        )

    def test_glucose(self):
        self.run_unit_check(
            "glucose_AA",
            {1: 0, 2: 1, 3: 2, 4: 3, 5: 4, 6: 5, 7: 6, 8: 7, 9: 8, 10: 9, 11: 10, 12: 11, 13: 12, 14: 13, 15: 14, 16: 15, 17: 16, 18: 17, 19: 18, 20: 19, 21: 20, 22: 21, 23: 22, 24: 23},
        )

    def test_truechiral(self):
        self.run_unit_check(
            "trueChiral",
            {1: 0, 2: 1, 3: 2, 4: 3, 5: 4, 6: 5, 7: 6, 8: 6, 9: 6},
        )

    def test_cnt(self):
        self.run_unit_check(
            "CNT",
            {1: 0, 2: 1, 3: 2, 4: 3, 5: 4, 6: 5, 7: 6, 8: 7, 9: 8, 10: 8, 11: 7, 12: 6, 13: 5, 14: 4, 15: 3, 16: 2, 17: 1, 18: 0, 19: 1, 20: 0, 21: 2, 22: 2, 23: 1, 24: 0, 25: 1, 26: 0, 27: 2, 28: 2, 29: 1, 30: 0, 31: 1, 32: 0, 33: 2, 34: 2, 35: 1, 36: 0, 37: 1, 38: 0, 39: 2, 40: 2, 41: 1, 42: 0, 43: 1, 44: 0, 45: 2, 46: 2, 47: 1, 48: 0, 49: 1, 50: 0, 51: 2, 52: 2, 53: 1, 54: 0, 55: 1, 56: 0, 57: 2, 58: 3, 59: 3, 60: 4, 61: 4, 62: 3, 63: 3, 64: 4, 65: 4, 66: 3, 67: 3, 68: 4, 69: 4, 70: 3, 71: 3, 72: 4, 73: 4, 74: 3, 75: 3, 76: 4, 77: 4, 78: 3, 79: 3, 80: 4, 81: 4, 82: 3, 83: 5, 84: 6, 85: 6, 86: 5, 87: 5, 88: 6, 89: 6, 90: 5, 91: 5, 92: 6, 93: 6, 94: 5, 95: 5, 96: 6, 97: 6, 98: 5, 99: 5, 100: 6, 101: 6, 102: 5, 103: 5, 104: 6, 105: 6, 106: 5, 107: 5, 108: 4, 109: 6, 110: 7, 111: 7, 112: 8, 113: 8, 114: 7, 115: 7, 116: 8, 117: 8, 118: 7, 119: 7, 120: 8, 121: 8, 122: 7, 123: 7, 124: 8, 125: 8, 126: 7, 127: 7, 128: 8, 129: 8, 130: 7, 131: 7, 132: 8, 133: 8, 134: 7, 135: 8, 136: 7, 137: 7, 138: 8, 139: 8, 140: 7, 141: 7, 142: 8, 143: 8, 144: 7, 145: 7, 146: 8, 147: 8, 148: 7, 149: 7, 150: 8, 151: 8, 152: 7, 153: 7, 154: 8, 155: 8, 156: 7, 157: 7, 158: 8, 159: 8, 160: 8, 161: 7, 162: 6, 163: 6, 164: 5, 165: 5, 166: 6, 167: 6, 168: 5, 169: 5, 170: 6, 171: 6, 172: 5, 173: 5, 174: 6, 175: 6, 176: 5, 177: 5, 178: 6, 179: 6, 180: 5, 181: 5, 182: 6, 183: 6, 184: 5, 185: 5, 186: 6, 187: 4, 188: 3, 189: 3, 190: 4, 191: 4, 192: 3, 193: 3, 194: 4, 195: 4, 196: 3, 197: 3, 198: 4, 199: 4, 200: 3, 201: 3, 202: 4, 203: 4, 204: 3, 205: 3, 206: 4, 207: 4, 208: 3, 209: 3, 210: 4, 211: 4, 212: 5, 213: 3, 214: 2, 215: 2, 216: 1, 217: 0, 218: 1, 219: 0, 220: 2, 221: 2, 222: 1, 223: 0, 224: 1, 225: 0, 226: 2, 227: 2, 228: 1, 229: 0, 230: 1, 231: 0, 232: 2, 233: 2, 234: 1, 235: 0, 236: 1, 237: 0, 238: 2, 239: 2, 240: 1, 241: 0, 242: 1, 243: 0, 244: 2, 245: 2, 246: 1, 247: 0, 248: 1, 249: 0, 250: 2, 251: 1, 252: 0},
        )

    def test_butadiene(self):
        self.run_unit_check(
            "butadiene",
            {1: 0, 2: 1, 3: 2, 4: 3, 5: 4, 6: 5, 7: 6, 8: 7, 9: 8, 10: 9},
        )

    def test_benzene(self):
        self.run_unit_check(
            "benzene",
            {1: 0, 2: 1, 3: 1, 4: 0, 5: 1, 6: 0, 7: 1, 8: 0, 9: 1, 10: 0, 11: 1, 12: 0},
        )

    def test_chlorobenzene(self):
        self.run_unit_check(
            "chlorobenzene",
            {1: 0, 2: 1, 3: 2, 4: 3, 5: 2, 6: 3, 7: 4, 8: 5, 9: 6, 10: 7, 11: 4, 12: 5},
        )

    def test_stresstest(self):
        self.run_unit_check(
            "stressTest",
            {1: 0, 2: 1, 3: 0, 4: 1, 5: 0, 6: 1, 7: 0, 8: 1, 9: 0, 10: 1, 11: 0, 12: 1},
        )

    def test_chloroethene(self):
        self.run_unit_check(
            "chloroethene",
            {1: 0, 2: 1, 3: 2, 4: 3, 5: 4, 6: 5},
        )

    def test_1_chloro_1_bromopropane(self):
        self.run_unit_check(
            "1-chloro-1-bromopropane",
            {1: 0, 2: 1, 3: 2, 4: 3, 5: 4, 6: 5, 7: 6, 8: 7, 9: 8, 10: 8, 11: 8},
        )

    def test_chlorocyclohexane(self):
        self.run_unit_check(
            "chlorocyclohexane",
            {1: 0, 2: 1, 3: 2, 4: 3, 5: 4, 6: 5, 7: 6, 8: 7, 9: 8, 10: 9, 11: 10, 12: 11, 13: 12, 14: 13, 15: 14, 16: 15, 17: 16, 18: 17},
        )

    def test_cyclohexane(self):
        self.run_unit_check(
            "cyclohexane",
            {1: 0, 2: 1, 3: 0, 4: 1, 5: 0, 6: 0, 7: 1, 8: 0, 9: 0, 10: 1, 11: 0, 12: 0, 13: 1, 14: 0, 15: 0, 16: 1, 17: 0, 18: 0},
        )

    def test_dichlorocyclohexane(self):
        self.run_unit_check(
            "1,1-dichlorocyclohexane",
            {1: 0, 2: 1, 3: 0, 4: 2, 5: 3, 6: 3, 7: 2, 8: 3, 9: 3, 10: 4, 11: 5, 12: 5, 13: 6, 14: 7, 15: 7, 16: 4, 17: 5, 18: 5},
        )

    def test_decalin(self):
        self.run_unit_check(
            "decalin",
            {1: 0, 2: 1, 3: 2, 4: 3, 5: 4, 6: 5, 7: 6, 8: 7, 9: 8, 10: 9, 11: 10, 12: 11, 13: 12, 14: 13, 15: 14, 16: 15, 17: 16, 18: 17, 19: 18, 20: 19, 21: 20, 22: 21, 23: 22, 24: 23, 25: 24, 26: 25, 27: 26, 28: 27},
        )

suite = unittest.TestLoader().loadTestsFromTestCase(ChemicalEquivalencyTest)
unittest.TextTestRunner(verbosity=4).run(suite)
