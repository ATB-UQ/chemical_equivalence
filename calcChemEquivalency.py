import logging
from optparse import OptionParser

from chemical_equivalence.log_helpers import print_stderr
from chemical_equivalence.NautyInterface import NautyInterface
from chemical_equivalence.molData import MolData, MolDataFailure
from chemical_equivalence.chiral import containsStereoheterotopicAtoms
from chemical_equivalence.doubleBonds import containsEquivalenceBreakingDoubleBond
from chemical_equivalence.rings import containsInversableRings

def getChemEquivGroups(molData, log=None):
    # for cases with stereogenic centers we need to add flavour (some additional degree of freedom)
    # to distinguish stereoheterotopic atoms
    class FlavourCounter(object):
        def __init__(self):
            self._i = 0
        def getNext(self):
            self._i += 1
            return self._i 

    nautyInterface = NautyInterface(molData)

    nautyInterface.calcEquivGroups(log)

    flavourCounter = FlavourCounter()
    if chemicalEquivalenceExceptions(molData, flavourCounter, log):
        clearEqGroupData(molData)
        nautyInterface.calcEquivGroups(log)

    return molData.equivalenceGroups

def chemicalEquivalenceExceptions(molData, flavourCounter, log):

    exceptionSearchingFunctions = [containsStereoheterotopicAtoms, containsEquivalenceBreakingDoubleBond, containsInversableRings]

    # If there is a chemical equivalence breaking groups then should_rerun = True
    should_rerun = any([func(molData, flavourCounter, log) for func in exceptionSearchingFunctions])

    if log:
        if should_rerun: log.debug("Molecule contains chemical equivalence breaking groups.")
        else:            log.debug("Molecule has NO chemical equivalence breaking groups.")

    return should_rerun

def clearEqGroupData(molData):
    molData.equivalenceGroups = {}
    for atom in list(molData.atoms.values()):
        del atom["equivalenceGroup"]

def partial_mol_data_for_pdbstr(pdb_string, united_atoms=True, debug=False):
    assert pdb_string, 'Empty PDB string'

    data = MolData(pdb_string)
    getChemEquivGroups(data)
    if united_atoms:
        if debug:
            print_stderr(
                "All atoms: {0}\n".format(
                    "".join([a["type"] for a in list(data.atoms.values())]),
                ),
            )
        data.unite_atoms()
        if debug:
            print_stderr(
                "United atoms: {0}\n".format(
                    "".join([a["type"] for a in list(data.atoms.values()) if "uindex" in a]),
                ),
            )
    return data

def parseCommandline():
    # run CGP and symmetrization from commandline arguments

    parser = OptionParser(usage="python %prog options filename [options [filename]]")
    # required
    parser.add_option("-p", "--pdb", dest="pdb", help="PDB structure file", action="store", type="string")

    (opts, _) = parser.parse_args()

    if not opts.pdb or not opts.mtb:
        parser.print_help()
    else:
        # Try and open pdb file
        try:
            fh = open(opts.pdb, "r")
            pdbStr = fh.read()
            fh.close()

        except:
            print_stderr("ERROR: there was a problem with the PDB file.")
            parser.print_help()
            return

    data = MolData(pdbStr)

    nautyInterface = NautyInterface(data)

    # Run symmetrization
    nautyInterface.calcSym()
    print(nautyInterface.data.equivalenceGroups)

if __name__=="__main__":
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - [%(levelname)s] - %(message)s  -->  (%(module)s.%(funcName)s: %(lineno)d)', datefmt='%d-%m-%Y %H:%M:%S')
    #parseCommandline()
    #data = MolData(open("testing/stressTest.pdb").read(), open("testing/stressTest.mtb").read())
    #data = MolData(open("testing/benzene.pdb").read())
    #data = MolData(open("testing/cyclohexane.pdb").read())
    data = MolData(open("testing/decalin.pdb").read(), logging.getLogger())
    #data = MolData(open("testing/1,1-dichlorocyclohexane.pdb").read(), open("testing/1,1-dichlorocyclohexane.mtb").read())
    #data = MolData(open("testing/cyclohexane.pdb").read(), open("testing/cyclohexane.mtb").read())
    #data = MolData(open("testing/chlorocyclohexane.pdb").read(), open("testing/chlorocyclohexane.mtb").read())
    #data = MolData(open("testing/1,1-dichlorocyclohexane.pdb").read(), open("testing/1,1-dichlorocyclohexane.mtb").read())
    #data = MolData(open("testing/cyclohexane.pdb").read(), open("testing/cyclohexane.mtb").read())
    #data = MolData(open("testing/chlorocyclohexane.pdb").read(), open("testing/chlorocyclohexane.mtb").read())
    #data = MolData(open("testing/decalin.pdb").read(), open("testing/decalin.mtb").read())
    #data = MolData(open("testing/cyclobutadiene.pdb").read(), open("testing/cyclobutadiene.mtb").read())
    #data = MolData(open("testing/(1S,4S)-1,4-dibromo-1,4-dichloro-2,2,3,3-tetramethylbutane.pdb").read(), open("testing/(1S,4S)-1,4-dibromo-1,4-dichloro-2,2,3,3-tetramethylbutane.mtb").read())

    #data = MolData(open("testing/dichloro-dibromobutane.pdb").read(), open("testing/dichloro-dibromobutane.mtb").read())
    #data = MolData(open("testing/(1S,3S)-1,3-dibromo-1,3-dichloropropane.pdb").read(), open("testing/(1S,3S)-1,3-dibromo-1,3-dichloropropane.mtb").read())
    #data = MolData(open("testing/(1R,3S)-1,3-dibromo-1,3-dichloropropane.pdb").read(), open("testing/(1R,3S)-1,3-dibromo-1,3-dichloropropane.mtb").read())

    print("Rings: {0}".format( data.rings ))
    getChemEquivGroups(data, log=logging.getLogger())
    print(list(data.equivalenceGroups.values()))
    print("\n".join([str((at["index"], at["equivalenceGroup"])) for at in list(data.atoms.values())]))
