from NautyInterface import NautyInterface
from molData import MolData
from optparse import OptionParser
from chiral import containsDiastereotopicAtoms
from doubleBonds import containsEquivalenceBreakingDoubleBond
import logging
symMaxDiff = 0.2 #Maximum absolute charge difference allowed for atoms
    

def getChemEquivGroups(molData, log=None):
    # for cases with stereogenic centers we need to add flavour (some additional degree of freedom)
    # to distinguish diastereotopic atoms
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
    
    exceptionSearchingFunctions = [containsDiastereotopicAtoms, containsEquivalenceBreakingDoubleBond]
    
    # If there is a chemical equivalence breaking groups then should_rerun = True
    should_rerun = any([func(molData, flavourCounter, log) for func in exceptionSearchingFunctions])
    
    if log:
        if should_rerun: log.info("Molecule contains chemical equivalence breaking groups.")
        else:            log.info("Molecule has NO chemical equivalence breaking groups.")
        
    return should_rerun

def clearEqGroupData(molData):
    molData.equivalenceGroups = {}
    for atom in molData.atoms.values():
        del atom["equivalenceGroup"]
    
def parseCommandline():
    # run CGP and symmetrization from commandline arguments
    
    parser = OptionParser(usage="python %prog options filename [options [filename]]")
    # required
    parser.add_option("-p", "--pdb", dest="pdb", help="PDB structure file", action="store", type="string")
    parser.add_option("-m", "--mtb", dest="mtb", help="MTB file", action="store", type="string")
    
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
            print "ERROR: there was a problem with the PDB file."
            parser.print_help()
            return
            
        # Try and open mtb file
        try:
            fh = open(opts.mtb, "r")
            mtbStr = fh.read()
            fh.close()
            
        except:
            print "ERROR: there was a problem with the MTB file."
            parser.print_help()
            return
        
    data = MolData(pdbStr, mtbStr)
    
    nautyInterface = NautyInterface(data)
    
    # Run symmetrization
    nautyInterface.calcSym()
    print nautyInterface.data.equivalenceGroups
    

if __name__=="__main__":
    #parseCommandline()
    data = MolData(open("testing/pseudoChiral.pdb").read(), open("testing/pseudoChiral.mtb").read())
    #data = MolData(open("testing/butadiene.pdb").read(), open("testing/butadiene.mtb").read())
    #data = MolData(open("testing/glucose_AA.pdb").read(), open("testing/glucose_AA.dat").read())
    #data = MolData(open("testing/trueChiral.pdb").read(), open("testing/trueChiral.mtb").read())
    #data = MolData(open("testing/1-chloro-1-bromopropane.pdb").read(), open("testing/1-chloro-1-bromopropane.mtb").read())
    #data = MolData(open("testing/(1S,4S)-1,4-dibromo-1,4-dichloro-2,2,3,3-tetramethylbutane.pdb").read(), open("testing/(1S,4S)-1,4-dibromo-1,4-dichloro-2,2,3,3-tetramethylbutane.mtb").read())
    #data = MolData(open("testing/CNT.pdb").read(), open("testing/CNT.mtb").read())
    #data = MolData(open("testing/stressTest.pdb").read(), open("testing/stressTest.mtb").read())
    #data = MolData(open("testing/dichloro-dibromobutane.pdb").read(), open("testing/dichloro-dibromobutane.mtb").read())
    #data = MolData(open("testing/(1S,3S)-1,3-dibromo-1,3-dichloropropane.pdb").read(), open("testing/(1S,3S)-1,3-dibromo-1,3-dichloropropane.mtb").read())
    #data = MolData(open("testing/(1R,3S)-1,3-dibromo-1,3-dichloropropane.pdb").read(), open("testing/(1R,3S)-1,3-dibromo-1,3-dichloropropane.mtb").read())
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - [%(levelname)s] - %(message)s  -->  (%(module)s.%(funcName)s: %(lineno)d)', datefmt='%d-%m-%Y %H:%M:%S')
    getChemEquivGroups(data, log=logging.getLogger())
    print "\n".join([str((at["index"], at["equivalenceGroup"])) for at in data.atoms.values()])
