from NautyInterface import NautyInterface
from molData import MolData
from optparse import OptionParser
from chiral import wrongChemicalEquivalencies
symMaxDiff = 0.2 #Maximum absolute charge difference allowed for atoms
    

def getChemEquivGroups(molData, log=None):
    
    nautyInterface = NautyInterface(molData)
    
    nautyInterface.calcEquivGroups(log)
    
    while wrongChemicalEquivalencies(molData):
        clearEqGroupData(molData)
        nautyInterface.calcEquivGroups(log)
    
    # for highly symmetric molecules we don't want any diastereotopic atoms
    # in any symmetry groups, so we'll remove them manually
    atomsWithFlavour = []
    for atomID, atom in molData.atoms.items():
        if atom.has_key("flavour"):
            atom["symGr"] = -1
            atomsWithFlavour.append(atomID)
    
    for eqGrp in molData.symgroups.values():
        matched = [atm for atm in atomsWithFlavour if atm in eqGrp] 
        if matched:
            for matchedAtom in matched:
                eqGrp.remove(matchedAtom)
    
    # check that none of the symmetry groups are empty
    molData.symgroups = dict([(k,v) for k,v in molData.symgroups.items() if v])
    
    
    print molData.symgroups
    print "\n".join([str((at["index"], at["symGr"])) for at in molData.atoms.values()])
    
def clearEqGroupData(molData):
    molData.symgroups = {}
    for atom in molData.atoms.values():
        del atom["symGr"]
    
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
    print nautyInterface.data.symgroups
    

if __name__=="__main__":
    #parseCommandline()
    #data = MolData(open("testing/pseudoChiral.pdb").read(), open("testing/pseudoChiral.mtb").read())
    #data = MolData(open("testing/glucose.pdb").read(), open("testing/glucose.dat").read())
    #data = MolData(open("testing/trueChiral.pdb").read(), open("testing/trueChiral.mtb").read())
    #data = MolData(open("testing/1-chloro-1-bromopropane.pdb").read(), open("testing/1-chloro-1-bromopropane.mtb").read())
    data = MolData(open("testing/(1S,4S)-1,4-dibromo-1,4-dichloro-2,2,3,3-tetramethylbutane.pdb").read(), open("testing/(1S,4S)-1,4-dibromo-1,4-dichloro-2,2,3,3-tetramethylbutane.mtb").read())
    getChemEquivGroups(data)