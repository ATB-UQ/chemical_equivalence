from NautyInterface import NautyInterface
from molData import MolData
from optparse import OptionParser

symMaxDiff = 0.2 #Maximum absolute charge difference allowed for atoms
    

def getChemEquivGroups(molData, log=None):
    
    nautyInterface = NautyInterface(molData)
    nautyInterface.calcEquivGroups(log)
    
    
    
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
    data = MolData(open("testing/manyEqGrps.pdb").read(), open("testing/manyEqGrps.mtb").read())
    #parseCommandline()
    nautyInterface = NautyInterface(data)
    nautyInterface.calcEquivGroups()
    print nautyInterface.data.symgroups
    #print "\n".join([str((at["index"], at["symGr"])) for at in nautyInterface.data.atoms.values()])
     

    
    