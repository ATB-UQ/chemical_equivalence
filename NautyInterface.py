import subprocess
import sys
from optparse import OptionParser
import tempfile

PATH_TO_SYM = "/usr/local/bin/symmetrize"

class NautyInterface(object):

    symMaxDiff = 0.2 #Maximum absolute charge difference allowed for atoms
    
    def __init__(self, pdbStr=None, mtbStr=None, atbDataObj=None):
        if atbDataObj is not None:
            self.data = atbDataObj
        else:
            self.data = CGPData(pdbStr, mtbStr)
        
    
    def calcSym(self, log=None, method=None):
        stdinput, args = self._genSymInput(log, method)
        lgfString = _run(args,stdinput,errorLog=log)
        
        if len(lgfString) == 0:
            if log is not None:
                log.warning("calcSym: symmetrize produced no output; could be due to timeout")
            return ""
        # update charges after sym
        nodeBlock=False
        lgf_lines = lgfString.splitlines()
        i = 0
        while i < len(lgf_lines):
            line = lgf_lines[i]
            if line.startswith("#"):    continue
            if not nodeBlock and line.strip() == "@nodes":
                nodeBlock=True
                i += 2 # skip the label block
                continue
            if line.strip() == "@edges":
                nodeBlock=False
                break
            if nodeBlock:
                #charge, atomID, atmName, atomType, coordX, coordY, coordZ, initColor, symGrp
                (atmCh, atmID, _, _, _, _, _, _, symGrp) = line.split()[0:9]
                symGrp = int(symGrp)
                atmID = int(atmID)
                # update atomic charge after sym
                self.data[atmID]["charge"] = float(atmCh)
                
                # store symmetry data
                self.data[atmID]["symGr"] = int(symGrp)
                if symGrp != -1:
                    if symGrp not in self.data.symgroups.keys():
                        self.data.symgroups[symGrp] = [atmID]
                    else:
                        self.data.symgroups[symGrp].append(atmID)
                        
                i += 1
        return self.data.symgroups
        
        
    def _genSymInput(self, log=None, method=None):
        
        lgfString = self._writeLGF()
        args = [PATH_TO_SYM]
        args.extend(["-symMaxDiff", str(self.symMaxDiff)])
        args.extend(["-n", "-1"])
        args.extend(["-deg1"])
        return lgfString, args
    
    def _writeLGF(self):
        
        # build lgf file with molecule data
        nodes = "@nodes\n" \
                "partial_charge    label    label2    atomType    coordX    coordY    coordZ    initColor\n"
        edges = "@edges\n" \
                "                    label\n"
        attributes = "@attributes\n"
        for atm in self.data.atoms.values():
            nodes += "%8.3f %8s %8s %8s %8.4f %8.4f %8.4f %8s\n" % (atm["charge"], atm["index"], atm["symbol"], atm["iacm"], atm["ocoord"][0],atm["ocoord"][1],atm["ocoord"][2],atm["index"]) 
        
        count = 0
        for bond in self.data.bonds:
            count += 1
            edges += "%s\t%s\t%s\n" % (bond['atoms'][0], bond['atoms'][1], count)
        
        lgfString = "%s\n%s\n%s" % (nodes, edges, attributes)
        return lgfString.strip()
    
            
def _run(args, stdin, errorLog=None):
    tmp = tempfile.TemporaryFile(bufsize=0)
    tmp.write(stdin)
    tmp.seek(0)

    proc = subprocess.Popen(args, stdin=tmp, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    stdout, stderr = proc.communicate()
    tmp.close()
    if errorLog is not None:
        errorLog.debug(stderr)
    return stdout

class CGPData(object):
    atoms      = {}
    bonds     = []
    symgroups = {}
    
    def __init__(self, pdbStr, mtbStr):
        self._readPDB(pdbStr)
        self._readMTB(mtbStr)
    
    def __getitem__(self, atomid):
        '''return an atom with atomid. '''
        assert type(atomid) == int, 'Atom identifiers are integers'
        try:
            return self.atoms[atomid]
        except:
            raise Exception, 'atom with id number %d not found.' % atomid
    
    def _addAtomData(self, index, symbol=None, aType=None, charge=None, coord=None, cg=None):
        
        index = int(index)
        if index not in self.atoms.keys():
            self.atoms[index] = {}
        
        self.atoms[index]["index"] = index
        if symbol is not None:
            self.atoms[index]["symbol"] = symbol
        if aType is not None:
            self.atoms[index]["iacm"] = aType
        if charge is not None:
            self.atoms[index]["charge"] = float(charge)
        if coord is not None:
            self.atoms[index]["ocoord"] = coord
        if cg is not None:
            self.atoms[index]["cg"] = cg
            
    def _addBondData(self, atm1, atm2):
        b = {"atoms":[int(atm1),int(atm2)]}
        if b not in self.bonds:
            self.bonds.append(b)
        
    
    def _readMTB(self, mtbStr):
        
        # fill data object with necessary data in mtbStr
        readRnme = False
        
        readAtoms = False;
        nAtoms = -1;
        readBonds = False;
        nBonds = -1
        cgCounter = 1
        
        for line in mtbStr.splitlines():
            if line.startswith("#"):
                continue
            elif "MTBUILDBLSOLUTE" in line:
                readRnme = True
            elif readRnme:
                readAtoms = True
                readRnme = False
            elif readAtoms and nAtoms < 0:
                nAtoms = int(line.split()[0])
            elif readAtoms and nAtoms > 0:
                cols = line.split()
                self._addAtomData(index=cols[0], symbol=cols[1], aType=cols[2], charge=cols[4], cg=cgCounter)
                if cols[5] == "1":
                    cgCounter += 1
                nAtoms -= 1
            elif not readBonds and nAtoms == 0:
                nBonds = int(line.strip())
                readBonds = True
            elif readBonds and nBonds > 0:
                cols = line.split()
                self._addBondData(cols[0], cols[1])
                nBonds -= 1
            elif readBonds and nBonds == 0:
                readBonds = False
                break
            
        if not nAtoms == nBonds == 0:
            sys.stderr.write("ERROR: There was a problem parsing mtb file!\n")
    
    def _readPDB(self, pdbStr):
        
        # fill data object with necessary data in PDBStr
        for line in pdbStr.splitlines():
            if len(line) > 6 and line[0:6] == "HETATM":
                cols = line.split()
                self._addAtomData(cols[1], coord=[float(k)/10.0 for k in cols[5:8]])
            elif len(line) > 2 and line[0:3] == "END":
                break

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
            sys.exit(1)
            
        # Try and open mtb file
        try:
            fh = open(opts.mtb, "r")
            mtbStr = fh.read()
            fh.close()
            
        except:
            print "ERROR: there was a problem with the MTB file."
            parser.print_help()
            sys.exit(1)
    
    nautyInterface = NautyInterface(pdbStr=pdbStr, mtbStr=mtbStr)
    # Run symmetrization
    print nautyInterface.calcSym()
    

if __name__=="__main__":
    parseCommandline()
