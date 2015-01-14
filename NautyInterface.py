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
            self.data = MolData(pdbStr, mtbStr)
        
    
    def calcEquivGroups(self, log=None):
        nautyInput = self._writeNautyInput()
        
        args = ["dreadnaut"]
        nautyStdout = _run(args, nautyInput, errorLog=log)
        
        if len(nautyStdout) == 0:
            if log is not None:
                log.warning("calcEquivGroups: dreadnaut produced no output")
            return ""
        
        self._procNautyOutput(nautyStdout)
        
        if log: log.info("Equivalence groups: {0}".format(self._getLogInfo()))
    
    def _getLogInfo(self):
        output = ""
        for grpID, atomsIndexs in self.data.symgroups.items():
            atmNames = [self.data[self.data.get_id(i)]["symbol"] for i in atomsIndexs]
            output += "\n{0}: {1}".format(str(grpID), " ".join(atmNames))
        return output
        
    
    def _procNautyOutput(self, nautyStdout):
        orbitalData = nautyStdout.split("seconds")[-1].strip()
        
        # last item in each group is the number of nodes and not needed 
        eqGroups = [grp.split()[:-1] for grp in orbitalData.split(";") if grp] 
        
        # expand ranges in each group
        for grp in eqGroups:
            expandedEqGroup = []
            for element in grp:
                if ":" in element:
                    start, stop = map(int,element.split(":"))             
                    expandedEqGroup.extend(range(start, stop + 1))
                else:
                    expandedEqGroup.append(int(element))
            
            # append sym group and shift indexes up by 1
            if len(expandedEqGroup) > 1:
                self.data.symgroups[len(self.data.symgroups)] = map(lambda x:x+1, expandedEqGroup)
        
        for atmID, atm in self.data.atoms.items():
            found = False
            for eqGrpID, eqGrp in self.data.symgroups.items():
                if atm["index"] in eqGrp:
                    self.data.atoms[atmID]["symGr"] = int(eqGrpID)
                    found = True
                    break
            if not found:
                self.data.atoms[atmID]["symGr"] = -1
                         
        

    def _writeNautyInput(self):
        return  'n={numAtoms} g {nautyGraph}.'\
                'f=[{partition}] xo'.format(**{"numAtoms":len(self.data.atoms),
                                             "nautyGraph":self.genNautyGraph(),
                                             "partition":self.genNautyPartition()}
                                          )
    
    def genNautyGraph(self):
        graphStr = ""
        for bond in self.data.bonds:
            graphStr += "{0}:{1};".format(*map(lambda x:x-1, bond['atoms']))
        return graphStr
    
    def genNautyPartition(self):
        # atomTypes is a dictionnary where keys are iacm's (ex:12 for C) and values are a list of matching atom indexes. 
        # Ex: {'12': [2, 4, 7, 10, 13, 16], '20': [1, 3, 5, 6, 8, 9, 11, 12, 14, 15, 17, 18]}
        atomTypes = {}

        # Accumulate atom indexes
        for atm in self.data.atoms.values():
            atomTypes.setdefault(atm['iacm'],[]).append(atm['index'])

        # Shift atom indexes by one to match dreadnaut's convention (starts at 0)
        for atomType in atomTypes.keys():
            atomTypes[atomType] = map(lambda x:x-1, atomTypes[atomType])

        # Format it in dreadnaut's partition format. Ex: 1,2,3|4,5,6
        return '|'.join( [ ','.join( map(str,v) ) for v in atomTypes.values() ] )
    
def _run(args, stdin, errorLog=None):
    tmp = tempfile.TemporaryFile(bufsize=0)
    tmp.write(stdin)
    tmp.seek(0)

    proc = subprocess.Popen(args, stdin=tmp, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    stdout, stderr = proc.communicate()
    
    tmp.close()
    if errorLog is not None:
        errorLog.debug(stderr)
    return stdout.strip()

class MolData(object):
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
    
    def _readPDB(self, string):
        '''Read lines of PDB files'''
        
        pdbDict = {}
        for line in string.splitlines():
            #lines with atom coordinates
            if line.startswith("ATOM") |line.startswith('HETATM'):      
                #split line for different fields
                #it = line.split()
                #in case of some fields are missing, read according to pdb standard
                #if len(it) < 11:        
                it = [line[0:6],line[6:11],line[12:16],line[17:20],
                        line[22:27], line[30:38],line[38:46],line[46:54],line[54:60],
                        line[60:66],line[76:78]]
                it = [ i.strip() for i in it] 
                #store information that we are interested in. Refer to MoleculeData
                #class for more details, also here we have a angstrom to nm conversion
                pdbDict[int(it[1])] = {\
                    'index':int(it[1]), 'symbol':it[2], 'group': it[3], 
                    'coord':[float(it[5])/10.,float(it[6])/10.,float(it[7])/10.],
                    'pdb':line.strip(),'type':it[-1].upper()}
            #connectivity records
            if line.startswith('CONECT'):
                it = line.split()
                for key, item in pdbDict.items():
                    if key == int(it[1]):
                        conn_list = []
                        for i in it[2:6]:
                            try:
                                num = int(i)
                                if num == 0:
                                    break
                                conn_list.append(num)
                            except Exception:
                                break
                        if item.has_key('conn'):
                            item['conn'].extend(conn_list)
                        else:
                            item['conn'] = conn_list
                        break
    
        # make sure connectivity information is symmetric by simply mirroring connections
        # flag any orphan connectivities for removal
        orphanAtomReference = {}
        for key in pdbDict: # for each atom in pdb
            if pdbDict[key].has_key('conn'): # if this atom lists connections to others
                for conn in pdbDict[key]['conn']: # for each connection to others
                    if not pdbDict.has_key(conn): # if it connects to a non-existant atom
                        print "connectivity made from atom %s to non-existant atom %s!" % (key, conn)
                        if not orphanAtomReference.has_key(key): orphanAtomReference[key] = [conn]
                        else: orphanAtomReference[key].append(conn)
                        continue
                    if pdbDict[conn].has_key('conn'): # if the other atom has a list of connections
                        pdbDict[conn]['conn'].append(key) # append this atom to the end of the other atom's list
                    else:
                        pdbDict[conn]['conn'] = [ key ] # create new list containing this atom
                        
        # remove any orphan connection records
        for k, connList in orphanAtomReference.items():
            for c in connList: 
                pdbDict[k]["conn"].remove(c)
                self.log.warn("connectivity made from %s to non-existant atom %s removed" % (k, c))
    
        #connectivity missing in pdb (i.e. one or more atoms not connected to molecule), terminate directly
        if sum([not i.has_key('conn') for i in pdbDict.values()]) != 0:
            missings = filter(lambda x:not pdbDict[x].has_key('conn'), pdbDict.keys())
            missings = [str(m) for m in missings]
    
        #sort and unique connectivities
        for i in pdbDict.values():
            i['conn'] = sorted(list(set(i['conn'])))
     
        
        self.atoms = pdbDict
        

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
    nautyInterface.calcSym()
    print nautyInterface.data.symgroups
    

if __name__=="__main__":
    #parseCommandline()
    nautyInterface = NautyInterface(open("testing/manyEqGrps.pdb").read(), open("testing/manyEqGrps.mtb").read())
    nautyInterface.calcEquivGroups()
    print nautyInterface.data.symgroups
    #print "\n".join([str((at["index"], at["symGr"])) for at in nautyInterface.data.atoms.values()])
     
