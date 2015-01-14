import subprocess
import tempfile

class NautyInterface(object):
    
    def __init__(self, molData):
        self.data = molData
    
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
