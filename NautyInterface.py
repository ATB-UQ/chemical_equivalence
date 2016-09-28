import subprocess
import tempfile

class NautyInterface(object):

    def __init__(self, molData):
        self.data = molData

    def calcEquivGroups(self, log=None):
        if log: log.debug("Running Nauty")
        nautyInput = self._writeNautyInput()

        args = ["/usr/local/bin/dreadnaut"]
        nautyStdout = _run(args, nautyInput, errorLog=log)

        if len(nautyStdout) == 0:
            if log is not None:
                log.warning("calcEquivGroups: dreadnaut produced no output")
            return ""
        self._procNautyOutput(nautyStdout)

        if log: 
            log.debug("Equivalence groups: {0}".format(self._getLogInfo()))

    def _getLogInfo(self):
        output = ""
        for grpID, atomsIndexs in list(self.data.equivalenceGroups.items()):
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
                    start, stop = list(map(int,element.split(":")))
                    expandedEqGroup.extend(list(range(start, stop + 1)))
                else:
                    expandedEqGroup.append(int(element))

            # append sym group and shift indexes up by 1
            if len(expandedEqGroup) > 1:
                self.data.equivalenceGroups[len(self.data.equivalenceGroups)] = [x+1 for x in expandedEqGroup]

        for atmID, atm in list(self.data.atoms.items()):
            found = False
            for eqGrpID, eqGrp in list(self.data.equivalenceGroups.items()):
                if atm["index"] in eqGrp:
                    self.data.atoms[atmID]["equivalenceGroup"] = int(eqGrpID)
                    found = True
                    break
            if not found:
                self.data.atoms[atmID]["equivalenceGroup"] = -1
        return


    def _writeNautyInput(self):
        return  'n={numAtoms} g {nautyGraph}.'\
                'f=[{partition}] xo'.format(**{"numAtoms":len(self.data.atoms),
                                             "nautyGraph":self.genNautyGraph(),
                                             "partition":self.genNautyPartition()}
                                          )

    def genNautyGraph(self):
        graphStr = ""
        for bond in self.data.bonds:
            graphStr += "{0}:{1};".format(*[x-1 for x in bond['atoms']])
        return graphStr

    def genNautyPartition(self):
        # atomTypes is a dictionnary where keys are iacm or element type (ex:12 for C) and values are a list of matching atom indexes. 
        # Ex: {'12': [2, 4, 7, 10, 13, 16], '20': [1, 3, 5, 6, 8, 9, 11, 12, 14, 15, 17, 18]}
        atomTypes = {}

        # Accumulate atom indexes
        for atm in list(self.data.atoms.values()):
            atom_class_identifier = 'iacm' if 'iacm' in atm else 'type'
            if "flavour" in atm:
                # append flavour to atom_class_identifier in order to distinguish between stereoheterotopic atoms,
                # the 1000 is chosen as a large number that is much greater than the 80 existing atom types
                distinguishingID = 1000+atm["flavour"]
                atomTypes.setdefault("{0}{1}".format(atm[atom_class_identifier], distinguishingID), []).append(atm['index'])
            else:
                atomTypes.setdefault(atm[atom_class_identifier],[]).append(atm['index'])

        # Shift atom indexes by one to match dreadnaut's convention (starts at 0)
        for atomType in list(atomTypes.keys()):
            atomTypes[atomType] = [x-1 for x in atomTypes[atomType]]

        # Format it in dreadnaut's partition format. Ex: 1,2,3|4,5,6
        return '|'.join( [ ','.join( map(str,v) ) for v in list(atomTypes.values()) ] )

def _run(args, stdin, errorLog=None):
    tmp = tempfile.TemporaryFile(buffering=0)
    tmp.write(stdin.encode())
    tmp.seek(0)

    proc = subprocess.Popen(args, stdin=tmp, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()

    tmp.close()
    if stderr and errorLog:
        errorLog.debug(stderr)
    return stdout.strip().decode()
