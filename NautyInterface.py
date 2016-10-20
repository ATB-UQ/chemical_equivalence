import subprocess
import tempfile
from typing import Union, Any, Optional, List, Dict
from itertools import groupby

from chemical_equivalence.helpers.types import Logger

NAUTY_EXECUTABLE = '/usr/local/bin/dreadnaut'

atb_to_nauty = lambda x: (x - 1)
nauty_to_atb = lambda x: (x + 1)

LARGE_NUMBER = 1000

NO_EQUIVALENCE_VALUE = -1

class NautyInterface(object):
    def __init__(self, molData: Any) -> None:
        self.data = molData

    def calcEquivGroups(self, log: Optional[Any] = None) -> Union[str, None]:
        if log: log.debug("Running Nauty")

        nauty_stdout = _run(
            [NAUTY_EXECUTABLE],
            self.nauty_input(log=log),
            log=log,
        )

        if len(nauty_stdout) == 0:
            if log is not None:
                log.warning("calcEquivGroups: dreadnaut produced no output")
            return ""

        self._procNautyOutput(nauty_stdout)

        if log:
            log.debug("Equivalence groups: {0}".format(self._getLogInfo()))

        return None

    def _getLogInfo(self) -> str:
        output = ""
        for grpID, atomsIndexs in list(self.data.equivalenceGroups.items()):
            atmNames = [self.data[self.data.get_id(i)]["symbol"] for i in atomsIndexs]
            output += "\n{0}: {1}".format(str(grpID), " ".join(atmNames))
        return output

    def _procNautyOutput(self, nauty_stdout: str) -> None:
        orbitalData = nauty_stdout.split("seconds")[-1].strip()

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
                self.data.equivalenceGroups[len(self.data.equivalenceGroups)] = [nauty_to_atb(x) for x in expandedEqGroup]

        for atmID, atm in list(self.data.atoms.items()):
            found = False
            for eqGrpID, eqGrp in list(self.data.equivalenceGroups.items()):
                if atm["index"] in eqGrp:
                    self.data.atoms[atmID]["equivalenceGroup"] = int(eqGrpID)
                    found = True
                    break
            if not found:
                self.data.atoms[atmID]["equivalenceGroup"] = NO_EQUIVALENCE_VALUE
        return None

    def nauty_input(self, log: Optional[Logger] = None) -> str:
        input_str = 'n={num_atoms} g {edges}.f=[{node_partition}] xo'.format(
            num_atoms=len(self.data.atoms),
            edges=self.nauty_edges(),
            node_partition=self.nauty_node_partition(),
        )

        if log:
            log.debug('Nauty input: {0}'.format(input_str))

        return input_str

    def nauty_edges(self) -> str:
        return ''.join(
            [
                "{0}:{1};".format(*[atb_to_nauty(index) for index in bond['atoms']])
                for bond in self.data.bonds
            ]
        )

    def nauty_node_partition(self) -> str:
        # atom_types is a dictionnary where keys are iacm or element type (ex:12 for C) and values are a list of matching atom indexes. 
        # Ex: {'12': [2, 4, 7, 10, 13, 16], '20': [1, 3, 5, 6, 8, 9, 11, 12, 14, 15, 17, 18]}

        def atom_descriptor_key_for(atom: Dict[str, Any]) -> str:
            return 'iacm' if 'iacm' in atom else 'type'

        def atom_descriptor_for(atom: Dict[str, Any]) -> str:
            base_atom_descriptor = atom[atom_descriptor_key_for(atom)]

            if "flavour" in atom:
                # Append flavour to the atom descriptor in order to distinguish between stereoheterotopic atoms.
                # LARGE_NUMBER is chosen as a large number that is much greater than the 80 existing atom types
                return ''.join([base_atom_descriptor, str(LARGE_NUMBER + atom["flavour"])])
            else:
                return base_atom_descriptor

        # Accumulate atom indexes
        atom_types = dict(
            [
                # Shift atom indexes by one to match dreadnaut's convention (starts at 0)
                # Atoms are sorted by key=atom_descriptor_for for canonical flovouring of the nauty nodes
                (group_key, [atb_to_nauty(atom['index']) for atom in group_iterator])
                for (group_key, group_iterator) in groupby(
                    sorted(self.data.atoms.values(), key=atom_descriptor_for),
                    key=atom_descriptor_for,
                )
            ]
        )

        # Format it in dreadnaut's partition format. Ex: "1,2,3|4,5,6"
        return '|'.join(
            [
                ','.join(map(str, indices))
                for (_, indices) in sorted(atom_types.items())
            ]
        )

def _run(args: List[str], stdin: str, log: Optional[Any] = None) -> str:
    tmp = tempfile.TemporaryFile(buffering=0)
    tmp.write(stdin.encode())
    tmp.seek(0)

    proc = subprocess.Popen(args, stdin=tmp, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()

    tmp.close()
    if stderr and log:
        log.debug(stderr)
    return stdout.strip().decode()

if __name__ == '__main__':
    pass
