import subprocess
import tempfile
from typing import Union, Any, Optional, List, Dict
from itertools import groupby
from functools import reduce

from chemical_equivalence.helpers.types import Logger
from chemical_equivalence.helpers.atoms import EQUIVALENCE_CLASS_KEY

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

        equivalence_for_atom = self.nauty_equivalence(nauty_stdout)
        for atom_id in self.data.atoms.keys():
            self.data.atoms[atom_id][EQUIVALENCE_CLASS_KEY] = equivalence_for_atom[atom_id]

        if log:
            log.debug("Equivalence groups:\n{0}".format(self._getLogInfo(equivalence_for_atom)))

        return None

    def _getLogInfo(self, equivalence_dict: Dict[int, int]) -> str:
        return '\n'.join(
            "{equivalence_class}: {atoms}".format(
                equivalence_class=equivalence_class,
                atoms = ' '.join([atom['symbol'] for atom in self.data.atoms.values() if atom[EQUIVALENCE_CLASS_KEY] == equivalence_class])
            )
            for equivalence_class in sorted(set(equivalence_dict.values()))
        )

    def nauty_equivalence(self, nauty_stdout: str) -> Dict[int, int]:
        orbital_data = nauty_stdout.split("seconds")[-1].strip()

        def eval_group_field(group_field: str) -> List[int]:
            if ':' in group_field:
                # Defines a range
                start, stop = map(int, group_field.split(":"))
                return [int(x) for x in range(start, stop + 1)]
            else:
                return [int(group_field)]

        def eval_group_str(group_str: str) -> List[int]:
            if '(' in group_str:
                # Last element is group size '(N)', not needed
                group_fields = group_str.split()[:-1]
            else:
                group_fields = group_str.split()
            return concat([eval_group_field(group_field) for group_field in group_fields])

        equivalence_groups = [
            eval_group_str(group_str)
            for group_str in orbital_data.split(";")
            if group_str
        ]

        equivalence_for_atom = dict(
            concat(
                [
                    [(nauty_to_atb(index), n) for index in list_of_indices]
                    for (n, list_of_indices) in enumerate(equivalence_groups)
                ],
            ),
        )

        return equivalence_for_atom

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

def concat(list_of_lists: List[List[Any]]) -> List[Any]:
    return reduce(
        lambda acc, e: acc + e,
        list_of_lists,
        [],
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
