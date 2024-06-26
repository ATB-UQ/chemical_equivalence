import subprocess
import tempfile
from typing import Union, Optional, List, Dict, Tuple
from itertools import groupby
from operator import itemgetter

from chemical_equivalence.helpers.types_helpers import Logger
from chemical_equivalence.helpers.atoms import EQUIVALENCE_CLASS_KEY
from chemical_equivalence.helpers.iterables import concat
from chemical_equivalence.config import NAUTY_EXECUTABLE

from atb_outputs.helpers.types_helpers import MolData

atb_to_nauty = lambda x: (x - 1)
nauty_to_atb = lambda x: (x + 1)

LARGE_NUMBER = 1000

Partition = Dict[int, List[int]]


def generate_nauty_edges_str(mol_data: MolData) -> str:
    return ''.join(
        [
            "{0}:{1};".format(*[atb_to_nauty(index) for index in bond['atoms']])
            for bond in mol_data.bonds
        ]
    )


def nauty_graph(mol_data: MolData, nauty_node_partition: Optional[Partition] = None) -> str:
    return 'n={num_atoms} g {edges}.f=[{node_partition}]'.format(
        num_atoms=len(mol_data.atoms),
        edges=generate_nauty_edges_str(mol_data),
        node_partition=nauty_partition_str_for(
            get_partition_for(mol_data) if nauty_node_partition is None else nauty_node_partition)
    )


def atom_descriptor_key_for(atom: Dict[str, Union[str, int, float]]) -> str:
    return 'iacm' if 'iacm' in atom else 'type'


def atom_descriptor_for(atom: Dict[str, Union[str, int, float]]) -> str:
    base_atom_descriptor = str(atom[atom_descriptor_key_for(atom)])

    if "flavour" in atom:
        # Append flavour to the atom descriptor in order to distinguish between stereoheterotopic atoms.
        # LARGE_NUMBER is chosen as a large number that is much greater than the 80 existing atom types
        return ''.join([base_atom_descriptor, str(LARGE_NUMBER + atom["flavour"])])
    else:
        return base_atom_descriptor


def nauty_partition_str_for(partition: Partition) -> str:
    # Format it in dreadnaut's partition format. Ex: "1,2,3|4,5,6"
    return '|'.join(
        [
            ','.join(map(str, indices))
            for (_, indices) in sorted(partition.items())
        ]
    )


def partition_for_chemical_equivalence_dict(chemical_equivalence_dict: Dict[int, int]) -> Partition:
    return {
        key: [atom_id for (atom_id, equivalence_class_id) in group]
        for (key, group) in groupby(
            sorted(
                chemical_equivalence_dict.items(),
                key=itemgetter(1)
            ),
            key=itemgetter(1),
        )
    }


def get_partition_for(mol_data: MolData) -> Partition:
    # A parition is a dictionnary where keys are iacm or element type (ex:12 for C) and values are a list of matching atom indexes. 
    # Ex: {'12': [2, 4, 7, 10, 13, 16], '20': [1, 3, 5, 6, 8, 9, 11, 12, 14, 15, 17, 18]}
    return {
        # Shift atom indexes by one to match dreadnaut's convention (starts at 0)
        # Atoms are sorted by key=atom_descriptor_for for canonical flovouring of the nauty nodes
        group_key: [atb_to_nauty(atom['id']) for atom in group_iterator]
        for (group_key, group_iterator) in groupby(
            sorted(mol_data.atoms.values(), key=atom_descriptor_for),
            key=atom_descriptor_for,
        )
    }


def get_nauty_node_partition(mol_data: MolData) -> str:
    return nauty_partition_str_for(
        get_partition_for(mol_data),
    )


def generate_nauty_input_from_moldata(mol_data: MolData, log: Optional[Logger] = None) -> str:
    input_str = '{nauty_graph} c xo'.format(
        nauty_graph=nauty_graph(mol_data),
    )

    if log:
        log.debug('Nauty input: {0}'.format(input_str))

    return input_str


def generate_nauty_output_from_inputstr(nauty_input_str: str, log: Optional[Logger] = None) -> str:
    nauty_stdout = _run(
        [NAUTY_EXECUTABLE],
        nauty_input_str,
        log=log,
    )

    return nauty_stdout


HAS_FOUND_ISOMORPHISM_MSG = "h and h' are identical."


def get_partition_from_nauty_output(nauty_output_str: str) -> list[tuple[int, ...]]:
    assert HAS_FOUND_ISOMORPHISM_MSG in nauty_output_str, nauty_output_str

    return [
        tuple(map(int, field.split('-')))
        for field in nauty_output_str.split(HAS_FOUND_ISOMORPHISM_MSG)[1].strip().split()
    ]


def calcEquivGroups(mol_data: MolData, log: Optional[Logger] = None) -> Union[str, Dict[int, int]]:
    if log:
        log.debug("Running Nauty")

    nauty_stdout = generate_nauty_output_from_inputstr(
        generate_nauty_input_from_moldata(mol_data, log=log),
        log=log,
    )

    if len(nauty_stdout) == 0:
        if log is not None:
            log.warning("calcEquivGroups: dreadnaut produced no output")
        return ""
    else:

        equivalence_for_atom = nauty_equivalence_dict(nauty_stdout)

        for atom_id in mol_data.atoms.keys():
            mol_data.atoms[atom_id][EQUIVALENCE_CLASS_KEY] = equivalence_for_atom[atom_id]

        # check for the case the object was initiated from a molecule3D then remap the originids on the output
        if "Molecule3D_id" in list(mol_data.atoms.values())[0]:
            equivalence_for_atom = {mol_data.atoms[a]['Molecule3D_id']: symmetry_group_id
                                    for a, symmetry_group_id in equivalence_for_atom.items()}
        if log:
            log.debug("Equivalence groups:\n{0}".format(pretty_equivalence_class_dict(mol_data, equivalence_for_atom)))

        return equivalence_for_atom


def nauty_equivalence_dict(nauty_stdout: str) -> Dict[int, int]:
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

    try:
        equivalence_groups = [
            eval_group_str(group_str)
            for group_str in orbital_data.split(";")
            if group_str
        ]
    except:
        print(orbital_data)
        raise

    equivalence_for_atom = {
        nauty_to_atb(index): n
        for n, list_of_indices in enumerate(equivalence_groups)
        for index in list_of_indices
    }

    return equivalence_for_atom


def pretty_equivalence_class_dict(mol_data: MolData, equivalence_dict: Dict[int, int]) -> str:
    return '\n'.join(
        "{equivalence_class}: {atoms}".format(
            equivalence_class=equivalence_class,
            atoms=' '.join([atom['symbol'] for atom in mol_data.atoms.values() if
                            atom[EQUIVALENCE_CLASS_KEY] == equivalence_class])
        )
        for equivalence_class in sorted(set(equivalence_dict.values()))
    )


def _run(args: List[str], stdin: str, log: Optional[Logger] = None) -> str:
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
