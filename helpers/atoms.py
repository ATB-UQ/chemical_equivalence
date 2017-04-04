from typing import List
from math import sqrt

from chemical_equivalence.helpers.types_helpers import Atom, FlavourCounter, MolData

EQUIVALENCE_CLASS_KEY = 'equivalenceGroup'

def are_atoms_chemically_equivalent(atom1: Atom, atom2: Atom) -> bool:
    return atom1[EQUIVALENCE_CLASS_KEY] == atom2[EQUIVALENCE_CLASS_KEY]

def atom_names(atom_list: List[Atom]) -> str:
    return ', '.join([x['symbol'] for x in atom_list])

def atoms_with_indices(atoms: List[Atom], indices: List[int]) -> List[Atom]:
    return [
        atom
        for atom in atoms
        if atom['id'] in indices
    ]

def neighbouring_atoms(atom: Atom, atoms: List[Atom]) -> List[Atom]:
    return atoms_with_indices(atoms, atom['conn'])

def is_sp2_C_or_N_atom(atom: Atom) -> bool:
    return is_sp2_atom(atom) and is_element(atom, ['C', 'N'])

def is_sp2_atom(atom: Atom) -> bool:
    return has_N_neighbours(atom, N=3)

def is_sp3_atom(atom: Atom) -> bool:
    return has_N_neighbours(atom, N=4)

def is_bonded_to_sp2_C_or_N(atom: Atom, atoms: List[Atom]) -> bool:
    return any(
        [
            is_sp2_C_or_N_atom(atoms[neighbourAtomID])
            for neighbourAtomID in atom["conn"]
        ]
    )

def is_carbon(atom: Atom) -> bool:
    return is_element(atom, ['C'])

def is_element(atom: Atom, elements: List[str]) -> bool:
    return atom["type"].upper() in elements

def has_N_neighbours(atom: Atom, N: int) -> bool:
    return len(atom["conn"]) == N

def are_neighbours(atom1: Atom, atom2: Atom) -> bool:
    return any([True for atom_id in atom1['conn'] if atom2["id"] == atom_id])

def atom_distance(atom1: Atom, atom2: Atom) -> float:
    coord_key = "ocoord" if "ocoord" in atom1 else "coord"
    p1, p2 = atom1[coord_key], atom2[coord_key]
    return sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2)

def is_sterogenic_atom(atom: Atom, molData: MolData) -> bool:
    return is_sp3_atom(atom) and has_all_different_neighbours(atom, molData.atoms.values())

def are_substituents(atom_1: Atom, atom_2: Atom) -> bool:
    return (set(atom_1['conn']) & set(atom_2['conn'])) != set()

def neighbour_equivalence_classes(atom: Atom, atoms: List[Atom]) -> List[int]:
    return [
        neighbour_atom[EQUIVALENCE_CLASS_KEY]
        for neighbour_atom in neighbouring_atoms(atom, atoms)
    ]

def has_all_different_neighbours(atom: Atom, atoms: List[Atom]) -> bool:
    all_equivalence_classes = neighbour_equivalence_classes(atom, atoms)
    return len(set(all_equivalence_classes)) == len(all_equivalence_classes)

def flavour_atoms(atoms: List[Atom], flavour_counter: FlavourCounter) -> bool:
    '''
    Flavours atoms.
    Returns whether or not a rerun is necessary.
    '''
    should_rerun = any(['flavour' not in atom for atom in atoms])

    for atom in atoms:
        atom["flavour"] = flavour_counter.getNext()

    return should_rerun
