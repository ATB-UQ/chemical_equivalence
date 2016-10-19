from typing import List
from math import sqrt

from chemical_equivalence.helpers.types import Atom

def are_atoms_chemically_equivalent(atom1: Atom, atom2: Atom) -> bool:
    atom1EquivalenceGroup = atom1["equivalenceGroup"]
    atom2EquivalenceGroup = atom2["equivalenceGroup"]
    if any([eqGroup == -1 for eqGroup in [atom1EquivalenceGroup, atom2EquivalenceGroup]]):
        return False
    else:
        return atom1EquivalenceGroup == atom2EquivalenceGroup

def atom_names(atom_list: List[Atom]) -> str:
    return ', '.join([x['symbol'] for x in atom_list])

def atoms_with_indices(atoms: List[Atom], indices: List[int]) -> List[Atom]:
    return [
        atom
        for atom in atoms
        if atom['index'] in indices
    ]

def neighbouring_atoms(atom: Atom, atoms: List[Atom]) -> List[Atom]:
    return atoms_with_indices(atoms, atom['conn'])

def is_sp2_carbon_atom(atom: Atom) -> bool:
    return is_carbon(atom) and has_N_neighbours(atom, N=3)

def is_sp3_atom(atom: Atom) -> bool:
    return has_N_neighbours(atom, N=4)

def is_bonded_to_sp2_carbon(atom: Atom, atoms: List[Atom]) -> bool:
    return any(
        [
            is_sp2_carbon_atom(atoms[neighbourAtomID])
            for neighbourAtomID in atom["conn"]
        ]
    )

def is_carbon(atom: Atom) -> bool:
    return atom["type"].upper() == "C"

def has_N_neighbours(atom: Atom, N: int) -> bool:
    return len(atom["conn"]) == N

def are_neighbours(atom1: Atom, atom2: Atom) -> bool:
    return any([True for index in atom1['conn'] if atom2["index"]==index ])

def atom_distance(atom1: Atom, atom2: Atom) -> float:
    coord_key = "ocoord" if "ocoord" in atom1 else "coord"
    p1, p2 = atom1[coord_key], atom2[coord_key]
    return sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2)
