from typing import List
from numpy import array, cross, dot, arctan2, degrees, sqrt
from numpy.linalg import norm

from chemical_equivalence.helpers.types_helpers import Atom

EQUIVALENCE_CLASS_KEY = 'equivalenceGroup'

def is_sp2_C_or_N_atom(atom: Atom) -> bool:
    return is_sp2_atom(atom) and is_element(atom, ['C', 'N'])

def is_sp2_atom(atom: Atom) -> bool:
    return has_N_neighbours(atom, N=3)

def is_element(atom: Atom, elements: List[str]) -> bool:
    return atom["type"].upper() in elements

def has_N_neighbours(atom: Atom, N: int) -> bool:
    return "conn" in atom and len(atom["conn"]) == N


def atom_coord_key(atom: Atom) -> str:
    return "ocoord" if "ocoord" in atom else "coord"

def atom_distance(atom1: Atom, atom2: Atom) -> float:
    p1, p2 = atom1[atom_coord_key(atom1)], atom2[atom_coord_key(atom1)]
    return sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2)


def dihedral_angle(*atoms: List[Atom]) -> float:
    q1, q2, q3, q4 = map(lambda atom: array(atom[atom_coord_key(atom)]), atoms)

    v1 = q2 - q1
    v2 = q3 - q2
    v3 = q4 - q3

    p1 = cross(v1, v2)
    p2 = cross(v2, v3)

    return float(
        degrees(
            arctan2(
                dot(cross(p1, p2), v2 / norm(v2)),
                dot(p1, p2),
            ),
        ),
    )

