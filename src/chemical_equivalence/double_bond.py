"""Double-bond detection for stereo perception. Bond orders are not available from a PDB,
so a double bond is two bonded sp2 C/N atoms (three neighbours each) closer than the
element pair's single/double threshold."""
from typing import Optional, List, Dict, Tuple
from logging import Logger

from chemical_equivalence.config import DOUBLE_BOND_LENGTH_CUTOFF
from chemical_equivalence.helpers.atoms import is_sp2_C_or_N_atom, atom_distance
from chemical_equivalence.helpers.types_helpers import Atom, MolData


def double_bonds(mol_data: MolData, log: Optional[Logger] = None) -> List[Tuple[int, int]]:
    """``[(atom id, atom id)]`` for every detected double bond, lower id first."""
    return [
        (atom_1['id'], atom_2['id'])
        for (atom_1, atom_2) in pairs_of_bonded_sp2_C_or_N_atoms(mol_data.atoms, log)
    ]


def pairs_of_bonded_sp2_C_or_N_atoms(atoms: Dict[int, Atom], log: Optional[Logger] = None) -> List[Tuple[Atom, Atom]]:
    pairs = [
        (atom_1, atom_2)
        for atom_1 in sorted(atoms.values(), key=lambda atom: atom['id'])
        if is_sp2_C_or_N_atom(atom_1)
        for atom_2 in get_connected_sp2_C_or_N_atoms_of_greater_id(atom_1, atoms)
    ]
    if log:
        log.debug(
            "Found the following sp2 [C,N] atoms in a double bond: {0}".format(
                " ".join("{0}=={1}".format(a1["symbol"], a2["symbol"]) for (a1, a2) in pairs),
            ),
        )
    return pairs


def get_connected_sp2_C_or_N_atoms_of_greater_id(atom: Atom, atoms: Dict[int, Atom]) -> List[Atom]:
    return [
        atoms[bonded_atom_id]
        for bonded_atom_id in sorted(atom["conn"])
        if is_sp2_C_or_N_atom(atoms[bonded_atom_id])
        and has_suitable_double_bond_length(atom, atoms[bonded_atom_id])
        and atom['id'] < bonded_atom_id
    ]


def has_suitable_double_bond_length(atom1: Atom, atom2: Atom) -> bool:
    return (atom_distance(atom1, atom2) < DOUBLE_BOND_LENGTH_CUTOFF[frozenset([atom1['type'], atom2['type']])])
