from typing import Any, Optional, List, Union, Tuple, Dict
from logging import Logger
from functools import reduce

from chemical_equivalence.config import DOUBLE_BOND_LENGTH_CUTOFF
from chemical_equivalence.helpers.atoms import are_atoms_chemically_equivalent, atom_names, atoms_with_indices, neighbouring_atoms, is_sp2_carbon_atom, is_carbon, is_bonded_to_sp2_carbon, atom_distance, flavour_atoms
from chemical_equivalence.helpers.types import Atom, FlavourCounter, MolData

def contains_equivalence_breaking_double_bond(molData: MolData, flavourCounter: FlavourCounter, log: Optional[Logger] = None) -> bool:
    connected_sp2_carbons = pairs_of_bonded_sp2_carbon_atoms(molData.atoms, log)

    def should_rerun_for_carbon_pair(atom_pair: Tuple[Atom, Atom]):
        (atom1, atom2) = atom_pair
        if log:
            log.debug('Found double bond that could disturb chemical equivalency: {0}'.format(atom_names([atom1, atom2])))

        # Get the neighbouring atoms to atom1 and atom2, excluding eachother
        neighbours = {
            atom1["index"]: getNeighboursExcludingOne(atom1, atom2, molData),
            atom2["index"]: getNeighboursExcludingOne(atom2, atom1, molData),
        }

        if log:
            log.debug("    Double bond neighbourhood: ({atom1Neighbours})--{atom1}={atom2}--({atom2Neighbours})".format(
                atom1=atom1["symbol"],
                atom2=atom2["symbol"],
                atom1Neighbours=",".join([n["symbol"] for n in neighbours[atom1["index"]]]),
                atom2Neighbours=",".join([n["symbol"] for n in neighbours[atom2["index"]]]),
            ))

        requires_rerun = correct_symmetry(neighbours, flavourCounter, log)
        if requires_rerun:
            should_rerun = True
        else:
            should_rerun = False

    if len(connected_sp2_carbons) == 0:
        should_rerun = False
    else:
        should_rerun = any([should_rerun_for_carbon_pair(atom_pair) for atom_pair in connected_sp2_carbons])

    return should_rerun

def correct_symmetry(neighbours: Any, flavourCounter: FlavourCounter, log: Logger) -> bool:
    neighbourListLeft, neighbourListRight = list(neighbours.values())

    # Try matching them two by two
    if are_atoms_chemically_equivalent(*neighbourListLeft) and not are_atoms_chemically_equivalent(*neighbourListRight):
        if log:
            log.debug(
                '    Found asymmetric substituents on one side of the double bond that will break the symmetry of the other side: {0}'.format(
                    atom_names(neighbourListRight),
                ),
            )
            log.debug(
                "    Removed chemical equivalence between {0} and {1} (other side of the double bond)".format(*[x["symbol"] for x in neighbourListLeft]),
            )
        should_rerun = flavour_atoms(neighbourListLeft, flavourCounter)
    elif are_atoms_chemically_equivalent(*neighbourListRight) and not are_atoms_chemically_equivalent(*neighbourListLeft):
        if log:
            log.debug(
                '    Found asymmetric substituents on one side of the double bond that will break the symmetry of the other side: {0}'.format(atom_names(neighbourListLeft)),
            )
            log.debug(
                "    Removed chemical equivalence between {0} and {1} (other side of the double bond)".format(*[x["symbol"] for x in neighbourListRight]),
            )
        should_rerun = flavour_atoms(neighbourListRight, flavourCounter)
    else:
        should_rerun = False
    # If they belong to the same groups, then they need to be colored
    # so that no face in more symetric than the other

    if log and not should_rerun:
        log.debug("    Double bond does NOT break chemical equivalence due to symmetry about double bond axis")

    return should_rerun

def getNeighboursExcludingOne(atom: Atom, excludedAtom: Atom, molData: MolData) -> List[Atom]:
    return [
        molData.atoms[neighbourID]
        for neighbourID in atom["conn"]
        if neighbourID != excludedAtom["index"]
    ]

def pairs_of_bonded_sp2_carbon_atoms(atoms: Dict[int, Atom], log: Logger) -> List[Tuple[Atom, Atom]]:
    def canonise_pair(atom_1: Atom, atom_2: Atom):
        return tuple(
            sorted(
                (atom_1, atom_2),
                key=lambda atom: atom['id'],
            ),
        )

    pairs = list(
        reduce(
            lambda acc, e: acc + e,
            [
                [canonise_pair(sp2_atom_1, sp2_atom_2) for sp2_atom_2 in get_connected_sp2_carbon_atoms_of_greater_id(sp2_atom_1, atoms)]
                for sp2_atom_1 in atoms.values()
                if is_sp2_carbon_atom(sp2_atom_1)
            ],
            [],
        )
    )
    if log:
        log.debug(
            "Found the following sp2 carbon atoms in a double bond: {0}".format(
                " ".join(
                    [
                        "{0}=={1}".format(a1["symbol"], a2["symbol"])
                        for (a1, a2) in pairs
                    ],
                ),
            ),
        )

    return pairs

def get_connected_sp2_carbon_atoms_of_greater_id(atom: Atom, atoms: Dict[int, Atom]) -> List[Atom]:
    return list(
        filter(
            lambda bonded_atom: is_sp2_carbon_atom(bonded_atom) and has_suitable_double_bond_length(atom, bonded_atom) and atom['id'] < bonded_atom['id'],
            [
                atoms[bonded_atom_id]
                for bonded_atom_id in atom["conn"]
            ],
        ),
    )

def has_suitable_double_bond_length(atom1: Atom, atom2: Atom) -> bool:
    return (atom_distance(atom1, atom2) < DOUBLE_BOND_LENGTH_CUTOFF)
