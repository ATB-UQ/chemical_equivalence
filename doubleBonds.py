from typing import Any, Optional, List, Union
from logging import Logger

from chemical_equivalence.config import DOUBLE_BOND_LENGTH_CUTOFF
from chemical_equivalence.helpers.atoms import are_atoms_chemically_equivalent, atom_names, atoms_with_indices, neighbouring_atoms, is_sp2_carbon_atom, is_carbon, is_bonded_to_sp2_carbon, atom_distance
from chemical_equivalence.helpers.types import Atom

def contains_equivalence_breaking_double_bond(molData: Any, flavourCounter: Any, log: Optional[Logger] = None) -> bool:
    should_rerun = False

    connected_sp2_carbons = connectedSp2Carbons(molData.atoms, log)
    for (atom1, atom2) in connected_sp2_carbons:
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

        requires_rerun = correctSymmetry(neighbours, flavourCounter, log)
        if requires_rerun: should_rerun = True

    return should_rerun

def correctSymmetry(neighbours: Any, flavourCounter: Any, log: Logger) -> bool:
    should_rerun = False
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
        neighbourListLeft[0]["flavour"] = flavourCounter.getNext()
        neighbourListLeft[1]["flavour"] = flavourCounter.getNext()
        should_rerun = True
    elif are_atoms_chemically_equivalent(*neighbourListRight) and not are_atoms_chemically_equivalent(*neighbourListLeft):
        if log:
            log.debug(
                '    Found asymmetric substituents on one side of the double bond that will break the symmetry of the other side: {0}'.format(atom_names(neighbourListLeft)),
            )
            log.debug(
                "    Removed chemical equivalence between {0} and {1} (other side of the double bond)".format(*[x["symbol"] for x in neighbourListRight]),
            )
        neighbourListRight[0]["flavour"] = flavourCounter.getNext()
        neighbourListRight[1]["flavour"] = flavourCounter.getNext()
        should_rerun = True
    # If they belong to the same groups, then they need to be colored
    # so that no face in more symetric than the other

    if log and not should_rerun:
        log.debug("    Double bond does NOT break chemical equivalence due to symmetry about double bond axis")

    return should_rerun

def getNeighboursExcludingOne(atom: Atom, excludedAtom: Atom, molData: Any) -> List[Atom]:
    return [
        molData.atoms[neighbourID]
        for neighbourID in atom["conn"]
        if neighbourID != excludedAtom["index"]
    ]

def connectedSp2Carbons(atoms: List[Atom], log: Logger) -> Any:
    connected_sp2_carbons = []
    for atom in list(atoms.values()):
        if is_sp2_carbon_atom(atom) and is_bonded_to_sp2_carbon(atom, atoms) and has_suitable_bond_length(atom, atoms):
            atom2 = getConnectedSp2Carbon(atom, atoms)
            if not atom2:
                continue
            doubleBondPair = (atom, atom2)
            if not alreadyAdded(doubleBondPair, connected_sp2_carbons):
                connected_sp2_carbons.append( doubleBondPair )
    if log:
        log.debug("Found the following sp2 carbon atoms in a double bond: {0}".format(" ".join(["{0}=={1}".format(a1["symbol"], a2["symbol"]) for a1, a2 in connected_sp2_carbons ])))

    return connected_sp2_carbons

def alreadyAdded(doubleBondPair: Any, connected_sp2_carbons: Any) -> bool:
    return (doubleBondPair in connected_sp2_carbons) or (doubleBondPair[::-1] in connected_sp2_carbons)

def getConnectedSp2Carbon(atom: Atom, atoms: List[Atom]) -> Union[None, Atom]:
    for neighbourAtomID in atom["conn"]:
        if is_sp2_carbon_atom(atoms[neighbourAtomID]) and suitable_bond_length(atom, atoms[neighbourAtomID]):
            return atoms[neighbourAtomID]
    return None

def has_suitable_bond_length(atom1: Atom, atoms: List[Atom]) -> bool:
    return any(
        [
            suitable_bond_length(atom1, atom2)
            for atom2 in atoms.values()
        ]
    )

def suitable_bond_length(atom1: Atom, atom2: Atom) -> bool:
    return (atom_distance(atom1, atom2) < DOUBLE_BOND_LENGTH_CUTOFF)
