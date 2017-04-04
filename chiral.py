from typing import List, Dict
from logging import Logger
from itertools import groupby

from chemical_equivalence.helpers.types_helpers import Atom, MolData, FlavourCounter, Tuple
from chemical_equivalence.helpers.atoms import is_sterogenic_atom, EQUIVALENCE_CLASS_KEY, are_substituents, neighbouring_atoms, flavour_atoms, is_sp3_atom
from chemical_equivalence.helpers.iterables import concat

MINIMUM_IDENTICAL_NEIGHBOUR_COUNT_FOR_CHIRAL = 2

def pair_of_stereo_heterotopic_atoms_for(atom: Atom, atoms: List[Atom]) -> List[Atom]:
    def by_equivalence_class(atom: Atom) -> int:
        return atom[EQUIVALENCE_CLASS_KEY]

    groups_of_atoms = [
        list(group_of_atoms)
        for (equivalence_class, group_of_atoms) in groupby(
            sorted(neighbouring_atoms(atom, atoms), key=by_equivalence_class),
            key=by_equivalence_class,
        )
    ]

    if sum([1 for group_of_atoms in groups_of_atoms if len(group_of_atoms) == MINIMUM_IDENTICAL_NEIGHBOUR_COUNT_FOR_CHIRAL]) == 1:
        return [
            group_of_atoms
            for group_of_atoms in groups_of_atoms
            if len(group_of_atoms) == MINIMUM_IDENTICAL_NEIGHBOUR_COUNT_FOR_CHIRAL
        ][0]
    else:
        return []

def get_pairs_of_stereo_heterotopic_atoms(molData: MolData) -> List[Tuple[Atom, Atom]]:
    atoms = molData.atoms.values()

    return filter(
        bool,
        [
            pair_of_stereo_heterotopic_atoms_for(atom, atoms)
            for atom in atoms
            if is_sp3_atom(atom)
        ],
    )

def contains_stereo_heterotopic_atoms(molData: MolData, flavourCounter: FlavourCounter, log: Logger) -> bool:
    sterogenic_atoms = [
        atom
        for atom in molData.atoms.values()
        if is_sterogenic_atom(atom, molData)
    ]
    if len(sterogenic_atoms) == 0:
        should_rerun = False
    else:
        # Has at least one stereogeneic atom
        pairs_of_stereo_heterotopic_atoms = get_pairs_of_stereo_heterotopic_atoms(molData)

        atoms_to_flavour = sterogenic_atoms + concat(map(list, pairs_of_stereo_heterotopic_atoms))

        if len(atoms_to_flavour) == 0:
            should_rerun = False
        else:
            should_rerun = flavour_atoms(atoms_to_flavour, flavourCounter)

            if log:
                if should_rerun:
                    log.debug(
                        "Stereogenic atoms: {0}".format(
                            [atom["symbol"] for atom in sterogenic_atoms],
                        ),
                    )
                    log.debug("HAS AT LEAST ONE STEREOGENIC ATOM. NOW LOOKING FOR STEREOHETEROTOPIC ATOMS.")

                    if len(sterogenic_atoms) >= 2:
                        log.warning("KNOWN ISSUES: if there are two or more groups of stereoheterotopic atoms then chemical equivalence depends on chiral configuration (R or S)")

                    for pair_of_atoms in pairs_of_stereo_heterotopic_atoms:
                        log.debug(
                            "FOUND 2 STEREOHETEROTOPIC ATOMS: {0}".format([atom["symbol"] for atom in pair_of_atoms])
                        )
                else:
                    log.debug('No rerun necessary.')

    return should_rerun
