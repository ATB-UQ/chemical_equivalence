from typing import List, Any, Dict
from logging import Logger
from itertools import groupby

from chemical_equivalence.helpers.types import Atom, MolData, FlavourCounter, Tuple
from chemical_equivalence.helpers.atoms import is_sterogenic_atom, EQUIVALENCE_CLASS_KEY, are_substituents

MINIMUM_IDENTICAL_NEIGHBOUR_COUNT_FOR_CHIRAL = 2

def has_stereogenic_atom(molData: MolData, log: Logger) -> bool:
    if log:
        log.debug(
            "Stereogenic atoms: {0}".format(
                [atom["symbol"] for atom in molData.atoms.values() if is_sterogenic_atom(atom, molData)],
            ),
        )
    return any(
        [
            is_sterogenic_atom(atom, molData)
            for atom in molData.atoms.values()
        ]
    )

def get_pairs_of_stereo_heterotopic_atoms(molData: MolData) -> List[Tuple[Atom, Atom]]:
    def by_equivalence_class(atom: Atom) -> int:
        return atom[EQUIVALENCE_CLASS_KEY]

    return list(
        filter(
        lambda group_of_atoms: len(group_of_atoms) == MINIMUM_IDENTICAL_NEIGHBOUR_COUNT_FOR_CHIRAL,# and are_substituents(*group_of_atoms),
            [
                list(group_of_atoms)
                for (equivalence_class, group_of_atoms) in groupby(
                    sorted(molData.atoms.values(), key=by_equivalence_class),
                    key=by_equivalence_class,
                )
            ],
        ),
    )

def contains_stereo_heterotopic_atoms(molData: MolData, flavourCounter: FlavourCounter, log: Logger) -> bool:
    if not has_stereogenic_atom(molData, log):
        should_rerun = False
    else:
        # Has at least one stereogeneic atom
        if log:
            log.debug("HAS AT LEAST ONE STEREOGENIC ATOM. NOW LOOKING FOR STEREOHETEROTOPIC ATOMS.")
            log.warning("KNOWN ISSUES: if there are two or more groups of stereoheterotopic atoms then chemical equivalence depends on chiral configuration (R or S)")

        pairs_of_stereo_heterotopic_atoms = get_pairs_of_stereo_heterotopic_atoms(molData)
        if len(pairs_of_stereo_heterotopic_atoms) == 0:
            should_rerun = False
        else:
            should_rerun = True
            for pair_of_atoms in pairs_of_stereo_heterotopic_atoms:
                if log:
                    log.debug("FOUND 2 STEREOHETEROTOPIC ATOMS: {0}".format([atom["symbol"] for atom in pair_of_atoms]))
                for stereo_heterotopic_atom in pair_of_atoms:
                    stereo_heterotopic_atom["flavour"] = flavourCounter.getNext()

    return should_rerun
