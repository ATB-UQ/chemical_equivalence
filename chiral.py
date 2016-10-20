from typing import List, Any, Dict
from logging import Logger

from chemical_equivalence.helpers.types import Atom, MolData, FlavourCounter
from chemical_equivalence.helpers.atoms import is_sterogenic_atom, EQUIVALENCE_CLASS_KEY

MINIMUM_NEIGHBOUR_COUNT_FOR_CHIRAL = 4
MINIMUM_IDENTICAL_NEIGHBOUR_COUNT_FOR_CHIRAL = 2

def get_neighbours_equivalence_groups(atom: Atom, molData: MolData) -> Dict[int, int]:
    '''Return dictionary of: id's -> equivalence groups'''
    # Get back the equivalence groups of the neighbours and their id's
    return dict(
        [
            (neighbour_id, molData.atoms[neighbour_id][EQUIVALENCE_CLASS_KEY])
            for neighbour_id in atom["conn"]
        ]
    )

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

def has_stereo_heterotopic_neighbours(neighbours_equivalence_groups: Any, log: Logger) -> bool:
    '''For molecules that have at least one stereogenic atom, then any neighbour group that has
     two, and only two, chemically equivalent atoms is distereotopic. In the case where two pairs 
    of neighbours are chemically equivalent, then they are only stereoheterotopic if the chiral configurations
    of the two atoms are different (this is an obscure and unlikely case but could occur). 
    '''
    countStereoheterotopicGroups = 0
    for (_, occurence) in list(count_value_groups(neighbours_equivalence_groups).items()):
        # If there are exactly two equivalent neighbours bonded to the 'atm' atom 
        if occurence == MINIMUM_IDENTICAL_NEIGHBOUR_COUNT_FOR_CHIRAL:
            countStereoheterotopicGroups += 1

    # This is the common case of stereoheterotopic atoms
    if countStereoheterotopicGroups == 1:
        return True

    # This case is incorrectly handled as we don't know whether the chiral atoms are
    # in S or R configuration. To be more accurate we would return true only if the chrality is different
    elif countStereoheterotopicGroups == 2:
        if log:
            log.warning("KNOWN ISSUES: if there are two groups of stereoheterotopic atoms then chemical equivalence depends on chiral configuration (R or S)")
        return True
    else:
        return False

def getStereoheterotopicAtomGroups(neighbours_equivalence_groups: Any) -> List[Any]:
    atomGroups = []
    for (equivGrpID, occurence) in list(count_value_groups(neighbours_equivalence_groups).items()):
        # If there are exactly two equivalent neighbours bonded to the 'atm' atom 
        if occurence == MINIMUM_IDENTICAL_NEIGHBOUR_COUNT_FOR_CHIRAL:
            atomGroups.append(
                [
                    atomID
                    for (atomID, grpID) in neighbours_equivalence_groups.items()
                    if equivGrpID == grpID
                ]
            )
    return atomGroups

def contains_stereo_heterotopic_atoms(molData: MolData, flavourCounter: FlavourCounter, log: Logger) -> bool:
    if not has_stereogenic_atom(molData, log):
        should_rerun = False
    else:
        should_rerun = False
        if log:
            log.debug("HAS AT LEAST ONE STEREOGENIC ATOM. NOW LOOKING FOR STEREOHETEROTOPIC ATOMS.")

        for atom in list(molData.atoms.values()):
            # Get back the equivalence groups of the neighbours
            neighbours_equivalence_groups = get_neighbours_equivalence_groups(atom, molData)
            if has_stereo_heterotopic_neighbours(neighbours_equivalence_groups, log):
                # Then these two neighbours are actually stereoheterotopic and are therefore not chemically equivalent
                # Therefore, they should manually be made non equivalent
                # And the chemical equivalency algorithm should be re-run to avoid atoms further down the graph be considered equivalent
                atomGroups = getStereoheterotopicAtomGroups(neighbours_equivalence_groups)
                for atomIDs in atomGroups:
                    if log: log.debug("FOUND 2 STEREOHETEROTOPIC ATOMS: {0}".format([molData[a]["symbol"] for a in atomIDs]))
                    for atomID in atomIDs:
                        molData[atomID]["flavour"] = flavourCounter.getNext()
                    should_rerun = True
    return should_rerun

# Counts the occurence of a dictionnary's keys
def count_value_groups(dictionary: Dict[Any, Any]) -> Dict[Any, int]:
    retDict = {}
    for v in list(dictionary.values()):
        retDict.setdefault(v, 0)
        retDict[v] += 1
    return retDict
