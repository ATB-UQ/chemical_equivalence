from typing import List, Any, Dict
from logging import Logger

from chemical_equivalence.helpers.types import Atom
from chemical_equivalence.helpers.atoms import is_sp3_atom
from chemical_equivalence.NautyInterface import NO_EQUIVALENCE_VALUE

MINIMUM_NEIGHBOUR_COUNT_FOR_CHIRAL = 4
MINIMUM_IDENTICAL_NEIGHBOUR_COUNT_FOR_CHIRAL = 2

def getNeighboursEquivalenceGroups(atom: Atom, molData: Any) -> Dict[int, str]:
    '''Return dictionary of: id's -> equivalence groups'''
    # Get back the equivalence groups of the neighbours and their id's
    return dict(
        [
            (neighbour, molData.atoms[neighbour]["equivalenceGroup"])
            for neighbour in atom["conn"]
            if molData.atoms[neighbour]["equivalenceGroup"] != NO_EQUIVALENCE_VALUE
        ]
    )

def hasStereogenicAtom(molData: Any, log: Logger) -> bool:
    hasCenter = False
    for atom in list(molData.atoms.values()):
        if isSterogenicAtom(atom, molData):
            hasCenter = True
            if log: log.debug("Stereogenic atom: {0}".format(atom["symbol"]))
    return hasCenter

def isSterogenicAtom(atom: Atom, molData: Any) -> bool:
    return is_sp3_atom(atom) and hasAllDifferentNeighbours(atom, molData)

def hasAllDifferentNeighbours(atom: Atom, molData: Any) -> bool:
    # Get the equivalence groups of the neighbours
    neighbours_equivalence_groups = getNeighboursEquivalenceGroups(atom, molData)

    # remove case where atoms are not in any equivalence group
    neighbours_equivalence_groups_filtered = dict(
        [
            x
            for x in neighbours_equivalence_groups.items()
            if x[1] != NO_EQUIVALENCE_VALUE
        ]
    )

    # check if all atoms are in different equivalence groups
    return all(
        [
            count == 1
            for (_, count) in
            countValueGroups(neighbours_equivalence_groups_filtered).items()
        ]
    )

def hasStereoheterotopicNeighbours(neighbours_equivalence_groups: Any, log: Logger) -> bool:
    '''For molecules that have at least one stereogenic atom, then any neighbour group that has
     two, and only two, chemically equivalent atoms is distereotopic. In the case where two pairs 
    of neighbours are chemically equivalent, then they are only stereoheterotopic if the chiral configurations
    of the two atoms are different (this is an obscure and unlikely case but could occur). 
    '''
    countStereoheterotopicGroups = 0
    for _, occurence in list(countValueGroups(neighbours_equivalence_groups).items()):
        # If there are exactly two equivalent neighbours bonded to the 'atm' atom 
        if occurence == MINIMUM_IDENTICAL_NEIGHBOUR_COUNT_FOR_CHIRAL:
            countStereoheterotopicGroups += 1

    # This is the common case of stereoheterotopic atoms
    if countStereoheterotopicGroups == 1:
        return True

    # This case is incorrectly handled as we don't know whether the chiral atoms are
    # in S or R configuration. To be more accurate we would return true only if the chrality is different
    elif countStereoheterotopicGroups == 2:
        if log: log.warning("KNOWN ISSUES: if there are two groups of stereoheterotopic atoms then chemical equivalence depends on chiral configuration (R or S)")
        return True
    else:
        return False

def getStereoheterotopicAtomGroups(neighbours_equivalence_groups: Any) -> List[Any]:
    atomGroups = []
    for equivGrpID, occurence in list(countValueGroups(neighbours_equivalence_groups).items()):
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

def containsStereoheterotopicAtoms(molData: Any, flavourCounter: Any, log: Logger) -> bool:
    should_rerun = False

    if not hasStereogenicAtom(molData, log):
        return should_rerun

    if log: log.debug("HAS AT LEAST ONE STEREOGENIC ATOM. NOW LOOKING FOR STEREOHETEROTOPIC ATOMS.")
    # For every atom 'atm'
    for atom in list(molData.atoms.values()):
        # Get back the equivalence groups of the neighbours
        neighbours_equivalence_groups = getNeighboursEquivalenceGroups(atom, molData)
        if hasStereoheterotopicNeighbours(neighbours_equivalence_groups, log):
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
def countValueGroups(dictionary: Dict[Any, Any]) -> Dict[Any, int]:
    retDict = {}
    for v in list(dictionary.values()):
        retDict.setdefault(v, 0)
        retDict[v] += 1
    return retDict
