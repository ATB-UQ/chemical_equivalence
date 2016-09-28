import math
from .config import DOUBLE_BOND_LENGTH_CUTOFF

def containsEquivalenceBreakingDoubleBond(molData, flavourCounter, log=None):
    should_rerun = False

    connected_sp2_carbons = connectedSp2Carbons(molData.atoms, log)
    for atom1, atom2 in connected_sp2_carbons:
        if log: log.debug('Found double bond that could disturb chemical equivalency: {0}'.format(atomNames([atom1, atom2])))
        # get the neighbouring atoms to atom1 and atom2, excluding eachother
        neighbours = {atom1["index"]: getNeighboursExcludingOne(atom1, atom2, molData),
                      atom2["index"]: getNeighboursExcludingOne(atom2, atom1, molData)}

        if log: log.debug("    Double bond neighbourhood: ({atom1Neighbours})--{atom1}={atom2}--({atom2Neighbours})"\
                         .format(atom1=atom1["symbol"], 
                                 atom2=atom2["symbol"],
                                 atom1Neighbours=",".join([n["symbol"] for n in neighbours[atom1["index"]]]),
                                 atom2Neighbours=",".join([n["symbol"] for n in neighbours[atom2["index"]]]),
                                )
                         )

        requiers_rerun = correctSymmetry(neighbours, flavourCounter, log)
        if requiers_rerun:  should_rerun = True

    return should_rerun

def correctSymmetry(neighbours, flavourCounter, log):
    should_rerun = False
    neighbourListLeft, neighbourListRight = list(neighbours.values())

    # Try matching them two by two
    if areAtomsChemicallyEquivalent(*neighbourListLeft) and not areAtomsChemicallyEquivalent(*neighbourListRight):
        if log:
            log.debug('    Found asymmetric substituents on one side of the double bond that will break the symmetry of the other side: {0}'.format(atomNames(neighbourListRight)
                                                                                                                                                    )
                      )
            log.debug("    Removed chemical equivalence between {0} and {1} (other side of the double bond)".format(*[x["symbol"] for x in neighbourListLeft] 
                                                                                                                    ) 
                      )
        neighbourListLeft[0]["flavour"] = flavourCounter.getNext()
        neighbourListLeft[1]["flavour"] = flavourCounter.getNext()
        should_rerun = True
    elif areAtomsChemicallyEquivalent(*neighbourListRight) and not areAtomsChemicallyEquivalent(*neighbourListLeft):
        if log:
            log.debug('    Found asymmetric substituents on one side of the double bond that will break the symmetry of the other side: {0}'.format( atomNames(neighbourListLeft)
                                                                                                                                                     )
                      )
            log.debug("    Removed chemical equivalence between {0} and {1} (other side of the double bond)".format( *[x["symbol"] for x in neighbourListRight]
                                                                                                                     ) 
                      )
        neighbourListRight[0]["flavour"] = flavourCounter.getNext()
        neighbourListRight[1]["flavour"] = flavourCounter.getNext()
        should_rerun = True
    # If they belong to the same groups, then they need to be colored
    # so that no face in more symetric than the other

    if log and not should_rerun: log.debug("    Double bond does NOT break chemical equivalence due to symmetry about double bond axis")
    return should_rerun

def areAtomsChemicallyEquivalent(atom1, atom2):
    atom1EquivalenceGroup = atom1["equivalenceGroup"]
    atom2EquivalenceGroup = atom2["equivalenceGroup"]
    if any([eqGroup == -1 for eqGroup in [atom1EquivalenceGroup, atom2EquivalenceGroup]]):
        return False
    else:
        return atom1EquivalenceGroup == atom2EquivalenceGroup

def atomNames(atom_list):
    return ', '.join([x['symbol'] for x in atom_list])

def atomsWithIndexes(atoms, indexes):
    return [ atom for atom in atoms if atom['index'] in indexes ]

def neighbouringAtoms(atom, atoms):
    return atomsWithIndexes(atoms, atom['conn'])

def getNeighboursExcludingOne(atom, excludedAtom, molData):
    return [ molData.atoms[neighbourID] for neighbourID in atom["conn"] if neighbourID != excludedAtom["index"] ]

def connectedSp2Carbons(atoms, log):
    connected_sp2_carbons = []
    for atom in list(atoms.values()):
        if isSp2CarbonAtom(atom) and isConnectedToSp2Carbon(atom, atoms) and hasSuitableBondLength(atom, atoms):
            atom2 = getConnectedSp2Carbon(atom, atoms)
            if not atom2:
                continue
            doubleBondPair = (atom, atom2)
            if not alreadyAdded(doubleBondPair, connected_sp2_carbons):
                connected_sp2_carbons.append( doubleBondPair )
    if log:
        log.debug("Found the following sp2 carbon atoms in a double bond: {0}".format(" ".join(["{0}=={1}".format(a1["symbol"], a2["symbol"]) for a1, a2 in connected_sp2_carbons ])))

    return connected_sp2_carbons

def alreadyAdded(doubleBondPair, connected_sp2_carbons):
    return doubleBondPair in connected_sp2_carbons or doubleBondPair[::-1] in connected_sp2_carbons  

def areNeighbours(atom1, atom2):
    return any([True for index in atom1['conn'] if atom2["index"]==index ])

def isConnectedToSp2Carbon(atom, atoms):
    for neighbourAtomID in atom["conn"]:
        if isSp2CarbonAtom(atoms[neighbourAtomID]):
            return True

def getConnectedSp2Carbon(atom, atoms):
    for neighbourAtomID in atom["conn"]:
        if isSp2CarbonAtom(atoms[neighbourAtomID]) and suitableBondLength(atom, atoms[neighbourAtomID]):
            return atoms[neighbourAtomID]

def isSp2CarbonAtom(atom):
    return isCarbon(atom) and has3Neighbours(atom)

def isCarbon(atom):
    return atom["type"].upper() == "C"

def has3Neighbours(atom):
    return len(atom["conn"]) == 3

def hasSuitableBondLength(atom1, atoms):
    for atom2 in list(atoms.values()):
        if suitableBondLength(atom1, atom2):
            return True

def suitableBondLength(atom1, atom2):
    coord_key = "ocoord" if "ocoord" in atom1 else "coord"
    if _dist(atom1[coord_key], atom2[coord_key]) < DOUBLE_BOND_LENGTH_CUTOFF:
        return True

def _dist(p1, p2):
    return math.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2]-p2[2])**2 )