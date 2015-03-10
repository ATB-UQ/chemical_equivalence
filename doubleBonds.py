

def containsEquivalenceBreakingDoubleBond(molData, flavourCounter, log=None):
    should_rerun = False
    
    connected_sp2_carbons = connectedSp2Carbons(molData.atoms, log=log)
    for atom1, atom2 in connected_sp2_carbons:
        # get the neighbouring atoms to atom1 and atom2, excluding eachother
        neighbours = {atom1["index"]: getNeighboursExcludingOne(atom1, atom2, molData),
                      atom2["index"]: getNeighboursExcludingOne(atom2, atom1, molData)}

        if log: log.info("Double bond neighbourhood: ({atom1Neighbours})--{atom1}={atom2}--({atom2Neighbours})"\
                         .format(atom1=atom1["symbol"], 
                                 atom2=atom2["symbol"],
                                 atom1Neighbours=",".join(map(lambda n:n["symbol"], neighbours[atom1["index"]])),
                                 atom2Neighbours=",".join(map(lambda n:n["symbol"], neighbours[atom2["index"]])),
                                )
                         )
        if atom1["equivalenceGroup"] == atom2["equivalenceGroup"] and atom1["equivalenceGroup"] != -1:
            if log: log.info("Double bond does NOT break chemical equivalence due to symmetry about double bond axis")
            continue
        else:
            if log: log.info("Double bond breaks chemical equivalence")
            should_rerun = True
            for neighbour in [n for neighbourList in neighbours.values() for n in neighbourList]:
                neighbour["flavour"] = flavourCounter.getNext()
                
    return should_rerun

def getNeighboursExcludingOne(atom, excludedAtom, molData):
    return [ molData.atoms[neighbourID] for neighbourID in atom["conn"] if neighbourID != excludedAtom["index"] ]

def connectedSp2Carbons(atoms, log=None):
    connected_sp2_carbons = []
    for atom in atoms.values():
        if isSp2CarbonAtom(atom) and isConnectedToSp2Carbon(atom, atoms):
            doubleBondPair = (atom, getConnectedSp2Carbon(atom, atoms))
            if not alreadyAdded(doubleBondPair, connected_sp2_carbons):
                connected_sp2_carbons.append( doubleBondPair )
    if log:
        log.info("Found the following sp2 carbon atoms in a double bond: {0}".format(" ".join(["{0}=={1}".format(a1["symbol"], a2["symbol"]) for a1, a2 in connected_sp2_carbons ])))
    
    return connected_sp2_carbons

def alreadyAdded(doubleBondPair, connected_sp2_carbons):
    flattered_connected_carbons = [a for ccPair in connected_sp2_carbons for a in ccPair]
    return any([atom in flattered_connected_carbons for atom in doubleBondPair])

def areNeighbours(atom1, atom2):
    return any([True for index in atom1['conn'] if atom2["index"]==index ])

def isConnectedToSp2Carbon(atom, atoms):
    for neighbourAtomID in atom["conn"]:
        if isSp2CarbonAtom(atoms[neighbourAtomID]):
            return True
        
def getConnectedSp2Carbon(atom, atoms):
    for neighbourAtomID in atom["conn"]:
        if isSp2CarbonAtom(atoms[neighbourAtomID]):
            return atoms[neighbourAtomID]

def isSp2CarbonAtom(atom):
    return isCarbon(atom) and has3Neighbours(atom)

def isCarbon(atom):
    return atom["type"].upper() == "C"

def has3Neighbours(atom):
    return len(atom["conn"]) == 3
    
