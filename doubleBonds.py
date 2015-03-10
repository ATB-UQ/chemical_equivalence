

def containsEquivalenceBreakingDoubleBond(molData, flavourCounter, log):
    connectedSp2Carbons(molData.atoms)
    return False

def connectedSp2Carbons(atoms):
    doubleBondCarbons = []
    for atom in atoms.values():
        if isSp2CarbonAtom(atom) and isConnectedToSp2Carbon(atom, atoms):
            doubleBondCarbons.append(atom)
    print " ".join([a["symbol"] for a in doubleBondCarbons ])
            
def isConnectedToSp2Carbon(atom, atoms):
    for neighbourAtomID in atom["conn"]:
        if isSp2CarbonAtom(atoms[neighbourAtomID]):
            return True
    
    
def isSp2CarbonAtom(atom):
    return isCarbon(atom) and has3Neighbours(atom)

def isCarbon(atom):
    return atom["type"].upper() == "C"

def has3Neighbours(atom):
    return len(atom["conn"]) == 3
    