

def containsEquivalenceBreakingDoubleBond(molData, flavourCounter, log=None):
    connected_sp2_carbons = connectedSp2Carbons(molData.atoms, log=log)
    grouped_sp2_carbon_pairs = groupedCarbonsByPairs(connected_sp2_carbons, molData.atoms)
    print "Grouped sp2 carbon pairs: {0}".format(grouped_sp2_carbon_pairs)
    return False

def connectedSp2Carbons(atoms, log=None):
    connected_sp2_carbons = []
    for atom in atoms.values():
        if isSp2CarbonAtom(atom) and isConnectedToSp2Carbon(atom, atoms):
            connected_sp2_carbons.append(atom)
    if log:
        log.info("Found the following carbons atoms in double bonds: {0}".format(" ".join([a["symbol"] for a in connected_sp2_carbons ])))
    print " ".join([a["symbol"] for a in connected_sp2_carbons ])
    print "Connected sp2 carbons: {0}".format(map(lambda x:x['symbol'],connected_sp2_carbons))
    return connected_sp2_carbons

def groupedCarbonsByPairs(connected_sp2_carbons, atoms):
    grouped_sp2_carbon_pairs = []
    for i, atom1 in enumerate(connected_sp2_carbons):
        for j,atom2 in enumerate(connected_sp2_carbons[i+1:]) :
            if areNeighbours(atom1,atom2):
                print "{0} and {1} are neightbours".format(atom1['symbol'], atom2['symbol'])
                grouped_sp2_carbon_pairs.append([atom1['symbol'], atom2['symbol']])
                break
    return grouped_sp2_carbon_pairs

def areNeighbours(atom1, atom2):
    return any([True for index in atom1['conn'] if atom2["index"]==index ])

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
    
