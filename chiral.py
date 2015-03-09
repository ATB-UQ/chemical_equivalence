def getNeighboursEquivalenceGroups(atom, molData):
        # Get back the equivalence groups of the neighbours and their id's
        listIDEq = [(neighbour, molData.atoms[neighbour]["equivalenceGroup"]) for neighbour in atom["conn"] if molData.atoms[neighbour]["equivalenceGroup"] != -1 ]
        # return dictionary of: id's -> equivalence groups
        return dict(listIDEq) 
        
def hasStereogenicCenters(molData, log):
    hasCenter = False
    for atom in molData.atoms.values():
        # Chiral centers should have at least three neighbours
        if len(atom['conn']) <= 2: continue
        
        # Get back the equivalence groups of the neighbours
        neighbours_equivalence_groups = getNeighboursEquivalenceGroups(atom, molData)
        # If not two atoms in the same equivalence group
        if not [x for x, y in counter(neighbours_equivalence_groups.values()).items() if y >= 2]: 
            hasCenter = True
            if log: log.info("Stereogenic atom: {0}".format(atom["symbol"]))
    return hasCenter

def containsDiastereotopicAtoms(molData, flavourCounter, log):
    should_rerun = False
    
    if not hasStereogenicCenters(molData, log): 
        return should_rerun
    
    if log: log.info("HAS STEREOGENIC CENTERS. NOW LOOKING FOR DIASTEREOTOPIC ATOMS.")
    # For every atom 'atm'
    for atom in molData.atoms.values():
        # Get back the equivalence groups of the neighbours
        neighbours_equivalence_groups = getNeighboursEquivalenceGroups(atom, molData)
        for equivGrpID, occurence in counter(neighbours_equivalence_groups.values()).items():
            # If there are exactly two equivalent neighbours bonded to the 'atm' atom 
            if occurence == 2:
                # Then these two neighbours are actually diasterotopic and are therefore not chemically equivalent
                # Therefore, they should manually be made non equivalent
                # And the chemical equivalency algorithm should be re-run to avoid atoms further down the graph be considered equivalent
                atomIDs = [atomID for atomID, grpID in neighbours_equivalence_groups.items() if equivGrpID==grpID]
                if log: log.info("FOUND 2 DIASTEROTOPIC ATOMS: {0}".format([molData[a]["symbol"] for a in atomIDs]))
                for atomID in atomIDs:
                    molData[atomID]["flavour"] = flavourCounter.getNext()
                should_rerun = True
    return should_rerun

def counter(iterable):
    retDict = {}
    for i in iterable:
        retDict.setdefault(i, 0)
        retDict[i] += 1
    return retDict
