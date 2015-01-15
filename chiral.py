def getNeighboursEquivalenceGroups(atom, molData):
        # Get back the equivalence groups of the neighbours and their id's
        listIDEq = [(neighbour, molData.atoms[neighbour]["symGr"]) for neighbour in atom["conn"] if molData.atoms[neighbour]["symGr"] != -1 ]
        # return dictionary of: id's -> equivalence groups
        return dict(listIDEq) 


        
def hasStereogenicCenters(molData):
    for atom in molData.atoms.values():
        # Chiral centers should have at least three neighbours
        if len(atom['conn']) <= 2: continue
        
        # Get back the equivalence groups of the neighbours
        neighbours_equivalence_groups = getNeighboursEquivalenceGroups(atom, molData)
        # If not two atoms in the same equivalence group
        if not [x for x, y in counter(neighbours_equivalence_groups.values()).items() if y >= 2]: 
            return True
    return False


# Pseudo-code for poisoning groups

def wrongChemicalEquivalencies(molData):
    should_rerun = False
    # If there is a chiral center in a molecule
    
    if hasStereogenicCenters(molData):
        print "HAS STEREOGENIC CENTERS. LOOKING FOR DIASTEREOTOPIC ATOMS."
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
                    print "FOUND 2 DIASTEROTOPIC ATOMS: {0}".format(atomIDs)
                    for i, atomID in enumerate(atomIDs):
                        molData[atomID]["flavour"] = i
                    should_rerun = True
        return should_rerun
    else:
        print "HAS NO STEREOGENIC CENTERS."
        return should_rerun

def counter(iterable):
    retDict = {}
    for i in iterable:
        retDict.setdefault(i, 0)
        retDict[i] += 1
    return retDict


# Called like this
#   run nauty()
#   while wrongChemicalEquivalencies()
#       run nauty()
