import collections

def is_chiral(atom):
    pass
    # Chiral centers should have at least three neighbours
    if len(atom['conn']) <= 2: return False;
    # Get back the equivalence groups of the neighbours
    neighbours_equivalence_groups = [-1,-1,1,2]
    # Filter out atoms in no eq group, e.g. -1's
    neighbours_equivalence_groups = filter(lambda x: x!=-1,neighbours_equivalence_groups)
    # If not two atoms in the same equivalence group
    return not [x for x, y in collections.Counter(neighbours_equivalence_groups).items() if y >= 2]


# Pseudo-code for poisoning groups

def wrongChemicalEquivalencies(nautyData):
    should_rerun = False
    # If there is a chiral center in a molecule
    print any( [is_chiral(atm) for atm in nautyData.atoms.values()] )
    if any( [is_chiral(atm) for atm in nautyData.atoms.values()] ) :
        pass
        # For every atom 'atm'
    #    for atm in atoms:
            # If there is exactly two equivalent neighbours bonded to the 'atm' atom
                # Then these two neighbours are actually diasterotopic and are therefore not chemically equivalent
                # Therefore, they should manually be made non equivalent
                # And the chemical equivalency algorithm should be re-run to avoid atoms further down the graph be considered equivalent
    #            should_rerun = True
    #    return should_rerun
    #else:
    #    return False




# Called like this
#   run nauty()
#   while wrongChemicalEquivalencies()
#       run nauty()
