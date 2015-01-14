import collections

def is_chiral(atom_index)
    # Chiral centers should have at least three neighbours
    return false if less than 3 neighbours
    # Get back the equivalence groups of the neighbours
    equivalent_groups_neighbours = ...
    # Filter out atoms in no eq group, e.g. -1's
    equivalent_neighbours = filter(lambda x: x!=-1,)
    # If not two atoms in the same equivalence group
    return true if not [x for x, y in collections.Counter(equivalent_neighbours).items() if y >= 2]


# Pseudo-code for poisoning groups

def wrongChemicalEquivalency()
    should_rerun = False
    # If there is a chiral center in a molecule
    if any( [is_chiral(atm_index) for atm_index in atoms] ) :
        # For every atom 'atm'
        for atm in atoms:
            # If there is exactly two equivalent neighbours bonded to the 'atm' atom
                # Then these two neighbours are actually diasterotopic and are therefore not chemically equivalent
                # Therefore, they should manually be made non equivalent
                # And the chemical equivalency algorithm should be re-run to avoid atoms further down the graph be considered equivalent
                should_rerun = True
        return should_rerun
    else:
        return False




# Called like this
run nauty()
while wrongChemicalEquivalency()
    run nauty()
