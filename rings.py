from doubleBonds import areAtomsChemicallyEquivalent, atomsWithIndexes, neighbouringAtoms, atomNames

def containsInversableRings(molData, flavourCounter, log=None):
    should_rerun = False
    all_atoms = molData.atoms.values()
    rings = molData.rings.values()
    for ring in rings:
        ring_atoms = atomsWithIndexes(all_atoms, ring['atoms'])
        
        if not is_inversable_ring(ring_atoms,log) :
            continue
        if log: log.debug('Found inversable ring that could disturb chemical equivalency: {0}'.format(atomNames(ring_atoms)))
        if any([ hasDifferentSubstituent(node, all_atoms, ring_atoms) for node in ring_atoms]):
            if log: log.debug('    Found asymmetrically-substituted atom(s) in the nodes of the ring that will break the axial/equatorial symmetry: {0}'.format( atomNames([ node for node in ring_atoms if hasDifferentSubstituent(node, all_atoms, ring_atoms)] ) ) )
            # Then break the axial/equatorial symmetry for everyone which has some
            symmetric_nodes = filter(lambda x: not hasDifferentSubstituent(x, all_atoms, ring_atoms) , ring_atoms)
            #print "Symmetric nodes: {0}".format(symmetric_nodes)
            for symmetric_node in symmetric_nodes:
                substituents = getCarbonSubstituents(symmetric_node, all_atoms, ring_atoms)
                if len(substituents) == 2 :
                    substituents[0]['flavour'] = flavourCounter.getNext()
                    substituents[1]['flavour'] = flavourCounter.getNext()
                    if log: log.debug("    Removed chemical equivalency between {0} and {1} (axial and equatorial substituents on inversable ring)".format(*map(lambda x:x['symbol'], substituents)))
                    should_rerun = True
                else :
                    if log: log.debug("    There were no axial and equatorial substituent, as there were only {0} of them.".format(len(substituents)))
        else:
            if log: log.debug('    Ring did not contain any asymmetrically-substituted nodes that will break the axial-equatorial symmetry')
    return should_rerun

def is_inversable_ring(ring_atoms, log):
    return len(ring_atoms) in [5,6,7,8]

def getCarbonSubstituents(node, atoms, ring_atoms):
    neighbours = neighbouringAtoms(node, atoms)
    return [neighbour for neighbour in neighbours if neighbour['index'] not in map(lambda x:x['index'], ring_atoms) ]

def hasDifferentSubstituent(node, atoms, ring_atoms):
    # First, get out the two ring members from the node neighbours
    substituents = getCarbonSubstituents(node, atoms, ring_atoms)
    # Then, return whether or not the are chemically equivalent
    if len(substituents) == 2 :
        # For sp3 atoms
        return not areAtomsChemicallyEquivalent(substituents[0], substituents[1])
    else:
        # Not sp3 atoms do not break symmetry
        return False
