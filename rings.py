from doubleBonds import areAtomsChemicallyEquivalent, atomsWithIndexes, neighbouringAtoms
from build_rings import build_rings

def atomNames(atom_list):
    return ', '.join(map(lambda x: x['symbol'], atom_list))

def containsInversableRings(molData, flavourCounter, log=None):
    should_rerun = False
    all_atoms = molData.atoms.values()
    rings = molData.rings.values()
    for ring in rings:
        ring_atoms = atomsWithIndexes(all_atoms, ring['atoms'])
        if not len(ring_atoms) in [5,6] :
            continue
        if log: log.info('Found inversable ring that could disturb chemical equivalency: {0}'.format(atomNames(ring_atoms)))
        if any([ hasDifferentSubstituent(node, all_atoms, ring_atoms) for node in ring_atoms]):
            if log: log.info('Found asymetricly-substituted atom in the nodes of the ring that will break the axial/equatorial symmetry: {0}'.format( atomNames([ node for node in ring_atoms if hasDifferentSubstituent(node, all_atoms, ring_atoms)] ) ) )
            # Then break the axial/equatorial symmetry for everyone which has some
            symmetric_nodes = filter(lambda x: not hasDifferentSubstituent(x, all_atoms, ring_atoms) , ring_atoms)
            #print "Symmetric nodes: {0}".format(symmetric_nodes)
            for symmetric_node in symmetric_nodes:
                substituents = getCarbonSubstituents(symmetric_node, all_atoms, ring_atoms)
                substituents[0]['flavour'] = flavourCounter.getNext()
                substituents[1]['flavour'] = flavourCounter.getNext()
                if log: log.info("Removed chemical equivalency between {0} and {1} (axial and equatorial substituent on inversable ring)".format(*map(lambda x:x['symbol'], substituents)))
                should_rerun = True
    return should_rerun

def inversable_rings(atoms, log):
    return [ [atom for atom in atoms if atom['type'].upper() == 'C' ] ]

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
