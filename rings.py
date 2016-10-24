from typing import Dict, Any, Optional, List

from chemical_equivalence.helpers.types import MolData, FlavourCounter, Logger, Atom
from chemical_equivalence.helpers.atoms import are_atoms_chemically_equivalent, atom_names, atoms_with_indices, neighbouring_atoms, flavour_atoms

INVERSABLE_RING_SIZES = [5, 6, 7, 8]

def is_inversable_ring(ring_atoms: List[Atom], log: Logger) -> bool:
    return len(ring_atoms) in INVERSABLE_RING_SIZES

def contains_inversable_rings(molData: MolData, flavourCounter: FlavourCounter, log: Optional[Logger] = None) -> bool:
    all_atoms = list(molData.atoms.values())
    rings = list(molData.rings.values())

    def should_rerun_for_ring(ring: List[Atom]) -> bool:
        ring_atoms = atoms_with_indices(all_atoms, ring['atoms'])

        if not is_inversable_ring(ring_atoms, log):
            should_rerun = False
        else:
            should_rerun = False
            if log: log.debug('Found inversable ring that could disturb chemical equivalency: {0}'.format(atom_names(ring_atoms)))

            if any([has_different_substituent(node, all_atoms, ring_atoms) for node in ring_atoms]):
                if log:
                    log.debug('    Found asymmetrically-substituted atom(s) in the nodes of the ring that will break the axial/equatorial symmetry: {0}'.format(atom_names([ node for node in ring_atoms if has_different_substituent(node, all_atoms, ring_atoms)])))

                # Then break the axial/equatorial symmetry for everyone
                for node in ring_atoms:
                    substituents = get_ring_substituents(node, all_atoms, ring_atoms)
                    if len(substituents) == 2 :
                        should_rerun = should_rerun or flavour_atoms(substituents, flavourCounter)
                        if log:
                            log.debug("    Make heterogeneous {0} and {1} (axial and equatorial substituents on inversable ring)".format(*[x['symbol'] for x in substituents]))
                    else :
                        if log: log.debug("    There were no axial and equatorial substituent, as there were only {0} of them.".format(len(substituents)))
            else:
                if log: log.debug('    Ring did not contain any asymmetrically-substituted nodes that will break the axial-equatorial symmetry')
            return should_rerun

    if len(rings) == 0:
        should_rerun = False
    else:
        should_rerun = any([should_rerun_for_ring(ring) for ring in rings])
    return should_rerun

def get_ring_substituents(node: Atom, atoms: List[Atom], ring_atoms: List[Atom]) -> List[Atom]:
    neighbours = neighbouring_atoms(node, atoms)
    return [neighbour for neighbour in neighbours if neighbour['index'] not in [x['index'] for x in ring_atoms] ]

def has_different_substituent(node: Atom, atoms: List[Atom], ring_atoms: List[Atom]) -> bool:
    # First, get out the two ring members from the node neighbours
    substituents = get_ring_substituents(node, atoms, ring_atoms)
    # Then, return whether or not the are chemically equivalent
    if len(substituents) == 2 :
        # For sp3 atoms
        return (not are_atoms_chemically_equivalent(*substituents))
    else:
        # Atoms which are not sp3 do not break symmetry
        return False
