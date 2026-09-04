"""Stereo perception for chemical equivalence: which graph automorphisms are real symmetries.

Nauty's orbits of the molecular graph give *constitutional* equivalence: two atoms are
equivalent when some relabelling preserves every bond. That is too coarse. The two
hydrogens of a CH2 next to a stereocentre are constitutionally equivalent but
diastereotopic: no rotation or reflection of the molecule, in any accessible conformer,
exchanges them. The rule that separates the two is:

    a graph automorphism g is a symmetry operation iff there is one global handedness
    e in {+1, -1} such that at every tetrahedral atom X the orientation of the mapped
    neighbours at g(X) equals e times the orientation at X, and every cis pair across a
    double bond maps onto a cis pair.

e = +1 automorphisms are proper (rotations, ring flips, bond rotations: anything reachable
by a continuous motion); e = -1 ones are improper (a mirror). Both are real symmetries for
scalar properties, so enantiotopic atoms (meso compounds, Cs molecules, the two ortho
carbons of chlorocyclohexane) stay equivalent, and only diastereotopic atoms split.
Tetrahedral orientation is a configurational invariant -- it cannot change under any
continuous deformation that keeps the centre tetrahedral -- so one static geometry is
enough, and conformational averaging (ring inversion, methyl rotation) falls out of the
rule without being modelled: the flip-plus-rotation of a chair is a proper motion whose
label permutation preserves every orientation.

This module computes the geometric data (orientations, cis pairs); NautyInterface encodes
them as gadgets on the graph and finds the orbits. The old approach -- stamping arbitrary
"flavour" integers on atoms near a stereocentre -- could not keep two symmetry-related
sites equivalent once both were corrected, which is why chlorocyclohexane used to come out
with 18 distinct atoms.

Not detected (same as before): axial chirality (allenes, atropisomers, spiro), pyramidal
inversion-locked centres with three neighbours (sulfoxides, phosphines).
"""
from collections import deque
from logging import Logger
from typing import Dict, Iterable, List, Optional, Set, Tuple

from numpy import array
from numpy.linalg import det, norm

from chemical_equivalence.config import (CIS_MAX_DIHEDRAL, TETRAHEDRAL_MIN_NORMALISED_VOLUME,
                                         TRANS_MIN_DIHEDRAL)
from chemical_equivalence.helpers.atoms import atom_coord_key, dihedral_angle
from chemical_equivalence.helpers.types_helpers import MolData

Orientation = int  # +1 or -1
CisPair = Tuple[int, int]


def _coordinates(mol_data: MolData, atom_id: int):
    atom = mol_data.atoms[atom_id]
    return array(atom[atom_coord_key(atom)], dtype=float)


def tetrahedral_orientations(mol_data: MolData, log: Optional[Logger] = None) -> Dict[int, Orientation]:
    """The handedness of every 4-coordinate atom: the sign of the volume spanned by its
    neighbours taken in increasing atom-id order. Atoms whose neighbours are (near)
    coplanar carry no orientation."""
    orientations = {}
    for (atom_id, atom) in sorted(mol_data.atoms.items()):
        neighbours = sorted(atom.get('conn', []))
        if len(neighbours) != 4:
            continue
        centre = _coordinates(mol_data, atom_id)
        positions = [_coordinates(mol_data, neighbour) for neighbour in neighbours]
        mean_bond_length = sum(norm(position - centre) for position in positions) / 4.0
        if mean_bond_length <= 0:
            continue
        volume = float(det(array([positions[0] - positions[3],
                                  positions[1] - positions[3],
                                  positions[2] - positions[3]])))
        if abs(volume) / mean_bond_length ** 3 < TETRAHEDRAL_MIN_NORMALISED_VOLUME:
            if log:
                log.debug('Atom {0} has four coplanar neighbours; no tetrahedral orientation'.format(atom['symbol']))
            continue
        orientations[atom_id] = 1 if volume > 0 else -1
    return orientations


def acyclic_double_bonds(mol_data: MolData, double_bonds: Iterable[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """The detected double bonds whose two atoms share no ring. A double bond inside a
    ring carries no stereo information beyond the graph -- each of its atoms has two ring
    neighbours and at most one exocyclic substituent, so there is no pair to be cis or
    trans -- and detection by length is unreliable there: aromatic bonds sit right at
    the single/double threshold, so one bond of a benzene ring gets picked and the other
    five do not, which would break the ring's symmetry for nothing."""
    rings = [set(ring['atoms']) for ring in getattr(mol_data, 'rings', {}).values()]
    return [
        (a, b) for (a, b) in double_bonds
        if not any(a in ring and b in ring for ring in rings)
    ]


def cis_pairs(mol_data: MolData, double_bonds: Iterable[Tuple[int, int]],
              log: Optional[Logger] = None) -> List[CisPair]:
    """Every (a, b) with a bonded to A, b bonded to B and a-A=B-b cis, for each double
    bond (A, B). A double bond with a substituent pair that is neither cis nor trans is
    twisted and contributes nothing."""
    atoms = mol_data.atoms
    pairs = []
    for (atom_a, atom_b) in double_bonds:
        substituents_a = [i for i in atoms[atom_a]['conn'] if i != atom_b]
        substituents_b = [i for i in atoms[atom_b]['conn'] if i != atom_a]
        bond_pairs = []
        twisted = False
        for a in substituents_a:
            for b in substituents_b:
                phi = abs(dihedral_angle(atoms[a], atoms[atom_a], atoms[atom_b], atoms[b]))
                if phi < CIS_MAX_DIHEDRAL:
                    bond_pairs.append((a, b))
                elif phi <= TRANS_MIN_DIHEDRAL:
                    twisted = True
        if twisted:
            if log:
                log.debug('Double bond {0}={1} is twisted; its cis/trans relations are ignored'.format(
                    atoms[atom_a]['symbol'], atoms[atom_b]['symbol']))
            continue
        pairs.extend(bond_pairs)
    return pairs


def stereo_sources(mol_data: MolData, classes: Dict[int, int], orientations: Dict[int, Orientation],
                   double_bonds: Iterable[Tuple[int, int]]) -> Set[int]:
    """The atoms that can make constitutionally equivalent atoms diastereotopic, judged on
    the constitutional classes: a tetrahedral atom with four differently-classed
    neighbours; both atoms of a double bond with differently-classed substituents on
    either side; a tetrahedral ring atom whose two neighbours outside some ring through
    it are differently classed (chlorocyclohexane C1, a decalin bridgehead)."""
    atoms = mol_data.atoms
    sources = set()

    for atom_id in orientations:
        neighbour_classes = [classes[n] for n in atoms[atom_id]['conn']]
        if len(set(neighbour_classes)) == 4:
            sources.add(atom_id)

    for (atom_a, atom_b) in double_bonds:
        for (this, other) in ((atom_a, atom_b), (atom_b, atom_a)):
            substituent_classes = [classes[n] for n in atoms[this]['conn'] if n != other]
            if len(substituent_classes) == 2 and substituent_classes[0] != substituent_classes[1]:
                sources.update((atom_a, atom_b))

    for ring in getattr(mol_data, 'rings', {}).values():
        ring_atoms = set(ring['atoms'])
        for atom_id in ring_atoms:
            if atom_id not in orientations:
                continue
            outside = [n for n in atoms[atom_id]['conn'] if n not in ring_atoms]
            if len(outside) == 2 and classes[outside[0]] != classes[outside[1]]:
                sources.add(atom_id)

    return sources


def bond_distances_from(mol_data: MolData, sources: Iterable[int]) -> Dict[int, int]:
    """Bond-count distance from every atom to the nearest source (BFS), with a terminal
    atom (a hydrogen, a halogen, a carbonyl oxygen) given the distance of the atom it
    hangs off: the group is the unit, so a methyl's hydrogens are as far from a
    stereocentre as the methyl carbon is, and cannot be merged back while it stays
    split. Atoms with no source in their connected component are absent."""
    atoms = mol_data.atoms
    distances = {}
    queue = deque()
    for source in sources:
        distances[source] = 0
        queue.append(source)
    while queue:
        atom_id = queue.popleft()
        for neighbour in atoms[atom_id]['conn']:
            if neighbour not in distances:
                distances[neighbour] = distances[atom_id] + 1
                queue.append(neighbour)
    for (atom_id, atom) in atoms.items():
        if len(atom['conn']) == 1 and atom['conn'][0] in distances and atom_id not in sources:
            distances[atom_id] = distances[atom['conn'][0]]
    return distances


def merge_distant_splits(constitutional: Dict[int, int], strict: Dict[int, int],
                         distances: Dict[int, int], max_distance: Optional[int]) -> Tuple[Dict[int, int], List[int]]:
    """The effective classes: strict where an atom is within max_distance bonds of a stereo
    source, constitutional beyond. A source set is invariant under every automorphism, so
    the distance is the same for every atom of a constitutional class and the decision is
    made per class. Returns (classes, constitutional classes that were merged back)."""
    if max_distance is None:
        return dict(strict), []
    effective = {}
    merged = set()
    for (atom_id, constitutional_class) in constitutional.items():
        distance = distances.get(atom_id)
        if distance is not None and distance <= max_distance:
            effective[atom_id] = ('strict', strict[atom_id])
        else:
            effective[atom_id] = ('constitutional', constitutional_class)
            if len({strict[a] for (a, c) in constitutional.items() if c == constitutional_class}) > 1:
                merged.add(constitutional_class)
    return renumber(effective), sorted(merged)


def renumber(classes: Dict[int, object]) -> Dict[int, int]:
    """Classes as consecutive integers, numbered in order of each class's smallest atom
    id -- the order nauty prints orbits in, so an unchanged molecule keeps its numbering."""
    first_atom = {}
    for (atom_id, label) in classes.items():
        first_atom[label] = min(first_atom.get(label, atom_id), atom_id)
    order = {label: n for (n, label) in enumerate(sorted(first_atom, key=first_atom.get))}
    return {atom_id: order[label] for (atom_id, label) in classes.items()}
