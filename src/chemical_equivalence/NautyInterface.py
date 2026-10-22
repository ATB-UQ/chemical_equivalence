"""Graph symmetry for chemical equivalence: nauty, in process, through ``atb_nauty``.

Two questions are asked of nauty, and both answers are invariants of the coloured graph,
not artefacts of a nauty version or of the vertex numbering: the orbits of the
automorphism group, and whether the stereo graph is isomorphic to its mirror image.
"""
from typing import Union, Optional, List, Dict, Tuple
from itertools import groupby

import atb_nauty
from atb_nauty import NautyError

from chemical_equivalence.helpers.types_helpers import Logger
from chemical_equivalence.helpers.atoms import EQUIVALENCE_CLASS_KEY

from atb_outputs.helpers.types_helpers import MolData

#: Until 2026-10 nauty ran as the dreadnaut program and failures raised DreadnautError;
#: kept as an alias for callers that still catch it.
DreadnautError = NautyError


atb_to_nauty = lambda x: (x - 1)
nauty_to_atb = lambda x: (x + 1)

Partition = Dict[str, List[int]]


def atom_descriptor_key_for(atom: Dict[str, Union[str, int, float]]) -> str:
    return 'iacm' if 'iacm' in atom else 'type'


def atom_descriptor_for(atom: Dict[str, Union[str, int, float]]) -> str:
    return str(atom[atom_descriptor_key_for(atom)])


def get_partition_for(mol_data: MolData) -> Partition:
    # A parition is a dictionnary where keys are iacm or element type (ex:12 for C) and values are a list of matching atom indexes.
    # Ex: {'12': [2, 4, 7, 10, 13, 16], '20': [1, 3, 5, 6, 8, 9, 11, 12, 14, 15, 17, 18]}
    return {
        # Shift atom indexes by one to match nauty's convention (starts at 0)
        # Atoms are sorted by key=atom_descriptor_for for canonical flovouring of the nauty nodes
        group_key: [atb_to_nauty(atom['id']) for atom in group_iterator]
        for (group_key, group_iterator) in groupby(
            sorted(mol_data.atoms.values(), key=atom_descriptor_for),
            key=atom_descriptor_for,
        )
    }


def molecule_graph(mol_data: MolData) -> atb_nauty.Graph:
    """The molecular graph: vertex = atom id - 1, one colour cell per atom descriptor
    (``iacm``, else ``type``), cells in descriptor-string order. Atom ids must be 1..N."""
    return atb_nauty.Graph(
        len(mol_data.atoms),
        edges=[tuple(atb_to_nauty(index) for index in bond['atoms']) for bond in mol_data.bonds],
        partition=[cell for (_, cell) in sorted(get_partition_for(mol_data).items())],
    )


def orbit_classes(orbits: Tuple[int, ...]) -> Dict[int, int]:
    """``{vertex: class}`` from nauty's orbits (``orbits[v]`` = least vertex of v's orbit),
    classes numbered by ascending least member: the order dreadnaut printed them in."""
    rank = {least: n for (n, least) in enumerate(sorted(set(orbits)))}
    return {vertex: rank[least] for (vertex, least) in enumerate(orbits)}


def calcEquivGroups(mol_data: MolData, log: Optional[Logger] = None) -> Dict[int, int]:
    if log:
        log.debug("Running Nauty")

    graph = molecule_graph(mol_data)
    if log:
        log.debug('Nauty graph: {0} vertices, {1} cells'.format(graph.n, len(graph.partition)))

    equivalence_for_atom = {
        nauty_to_atb(vertex): n
        for (vertex, n) in orbit_classes(atb_nauty.automorphisms(graph).orbits).items()
    }

    for atom_id in mol_data.atoms.keys():
        mol_data.atoms[atom_id][EQUIVALENCE_CLASS_KEY] = equivalence_for_atom[atom_id]

    # check for the case the object was initiated from a molecule3D then remap the originids on the output
    if "Molecule3D_id" in list(mol_data.atoms.values())[0]:
        equivalence_for_atom = {mol_data.atoms[a]['Molecule3D_id']: symmetry_group_id
                                for a, symmetry_group_id in equivalence_for_atom.items()}
    if log:
        log.debug("Equivalence groups:\n{0}".format(pretty_equivalence_class_dict(mol_data, equivalence_for_atom)))

    return equivalence_for_atom


def pretty_equivalence_class_dict(mol_data: MolData, equivalence_dict: Dict[int, int]) -> str:
    return '\n'.join(
        "{equivalence_class}: {atoms}".format(
            equivalence_class=equivalence_class,
            atoms=' '.join([atom['symbol'] for atom in mol_data.atoms.values() if
                            atom[EQUIVALENCE_CLASS_KEY] == equivalence_class])
        )
        for equivalence_class in sorted(set(equivalence_dict.values()))
    )


# ---------------------------------------------------------------------------------------
# Stereo-aware orbits (see stereo.py for the rule being encoded)
#
# The molecular graph is extended with gadgets so that nauty's automorphisms are exactly
# the label permutations that respect every tetrahedral orientation and cis relation:
#
#  * a tetrahedral atom X with neighbours n1<n2<n3<n4 gets three "pairing" vertices
#    P(12|34), P(13|24), P(14|23), each hung off two "pair" vertices that join it to its
#    two atom pairs, and the three P's are joined in a directed 3-cycle whose direction is
#    the orientation. S4 acts on the three pairings with kernel V4; an odd permutation of
#    the neighbours induces a transposition of the P's, which a directed cycle forbids, so
#    exactly the even permutations survive -- orientation preserved;
#  * a cis pair (a, b) across a double bond gets one vertex adjacent to both, so cis pairs
#    can only map onto cis pairs.
#
# Reversing every P-cycle gives the mirror image graph; an isomorphism from the graph to
# its mirror is an improper symmetry (handedness -1 everywhere). The orbits of all valid
# automorphisms are the orbits of the proper ones merged along one such isomorphism.

PAIRINGS = ((0, 1, 2, 3), (0, 2, 1, 3), (0, 3, 1, 2))  # (i, j | k, l) as index pairs


def atom_vertex_map(mol_data: MolData) -> Dict[int, int]:
    """Atom id -> nauty vertex, atoms in increasing id order (so it coincides with
    atb_to_nauty whenever ids are 1..N)."""
    return {atom_id: vertex for (vertex, atom_id) in enumerate(sorted(mol_data.atoms))}


def stereo_graph(mol_data: MolData, orientations: Dict[int, int], cis_pairs: List[Tuple[int, int]],
                 invert: bool = False) -> atb_nauty.Graph:
    """The molecule plus stereo gadgets, as a digraph; ``invert`` builds the mirror image.
    Every atom keeps its constitutional colour cell; gadget vertices get cells of their
    own, which sort after the numeric atom descriptors."""
    vertex_of = atom_vertex_map(mol_data)
    out_edges: Dict[int, set] = {vertex: set() for vertex in vertex_of.values()}
    cells: Dict[str, List[int]] = {}
    for (atom_id, atom) in mol_data.atoms.items():
        cells.setdefault(atom_descriptor_for(atom), []).append(vertex_of[atom_id])
    for bond in mol_data.bonds:
        (a, b) = (vertex_of[bond['atoms'][0]], vertex_of[bond['atoms'][1]])
        out_edges[a].add(b)
        out_edges[b].add(a)

    def new_vertex(cell: str) -> int:
        vertex = len(out_edges)
        out_edges[vertex] = set()
        cells.setdefault(cell, []).append(vertex)
        return vertex

    def undirected(u: int, v: int) -> None:
        out_edges[u].add(v)
        out_edges[v].add(u)

    for (atom_id, orientation) in sorted(orientations.items()):
        neighbours = [vertex_of[n] for n in sorted(mol_data.atoms[atom_id]['conn'])]
        pairing_vertices = []
        for (i, j, k, l) in PAIRINGS:
            pairing = new_vertex('zP')
            for (first, second) in ((i, j), (k, l)):
                pair = new_vertex('zQ')
                undirected(pair, pairing)
                undirected(pair, neighbours[first])
                undirected(pair, neighbours[second])
            pairing_vertices.append(pairing)
        forward = (orientation > 0) != invert
        for n in range(3):
            (u, v) = (pairing_vertices[n], pairing_vertices[(n + 1) % 3])
            out_edges[u].add(v) if forward else out_edges[v].add(u)

    for (a, b) in cis_pairs:
        cis = new_vertex('zCIS')
        undirected(cis, vertex_of[a])
        undirected(cis, vertex_of[b])

    return atb_nauty.Graph(
        len(out_edges),
        adjacency={vertex: sorted(targets) for (vertex, targets) in out_edges.items()},
        directed=True,
        partition=[sorted(cells[cell]) for cell in sorted(cells)],
    )


def stereo_orbits(mol_data: MolData, orientations: Dict[int, int], cis_pairs: List[Tuple[int, int]],
                  log: Optional[Logger] = None) -> Tuple[Dict[int, int], bool]:
    """Equivalence classes under every orientation-respecting automorphism, proper or
    improper, plus whether an improper one (a mirror) exists. Two nauty questions: the
    orbits of the proper automorphisms, then an isomorphism onto the mirror image."""
    vertex_of = atom_vertex_map(mol_data)
    atom_of = {vertex: atom_id for (atom_id, vertex) in vertex_of.items()}
    graph = stereo_graph(mol_data, orientations, cis_pairs)
    mirror = stereo_graph(mol_data, orientations, cis_pairs, invert=True)

    orbit_of_vertex = orbit_classes(atb_nauty.automorphisms(graph).orbits)
    orbit = {atom_of[vertex]: orbit_of_vertex[vertex] for vertex in atom_of}

    to_mirror = atb_nauty.isomorphism(graph, mirror)
    has_mirror = to_mirror is not None
    if has_mirror:
        # Any isomorphism onto the mirror will do: the improper symmetries are one coset
        # of the proper ones, so the merged orbits do not depend on which one nauty finds.
        parent = {}

        def find(x):
            while parent.get(x, x) != x:
                x = parent[x]
            return x

        for vertex in atom_of:
            image = to_mirror[vertex]
            if image in atom_of:
                (a, b) = (find(orbit[atom_of[vertex]]), find(orbit[atom_of[image]]))
                if a != b:
                    parent[a] = b
        orbit = {atom_id: find(label) for (atom_id, label) in orbit.items()}

    if log:
        log.debug('Stereo orbits: {0} proper, mirror image {1}'.format(
            len(set(orbit_of_vertex.values())), 'found' if has_mirror else 'absent'))
    return orbit, has_mirror
