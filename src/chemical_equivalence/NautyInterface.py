import subprocess
import tempfile
from typing import Union, Optional, List, Dict, Tuple
from itertools import groupby
from operator import itemgetter

from chemical_equivalence.helpers.types_helpers import Logger
from chemical_equivalence.helpers.atoms import EQUIVALENCE_CLASS_KEY
from chemical_equivalence.helpers.iterables import concat
from chemical_equivalence.config import DREADNAUT_EXECUTABLE

from atb_outputs.helpers.types_helpers import MolData

atb_to_nauty = lambda x: (x - 1)
nauty_to_atb = lambda x: (x + 1)

Partition = Dict[int, List[int]]


def generate_nauty_edges_str(mol_data: MolData) -> str:
    return ''.join(
        [
            "{0}:{1};".format(*[atb_to_nauty(index) for index in bond['atoms']])
            for bond in mol_data.bonds
        ]
    )


def nauty_graph(mol_data: MolData, nauty_node_partition: Optional[Partition] = None) -> str:
    return 'n={num_atoms} g {edges}.f=[{node_partition}]'.format(
        num_atoms=len(mol_data.atoms),
        edges=generate_nauty_edges_str(mol_data),
        node_partition=nauty_partition_str_for(
            get_partition_for(mol_data) if nauty_node_partition is None else nauty_node_partition)
    )


def atom_descriptor_key_for(atom: Dict[str, Union[str, int, float]]) -> str:
    return 'iacm' if 'iacm' in atom else 'type'


def atom_descriptor_for(atom: Dict[str, Union[str, int, float]]) -> str:
    return str(atom[atom_descriptor_key_for(atom)])


def nauty_partition_str_for(partition: Partition) -> str:
    # Format it in dreadnaut's partition format. Ex: "1,2,3|4,5,6"
    return '|'.join(
        [
            ','.join(map(str, indices))
            for (_, indices) in sorted(partition.items())
        ]
    )


def partition_for_chemical_equivalence_dict(chemical_equivalence_dict: Dict[int, int]) -> Partition:
    return {
        key: [atom_id for (atom_id, equivalence_class_id) in group]
        for (key, group) in groupby(
            sorted(
                chemical_equivalence_dict.items(),
                key=itemgetter(1)
            ),
            key=itemgetter(1),
        )
    }


def get_partition_for(mol_data: MolData) -> Partition:
    # A parition is a dictionnary where keys are iacm or element type (ex:12 for C) and values are a list of matching atom indexes. 
    # Ex: {'12': [2, 4, 7, 10, 13, 16], '20': [1, 3, 5, 6, 8, 9, 11, 12, 14, 15, 17, 18]}
    return {
        # Shift atom indexes by one to match dreadnaut's convention (starts at 0)
        # Atoms are sorted by key=atom_descriptor_for for canonical flovouring of the nauty nodes
        group_key: [atb_to_nauty(atom['id']) for atom in group_iterator]
        for (group_key, group_iterator) in groupby(
            sorted(mol_data.atoms.values(), key=atom_descriptor_for),
            key=atom_descriptor_for,
        )
    }


def get_nauty_node_partition(mol_data: MolData) -> str:
    return nauty_partition_str_for(
        get_partition_for(mol_data),
    )


def generate_nauty_input_from_moldata(mol_data: MolData, log: Optional[Logger] = None) -> str:
    input_str = '{nauty_graph} c xo'.format(
        nauty_graph=nauty_graph(mol_data),
    )

    if log:
        log.debug('Nauty input: {0}'.format(input_str))

    return input_str


def generate_nauty_output_from_inputstr(nauty_input_str: str, log: Optional[Logger] = None) -> str:
    nauty_stdout = _run(
        [DREADNAUT_EXECUTABLE],
        nauty_input_str,
        log=log,
    )

    return nauty_stdout


HAS_FOUND_ISOMORPHISM_MSG = "h and h' are identical."


def get_partition_from_nauty_output(nauty_output_str: str) -> list[tuple[int, ...]]:
    assert HAS_FOUND_ISOMORPHISM_MSG in nauty_output_str, nauty_output_str

    return [
        tuple(map(int, field.split('-')))
        for field in nauty_output_str.split(HAS_FOUND_ISOMORPHISM_MSG)[1].strip().split()
    ]


def calcEquivGroups(mol_data: MolData, log: Optional[Logger] = None) -> Union[str, Dict[int, int]]:
    if log:
        log.debug("Running Nauty")

    nauty_stdout = generate_nauty_output_from_inputstr(
        generate_nauty_input_from_moldata(mol_data, log=log),
        log=log,
    )

    if len(nauty_stdout) == 0:
        if log is not None:
            log.warning("calcEquivGroups: dreadnaut produced no output")
        return ""
    else:

        equivalence_for_atom = nauty_equivalence_dict(nauty_stdout)

        for atom_id in mol_data.atoms.keys():
            mol_data.atoms[atom_id][EQUIVALENCE_CLASS_KEY] = equivalence_for_atom[atom_id]

        # check for the case the object was initiated from a molecule3D then remap the originids on the output
        if "Molecule3D_id" in list(mol_data.atoms.values())[0]:
            equivalence_for_atom = {mol_data.atoms[a]['Molecule3D_id']: symmetry_group_id
                                    for a, symmetry_group_id in equivalence_for_atom.items()}
        if log:
            log.debug("Equivalence groups:\n{0}".format(pretty_equivalence_class_dict(mol_data, equivalence_for_atom)))

        return equivalence_for_atom


def nauty_equivalence_dict(nauty_stdout: str) -> Dict[int, int]:
    orbital_data = nauty_stdout.split("seconds")[-1].strip()

    def eval_group_field(group_field: str) -> List[int]:
        if ':' in group_field:
            # Defines a range
            start, stop = map(int, group_field.split(":"))
            return [int(x) for x in range(start, stop + 1)]
        else:
            return [int(group_field)]

    def eval_group_str(group_str: str) -> List[int]:
        if '(' in group_str:
            # Last element is group size '(N)', not needed
            group_fields = group_str.split()[:-1]
        else:
            group_fields = group_str.split()
        return concat([eval_group_field(group_field) for group_field in group_fields])

    try:
        equivalence_groups = [
            eval_group_str(group_str)
            for group_str in orbital_data.split(";")
            if group_str
        ]
    except:
        print(orbital_data)
        raise

    equivalence_for_atom = {
        nauty_to_atb(index): n
        for n, list_of_indices in enumerate(equivalence_groups)
        for index in list_of_indices
    }

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


def _run(args: List[str], stdin: str, log: Optional[Logger] = None) -> str:
    tmp = tempfile.TemporaryFile(buffering=0)
    tmp.write(stdin.encode())
    tmp.seek(0)

    proc = subprocess.Popen(args, stdin=tmp, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()

    tmp.close()
    if stderr and log:
        log.debug(stderr)
    return stdout.strip().decode()


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
                 invert: bool = False) -> str:
    """dreadnaut input (digraph) for the molecule plus stereo gadgets; ``invert`` builds
    the mirror image. Every atom keeps its constitutional colour cell; gadget vertices
    get cells of their own, which sort after the numeric atom descriptors."""
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

    n_vertices = len(out_edges)
    edges = ';'.join(
        '{0}:{1}'.format(vertex, ' '.join(str(target) for target in sorted(out_edges[vertex])))
        for vertex in range(n_vertices) if out_edges[vertex]
    )
    partition = '|'.join(
        ','.join(str(vertex) for vertex in sorted(cells[cell]))
        for cell in sorted(cells)
    )
    return 'n={0} d g {1}. f=[{2}]'.format(n_vertices, edges, partition)


def stereo_orbits(mol_data: MolData, orientations: Dict[int, int], cis_pairs: List[Tuple[int, int]],
                  log: Optional[Logger] = None) -> Tuple[Dict[int, int], bool]:
    """Equivalence classes under every orientation-respecting automorphism, proper or
    improper, plus whether an improper one (a mirror) exists. Two dreadnaut runs: the
    orbits of the proper automorphisms, then the isomorphism onto the mirror image."""
    vertex_of = atom_vertex_map(mol_data)
    atom_of = {vertex: atom_id for (atom_id, vertex) in vertex_of.items()}
    graph = stereo_graph(mol_data, orientations, cis_pairs)
    mirror = stereo_graph(mol_data, orientations, cis_pairs, invert=True)

    proper_output = generate_nauty_output_from_inputstr(graph + ' c xo', log=log)
    if not proper_output:
        raise Exception('dreadnaut produced no output for the stereo graph')
    orbit_of_vertex = nauty_equivalence_dict(proper_output)  # keys are vertex + 1
    orbit = {atom_of[vertex]: orbit_of_vertex[vertex + 1] for vertex in atom_of}

    mirror_output = generate_nauty_output_from_inputstr(graph + ' c x @ ' + mirror + ' x ##', log=log)
    has_mirror = HAS_FOUND_ISOMORPHISM_MSG in mirror_output
    if has_mirror:
        parent = {}

        def find(x):
            while parent.get(x, x) != x:
                x = parent[x]
            return x

        for (vertex, image) in get_partition_from_nauty_output(mirror_output):
            if vertex in atom_of and image in atom_of:
                (a, b) = (find(orbit[atom_of[vertex]]), find(orbit[atom_of[image]]))
                if a != b:
                    parent[a] = b
        orbit = {atom_id: find(label) for (atom_id, label) in orbit.items()}

    if log:
        log.debug('Stereo orbits: {0} proper, mirror image {1}'.format(
            len(set(orbit_of_vertex.values())), 'found' if has_mirror else 'absent'))
    return orbit, has_mirror


if __name__ == '__main__':
    pass
